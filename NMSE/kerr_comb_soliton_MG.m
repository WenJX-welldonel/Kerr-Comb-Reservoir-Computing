function kerr_comb_application()
% =========================================================
% kerr_comb_application.m
%
% Kerr 微环光梳 + RK4IP 数值解 LLE，用于 Mackey-Glass(MG)
% 序列的一步预测（Reservoir Computing）。
%
% 使用的是论文中的 LLE 形式：
%
%   ∂E/∂t = -(α + iδ0)/tR * E
%           - i (L/tR) (β2/2) ∂^2E/∂τ^2
%           + i (γL/tR) |E|^2 E
%           + sqrt(θ)/tR* Ein
%
% 数值求解方法：4 阶 Runge-Kutta Interaction Picture (RK4IP)。
%
% 关键物理参数（MG 任务默认，对齐论文）：
%   Q_int   = 1.0e7
%   Q_ext1  = 1.0e7
%   Q_ext2  = 1.0e7
%   Q_tot   = (1/Q_int + 1/Q_ext1 + 1/Q_ext2)^(-1) ≈ 3.33e6
%   γ       = 1.6e-2  W^-1 m^-1
%   β2      = -4.0e-9 ns^2 / m  (→ -4.0e-27 s^2/m)
%   δ0      = 1.6e-3
%   Pav     = 0.325 W
%   Δ       = 0.05 W
%   data rate = 1 GSa/s
%
% 读出层：微环 drop port 光谱功率，经 1 THz 四阶 super-Gaussian 滤波，
% 在线性回归前做特征标准化，使用岭回归预测 MG(n+1)。
%
% =========================================================

clc; close all;

%% ---------- 1. 读取 Mackey-Glass 序列 ----------
data = load('MG.mat');
if ~isfield(data, 'MG')  %判断结构体 data 里有没有字段 MG
    error('MG.mat 中未找到变量 MG，请确认变量名 MG。');
end
mg = data.MG(:);               % 列向量
N_total = length(mg);

% 预测任务：一步预测 MG(n+1)
tau_pred   = 1;    %单步
N_washout  = 800;  %热身数据
N_train    = 5360; %训练数据
N_test     = 800;  %测试数据

N_steps = N_washout + N_train + N_test;  % 总的用到的符号数，2800也就意味着还有301个用于预测

if N_steps + tau_pred > N_total
    error('MG 序列长度不足，请减小 N_train/N_test。');
end

%% ---------- 2. MG 归一化 → 映射到泵浦功率 Pin ∈ [Pav±Δ/2] ----------
u = mg;
u_min  = min(u);
u_max  = max(u);
u_norm = (u - u_min) / (u_max - u_min + eps);  % [0,1]    eps 是 MATLAB 的机器精度级别的极小正数，为的是防止分母为0
u_drive = 2*u_norm - 1;                        % [-1,1]
u_drive_used = u_drive(1:N_steps);             % 真正用于驱动的部分，取前2800个

% 泵浦参数
Pav   = 0.325;   % W, 平均泵浦功率
Delta = 0.05;    % W, 变化幅度

% 线性映射到功率：Pin ∈ [Pav - Δ/2, Pav + Δ/2] = [0.30, 0.35] W
Pin_series = Pav + (Delta/2) * u_drive_used;
Pin_series = max(Pin_series, 0);    % 防止出现负功率，在我们这个任务中其实不太可能，因为归一化到-1和1
Ein_series = sqrt(Pin_series);      % 注入场幅：这里的E是定义为单位为sqr（w）的E  Pin = |Ein|^2

%% ---------- 3. 微环参数 & LLE 网格 ----------
c0 = 2.998e8;           % 光速 (m/s)
nn = 2^8;               % 一半采样点数 256
N  = 2*nn;              % 总采样点（512）

% 自由光谱范围 FSR & roundtrip time
FSR = 200e9;            % Hz
t_R = 1/FSR;            % s, 单圈时间 ≈ 5 ps

% fast-time 网格（沿腔一圈）
dtau = t_R / N;                             % fast-time 步长
tau  = (-nn:nn-1) * dtau;                   % fast-time 轴（只用于画图）

% 频率网格：间隔 = FSR（rad/s）
D1  = 2*pi*FSR;                             % rad/s 计算角频率
k   = -nn:nn-1;
OM  = D1 * k;                               % rad/s，对应各梳齿频率偏移

% 中心波长 / 频率，用来构造 1 THz super-Gaussian 滤波器
lam_center = 1.600e-6;                      % m (~1600 nm)
f_center   = c0 / lam_center;               % Hz
f_axis     = f_center + OM/(2*pi);          % Hz，实际频率轴

% super-Gaussian 滤波器：阶数 4，宽度 1 THz
f_width      = 4e12;                        % 1 THz
filter_order = 4;
f_half       = f_width/2;
H_filter = exp( - ((f_axis - f_center)/f_half).^(2*filter_order) );  % 1×N 512 只有11个数有值

% LLE 物理参数
Beta_2 = -4.0e-9 * 1e-18;    % ns^2/m → s^2/m 正确的单位换算
gamma  = 1.6e-2;             % W^-1 m^-1
delta0_target = 4e-3;%1.6e-3;      % 目标失谐 δ0

% 环长：假设有效折射率 neff ≈ 1.44这是一个经验值，论文中并为给出
neff = 1.444;
L    = t_R * c0 / neff;      % m

% Q-factor → 光子寿命 → α, θ（对齐论文：Q_int = Q_ext1 = Q_ext2 = 1e7）
Q_int  = 1.0e7;
Qext1  = 1.0e7;
Qext2  = 1.0e7;
Q_tot  = 1 / (1/Q_int + 1/Qext1 + 1/Qext2);    % loaded Q ≈ 3.33e6

nu0   = f_center;                     % 载频（中心频率）f_center   = c0 / lam_center;
t_ph  = Q_tot / (2*pi*nu0);           % photon lifetime（光子寿命） τp = Q_tot/ω0

alpha = t_R / (2*t_ph);              % 振幅损耗系数（无量纲）

% 输入耦合对应的时间常数 → θ
t_e1  = Qext1 / (2*pi*nu0);
theta = t_R / t_e1;                  % 耦合系数（无量纲）

fprintf('==== 物理参数诊断 ====\n');
fprintf('Q_int = Q_ext1 = Q_ext2 = %.2e, Q_tot = %.2e\n', Q_int, Q_tot);
fprintf('FSR        = %.3f GHz\n', FSR*1e-9);
fprintf('t_R        = %.3f ps\n', t_R*1e12);
fprintf('t_ph       = %.3f ns (Q_tot = %.2e)\n', t_ph*1e9, Q_tot);
fprintf('alpha      = %.3e\n', alpha);
fprintf('theta      = %.3e\n', theta);
fprintf('Ring length L = %.3f mm\n\n', L*1e3);

% 数据率：1 GSa/s → 每个符号持续 1 ns
symbol_rate = 1.0e9;                  % 1G Sa/s
T_sym       = 1 / symbol_rate;        % s 一个符号的持续时间  1ns

% 选一个固定的 LLE 时间步（slow time），例如每个符号均匀切分为 20 步
N_substeps_sym =1000;
dt_slow = T_sym / N_substeps_sym; %每个符号只积分 1 次（步长 = 1 ns），步长通常太大，数值会不准甚至不稳定。所以把一个符号拆成 N_substeps_sym 段，用更小的时间步 dt_slow 去积分。

fprintf('Data rate  = %.1f GSa/s, T_sym = %.2f ns\n', symbol_rate/1e9, T_sym*1e9);
fprintf('dt_slow    = %.3f ps, 每符号 %d 步（RK4IP）\n\n', dt_slow*1e12, N_substeps_sym);

% RK4IP 需要的一些参数打包
params.t_R   = t_R;
params.L     = L;
params.Beta2 = Beta_2;
params.gamma = gamma;
params.alpha = alpha;
params.theta = theta;

%% ---------- 4. 初始腔内场（小噪声） ----------
rng(1);  % 为可重复性固定随机种子
E = 1e-6 * (randn(1,N) + 1i*randn(1,N));   % 初始为弱噪声

%% ---------- 5. 阶段1：扫描失谐到工作点 δ0_target ----------
N_scan_sym  = 600;         % 扫描时长（单位：符号数，可自行调整）1个符号是1ns
delta_start = -2e-2;          % 从 0 扫到 δ0_target
Ein_bias    = sqrt(Pav);    % 扫描阶段用平均功率 Pav

fprintf('阶段1：扫描失谐从 %.4g 到 %.4g，历时 %d 个符号。\n',...
    delta_start, delta0_target, N_scan_sym);

for n = 1:N_scan_sym
    % 当前符号的失谐（线性扫描）
    frac   = (n-1) / max(N_scan_sym-1,1);%max(N_scan_sym-1, 1) 的作用，这是一个“防止除零”的保护：为的是代码可以调整
    delta0 = delta_start + (delta0_target - delta_start) * frac;

    Ein_n = Ein_bias;     %Ein_bias    = sqrt(Pav);   扫描阶段用平均功率 Pav

    % 在一个符号周期内，Ein_n 和 delta0 保持不变，做 N_substeps_sym 个 RK4IP 步
    for kk = 1: N_substeps_sym
        E = rk4ip_lle_step(E, Ein_n, delta0, params, dt_slow, OM);
       
    end
end

fprintf('阶段1结束：当前失谐 delta0 = %.4g。\n\n', delta0_target);

N_hold_sym = 400;
fprintf('开始 hold：固定失谐 delta0 = %.4g，固定泵浦 Ein_bias，持续 %d 个符号。\n', ...
    delta0_target, N_hold_sym);

for n = 1:N_hold_sym
    for kk = 1:N_substeps_sym
        E = rk4ip_lle_step(E, Ein_bias, delta0_target, params, dt_slow, OM);
    end
end

fprintf('hold 结束。\n\n');


% 保存扫描结束时的场分布        图上可以看出来不够混沌
E_after_scan = E;
F_after_scan = fftshift(fft(E_after_scan));
E_after_scan = fftshift(E_after_scan);
I_t_scan = abs(E_after_scan).^2;
I_f_scan = abs(F_after_scan).^2;







%% ---------- 6. 阶段2：固定 δ0 + MG 调制泵浦，采集功率谱矩阵 X ----------
X = zeros(N_steps, N);   % 每行 = 一个符号对应的功率谱

fprintf('阶段2：MG 驱动 + LLE (RK4IP) 演化，共 %d 个符号。\n', N_steps);

for n = 1:N_steps
    Ein_n = Ein_series(n);     % 本符号对应的注入泵浦场幅

    % 一个符号周期内，Ein_n 固定，演化 N_substeps_sym 步
    for kk = 1:N_substeps_sym          
       E = rk4ip_lle_step(E, Ein_n, delta0_target, params, dt_slow, OM);
   
    end

    % 取一次输出频谱，经 1 THz super-Gaussian 滤波后求功率
   F_spec = fftshift(fft(E));
    F_filt = F_spec;%.* H_filter;
    X(n,:) = abs( F_filt).^2;
end

fprintf('阶段2结束：已获得 N_steps = %d 个功率谱样本。\n\n', N_steps);

%% ---------- 7. 构造训练/测试数据 ----------
% 舍弃前 N_washout 的状态，只保留后面 N_train+N_test
X_eff = X(N_washout+1 : N_washout+N_train+N_test, :); %x来自X

% 目标是原始 MG(n+1)
y_eff = mg(N_washout+1+tau_pred : ...
           N_washout+N_train+N_test+tau_pred);  %y来自mg         215进行了修改mg改成Pin_series

% 划分训练/测试
X_train = X_eff(1:N_train, :); %取第 1 行到第 N_train 行
y_train = y_eff(1:N_train);  

X_test  = X_eff(N_train+1 : N_train+N_test, :);
y_test  = y_eff(N_train+1 : N_train+N_test);


%% ---------- 8. 特征标准化 ----------?
[X_train_n, mu_X, sigma_X] = zscore_features(X_train);
X_test_n = (X_test - mu_X) ./ sigma_X;

Xtr = [X_train_n, ones(N_train,1)];  % 加偏置
Xte = [X_test_n,  ones(N_test,1)];

%% ---------- 9. λ 网格扫描，选择最优岭回归正则 ----------
lambdas   = logspace(-9, 2, 200);   % 1e-6 ~ e0  7个数  加个图吧    要改        
nmse_scan = zeros(size(lambdas));

for ii = 1:length(lambdas)
    lam = lambdas(ii);
    Ireg = eye(size(Xtr,2));             %eye(p) 生成PXP单位矩阵 I。size(Xtr,2) 是 Xtr 的列数
    Ireg(end,end) = 0;      % 不正则化最后一维（bias）
    beta_tmp = (Xtr.'*Xtr + lam*Ireg) \ (Xtr.'*y_train);
    y_test_pred_tmp = Xtr * beta_tmp;      %预测值
    nmse_scan(ii) = mean((y_test_pred_tmp - y_train).^2) / var(y_train,1);
end

[best_nmse, best_idx] = min(nmse_scan);       % min 返回最小值（best_nmse） ，及其索引（best_idx）：
lambda_opt = lambdas(best_idx);               %lambda_opt最小的lambdas
fprintf('lambda 扫描: 最优 lambda = %.2e, 对应测试 NMSE = %.3e\n',...
        lambda_opt, best_nmse);

 % ==================== 图：lambda vs NMSE（扫描曲线） ====================
figure('Color','w','Name','Lambda vs NMSE (Ridge Scan)');
semilogx(lambdas, nmse_scan, '-o', 'LineWidth', 1.2, 'MarkerSize', 6); 
grid on; hold on;  %hold on是因为下面还要调用这个semilogx函数，避免清空

% 标出最优点（可选：用不同标记强调）
semilogx(lambda_opt, best_nmse, 's', 'MarkerSize', 8, 'LineWidth', 1.5);

% 竖虚线标出最优 lambda
yl = ylim;       %ylim 返回当前 y 轴显示范围，例如 [ymin ymax]
plot([lambda_opt lambda_opt], yl, 'k--', 'LineWidth', 1.0);  %x 坐标是 [lambda_opt lambda_opt]，两个点的 x 相同，所以是一条竖线；
ylim(yl);  %避免画竖线导致 y 轴缩放变化。单纯是因为画图后matlab可能会自己调整

xlabel('\lambda');
ylabel('Test NMSE');
title('Ridge: Test NMSE vs \lambda');  %\lambda 会显示为希腊字母
legend('Scan NMSE','Best point','Best \lambda','Location','best');  %解释图中不同曲线/标记/线段分别代表什么
% =======================================================================   
    
% 用最优 lambda 重新训练
Ireg = eye(size(Xtr,2));
 Ireg(end,end) = 0; 
beta = (Xtr.'*Xtr + lambda_opt*Ireg) \ (Xtr.'*y_train);

y_train_pred = Xtr * beta;
y_test_pred  = Xte * beta;

% NMSE
nmse_train = mean((y_train_pred - y_train).^2) / var(y_train,1);
nmse_test  = mean((y_test_pred  - y_test ).^2) / var(y_test,1);      

fprintf('最终岭回归结果：\n');
fprintf('  训练 NMSE = %.3e\n', nmse_train);
fprintf('  测试  NMSE = %.3e\n', nmse_test);

%% ---------- 10. 画图：扫描结束时的时域 & 频域梳谱 ----------
lam_axis = c0 ./ f_axis * 1e9;    % nm

figure('Name','阶段1结束：时域和频域');
subplot(2,1,1);
plot(tau*1e12,  10*log10(I_t_scan + 1e-12) , 'b','LineWidth',1);
xlabel('Fast time \tau (ps)');
ylabel('Intensity');
title('阶段1结束：腔内时域强度');
grid on;

subplot(2,1,2);
plot(lam_axis, 10*log10(I_f_scan + 1e-12) ,'r','LineWidth',1);
xlabel('Wavelength (nm)');
ylabel('Power (dB)');
title('阶段1结束：频域光梳');
xlim([1500,1700]);
grid on;

%% ---------- 11. 画图：MG 测试段的真实 vs 预测（Baseline） ----------
figure('Name','MG 一步预测：测试段 (Baseline)');
plot(y_test, 'k','LineWidth',1); hold on;
plot(y_test_pred, 'r--','LineWidth',1);
xlabel('Time step (test)');
ylabel('Mackey-Glass');
legend('Ground truth','Prediction');
title(sprintf('Baseline - MG 一步预测, NMSE=%.3e', nmse_test));
grid on;

%% ---------- 12. 画图：功率谱矩阵（时-频演化） ----------
figure('Name','功率谱矩阵');
imagesc(10*log10(X_eff.' + 1e-12));
axis xy;
xlabel('Time index (after washout)');
ylabel('Spectral channel index');
title('Microcomb Reservoir: Filtered Power Spectral Densities (dB)');
colorbar;
end



%% 子函数：单步 RK4IP 求解 LLE
function E_next = rk4ip_lle_step(E, Ein, delta0, params, dt, OM)

t_R   = params.t_R;
L     = params.L;
Beta2 = params.Beta2;
gamma = params.gamma;
alpha = params.alpha;
theta = params.theta;

% 频域线性算符 L_op (作用在 E 上)：频域对角元素 D(Ω)
% D(Ω) = -(α + iδ0)/t_R + i(L/t_R)*(β2/2)*Ω^2
D_lin = (-(alpha + 1i*delta0)/t_R) + ...
        1i*(L/t_R)*(Beta2/2) .* (OM.^2);

% 预计算对角指数 e^{D dt}、e^{±D dt/2}
expD_h2      = exp(D_lin * (dt/2));
expD_h       = exp(D_lin * dt);
expMinusD_h2 = exp(-D_lin * (dt/2));
expMinusD_h  = exp(-D_lin * dt);

% 非线性 + 注入部分 N(E)：
% N(E) = i (γL/t_R) |E|^2 E + sqrt(θ/t_R) Ein
c_nl  = 1i * gamma * L / t_R;
drive = (sqrt(theta)/t_R) * Ein;

% 定义一个工具函数：在频域乘以 expD 再回到时域
apply_expD = @(X,expD) ifft( ifftshift( fftshift(fft(X)) .* expD ) );

% 交互表象变量 U(s)，初值 U(0) = E(0)
U0 = E;
% ---------- K1 ----------
E0 = U0;                          % s=0 时，E = U
N1 = c_nl .* abs(E0).^2 .* E0 + drive;
K1 = N1;                          % tildeN(U0,0) = N1
% ---------- K2 ----------
Utemp = U0 + (dt/2) * K1;
E2    = apply_expD(Utemp, expD_h2);         % 先乘 e^{D dt/2}
N2    = c_nl .* abs(E2).^2 .* E2 + drive;
K2    = apply_expD(N2, expMinusD_h2);       % 乘 e^{-D dt/2}
% ---------- K3 ----------
Utemp = U0 + (dt/2) * K2;
E3    = apply_expD(Utemp, expD_h2);
N3    = c_nl .* abs(E3).^2 .* E3 + drive;
K3    = apply_expD(N3, expMinusD_h2);
% ---------- K4 ----------
Utemp = U0 + dt * K3;
E4    = apply_expD(Utemp, expD_h);
N4    = c_nl .* abs(E4).^2 .* E4 + drive;
K4    = apply_expD(N4, expMinusD_h);

% ---------- 更新 U，再回到物理场 E ----------
U1 = U0 + (dt/6) * (K1 + 2*K2 + 2*K3 + K4);
E_next = apply_expD(U1, expD_h);

end

%% 子函数：特征标准化（逐列 z-score）
function [Xn, mu_X, sigma_X] = zscore_features(X)
mu_X    = mean(X, 1);        %mean(X,1) 表示沿着第 1 维（行方向）求均值，得到每一列的均值。
sigma_X = std(X, 0, 1);      %std(X,0,1) 表示对每一列计算标准差。
sigma_X(sigma_X < 1e-12) = 1e-12;  % 防止除 0
Xn = (X - mu_X) ./ sigma_X;
end 

