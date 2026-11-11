function kerr_comb_application_narma10()
% =========================================================
% kerr_comb_application_narma10.m
%
% Kerr 微环光梳 + RK4IP 数值解 LLE
% 用于 NARMA10 任务的 Reservoir Computing
%
% 使用的 LLE 形式：
%
%   ∂E/∂t = -(α + iδ0)/tR * E
%           - i (L/tR) (β2/2) ∂^2E/∂τ^2
%           + i (γL/tR) |E|^2 E
%           + sqrt(θ)/tR * Ein
%
% 数值方法：4 阶 Runge-Kutta Interaction Picture (RK4IP)
%
% 任务：
%   - 输入 u(k) ~ U(0, 0.5)
%   - 目标 y(k) 由标准 NARMA10 公式生成
%   - 微环输出频谱功率作为 reservoir states
%   - 岭回归预测 NARMA10 输出
% =========================================================

clc;
clear;
close all;

%% ---------- 1. NARMA10 任务参数 ----------
rng(1);   % 固定随机种子，保证可重复

tau_pred   = 1;      % 一步预测
N_washout  = 800;      % 热身数据
N_train    = 5360;%1360;   % 训练长度
N_test     = 800;%340;    % 测试长度

N_steps = N_washout + N_train + N_test;

% 为了生成 NARMA10，前面多留一些点
N_total = N_steps + tau_pred + 50;

% 标准 NARMA10 输入：u(k) ∈ [0, 0.5]
u = 0.5 * rand(N_total, 1);

% 生成 NARMA10 目标序列
y = zeros(N_total, 1);
for k = 10:N_total-1
    y(k+1) = 0.3 * y(k) ...
           + 0.05 * y(k) * sum(y(k-9:k)) ...
           + 1.5 * u(k-9) * u(k) ...
           + 0.1;
end

% 丢掉前面一段初始过渡，避免开头全零影响
discard = 50;
u = u(discard+1:end);
y = y(discard+1:end);

if length(u) < N_steps + tau_pred
    error('NARMA10 序列长度不足，请增大 N_total。');
end

%% ---------- 2. 输入映射到泵浦功率 ----------
% 对于 NARMA10，u ∈ [0,0.5]
% 线性映射到 [-1,1]
u_drive = 4*u - 1;

u_drive_used = u_drive(1:N_steps);

% 泵浦参数
Pav   = 0.325;   % W，平均泵浦功率
Delta = 0.05;    % W，变化幅度

% Pin ∈ [Pav-Delta/2, Pav+Delta/2]
Pin_series = Pav + (Delta/2) * u_drive_used;
Pin_series = max(Pin_series, 0);   % 防止负功率
Ein_series = sqrt(Pin_series);     % Pin = |Ein|^2

%% ---------- 3. 微环参数 & LLE 网格 ----------
c0 = 2.998e8;          % 光速 (m/s)
nn = 2^8;              % 256
N  = 2 * nn;           % 512

FSR = 200e9;           % Hz
t_R = 1 / FSR;         % s

dtau = t_R / N;
tau  = (-nn:nn-1) * dtau;

D1 = 2*pi*FSR;
k  = -nn:nn-1;
OM = D1 * k;           % rad/s

lam_center = 1.600e-6;     % m
f_center   = c0 / lam_center;
f_axis     = f_center + OM/(2*pi);

% super-Gaussian 滤波器
f_width      = 64e12;%4e12;
filter_order = 4;
f_half       = f_width / 2;
H_filter = exp(-((f_axis - f_center)/f_half).^(2*filter_order));

% LLE 参数
Beta_2 = -4.0e-9 * 1e-18;   % ns^2/m -> s^2/m
gamma  = 1.6e-2;            % W^-1 m^-1
delta0_target =4e-3;

neff = 1.444;
L    = t_R * c0 / neff;

Q_int = 1.0e7;
Qext1 = 1.0e7;
Qext2 = 1.0e7;
Q_tot = 1 / (1/Q_int + 1/Qext1 + 1/Qext2);

nu0  = f_center;
t_ph = Q_tot / (2*pi*nu0);

alpha = t_R / (2*t_ph);

t_e1  = Qext1 / (2*pi*nu0);
theta = t_R / t_e1;

fprintf('==== 物理参数诊断 ====\n');
fprintf('Q_int = Q_ext1 = Q_ext2 = %.2e, Q_tot = %.2e\n', Q_int, Q_tot);
fprintf('FSR        = %.3f GHz\n', FSR*1e-9);
fprintf('t_R        = %.3f ps\n', t_R*1e12);
fprintf('t_ph       = %.3f ns\n', t_ph*1e9);
fprintf('alpha      = %.3e\n', alpha);
fprintf('theta      = %.3e\n', theta);
fprintf('Ring length L = %.3f mm\n\n', L*1e3);

symbol_rate = 1.0e9;
T_sym       = 1 / symbol_rate;

N_substeps_sym = 1000;
dt_slow = T_sym / N_substeps_sym;

fprintf('Data rate  = %.1f GSa/s, T_sym = %.2f ns\n', symbol_rate/1e9, T_sym*1e9);
fprintf('dt_slow    = %.3f ps, 每符号 %d 步（RK4IP）\n\n', dt_slow*1e12, N_substeps_sym);

params.t_R   = t_R;
params.L     = L;
params.Beta2 = Beta_2;
params.gamma = gamma;
params.alpha = alpha;
params.theta = theta;

%% ---------- 4. 初始腔内场 ----------
rng(1);
E = 1e-6 * (randn(1,N) + 1i*randn(1,N));

%% ---------- 5. 阶段1：扫描失谐到工作点 ----------
N_scan_sym  = 320;
delta_start = -2e-2;
Ein_bias    = sqrt(Pav);

fprintf('阶段1：扫描失谐从 %.4g 到 %.4g，历时 %d 个符号。\n', ...
    delta_start, delta0_target, N_scan_sym);

for n = 1:N_scan_sym
    frac   = (n-1) / max(N_scan_sym-1, 1);
    delta0 = delta_start + (delta0_target - delta_start) * frac;

    for kk = 1:N_substeps_sym
        E = rk4ip_lle_step(E, Ein_bias, delta0, params, dt_slow, OM);
    end
end

fprintf('阶段1结束：当前失谐 delta0 = %.4g。\n\n', delta0_target);


fprintf('阶段1结束：当前失谐 delta0 = %.4g。\n\n', delta0_target);



E_after_scan = E;
F_after_scan = fftshift(fft(E_after_scan));
I_t_scan = abs(E_after_scan).^2;
I_f_scan = abs(F_after_scan).^2;

%% ---------- 6. 阶段2：固定失谐 + NARMA10 输入驱动 ----------
X = zeros(N_steps, N);

fprintf('阶段2：NARMA10 输入驱动 + LLE (RK4IP) 演化，共 %d 个符号。\n', N_steps);

for n = 1:N_steps
    Ein_n = Ein_series(n);

    for kk = 1:N_substeps_sym
        E = rk4ip_lle_step(E, Ein_n, delta0_target, params, dt_slow, OM);
    end

    F_spec = fftshift(fft(E));
    F_filt =F_spec ;   % 如果你想启用滤波，改成：F_spec .* H_filter
    X(n,:) = abs(F_filt).^2;
end

fprintf('阶段2结束：已获得 N_steps = %d 个功率谱样本。\n\n', N_steps);

%% ---------- 7. 构造训练/测试数据 ----------
X_eff = X(N_washout+1 : N_washout+N_train+N_test, :);

% 目标改成 NARMA10 输出
y_eff = y(N_washout+1+tau_pred : N_washout+N_train+N_test+tau_pred);

X_train = X_eff(1:N_train, :);
y_train = y_eff(1:N_train);

X_test  = X_eff(N_train+1 : N_train+N_test, :);
y_test  = y_eff(N_train+1 : N_train+N_test);

%% ---------- 8. 特征标准化 ----------
[X_train_n, mu_X, sigma_X] = zscore_features(X_train);
X_test_n = (X_test - mu_X) ./ sigma_X;

Xtr = [X_train_n, ones(N_train,1)];
Xte = [X_test_n,  ones(N_test,1)];

%% ---------- 9. λ 扫描，选择最优岭回归 ----------
lambdas   = logspace(-9, 2, 200);
nmse_scan = zeros(size(lambdas));

for ii = 1:length(lambdas)
    lam = lambdas(ii);
    Ireg = eye(size(Xtr,2));
    Ireg(end,end) = 0;   % bias 不正则化

    beta_tmp = (Xtr.'*Xtr + lam*Ireg) \ (Xtr.'*y_train);
    y_test_pred_tmp = Xtr * beta_tmp;

    nmse_scan(ii) = mean((y_test_pred_tmp - y_train).^2) / var(y_train, 1);
end

[best_nmse, best_idx] = min(nmse_scan);
lambda_opt = lambdas(best_idx);

fprintf('lambda 扫描: 最优 lambda = %.2e, 对应测试 NMSE = %.3e\n', ...
    lambda_opt, best_nmse);

figure('Color','w','Name','Lambda vs NMSE (Ridge Scan)');
semilogx(lambdas, nmse_scan, '-o', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
hold on;
semilogx(lambda_opt, best_nmse, 's', 'MarkerSize', 8, 'LineWidth', 1.5);
yl = ylim;
plot([lambda_opt lambda_opt], yl, 'k--', 'LineWidth', 1.0);
ylim(yl);
xlabel('\lambda');
ylabel('Test NMSE');
title('Ridge: Test NMSE vs \lambda');
legend('Scan NMSE','Best point','Best \lambda','Location','best');

% 用最优 lambda 重训
Ireg = eye(size(Xtr,2));
Ireg(end,end) = 0;
beta = (Xtr.'*Xtr + lambda_opt*Ireg) \ (Xtr.'*y_train);

y_train_pred = Xtr * beta;
y_test_pred  = Xte * beta;

nmse_train = mean((y_train_pred - y_train).^2) / var(y_train,1);
nmse_test  = mean((y_test_pred  - y_test ).^2) / var(y_test,1);

fprintf('最终岭回归结果：\n');
fprintf('  训练 NMSE = %.3e\n', nmse_train);
fprintf('  测试  NMSE = %.3e\n', nmse_test);

%% ---------- 10. 画图：阶段1结束时域 & 频域 ----------
lam_axis = c0 ./ f_axis * 1e9;

figure('Name','阶段1结束：时域和频域');
subplot(2,1,1);
plot(tau*1e12, 10*log10(I_t_scan + 1e-12), 'b', 'LineWidth', 1);
xlabel('Fast time \tau (ps)');
ylabel('Intensity (dB)');
title('阶段1结束：腔内时域强度');
grid on;

subplot(2,1,2);
plot(lam_axis, 10*log10(I_f_scan + 1e-12), 'r', 'LineWidth', 1);
xlabel('Wavelength (nm)');
ylabel('Power (dB)');
title('阶段1结束：频域光梳');
xlim([1500,1700]);
grid on;

%% ---------- 11. 画图：NARMA10 测试段真实值 vs 预测值 ----------
figure('Name','NARMA10 一步预测：测试段');
plot(y_test, 'k', 'LineWidth', 1);
hold on;
plot(y_test_pred, 'r--', 'LineWidth', 1);
xlabel('Time step (test)');
ylabel('NARMA10 output');
legend('Ground truth', 'Prediction');
title(sprintf('NARMA10 一步预测, NMSE = %.3e', nmse_test));
grid on;

%% ---------- 12. 画图：功率谱矩阵 ----------
figure('Name','功率谱矩阵');
imagesc(10*log10(X_eff.' + 1e-12));
axis xy;
xlabel('Time index (after washout)');
ylabel('Spectral channel index');
title('Microcomb Reservoir: Power Spectral Densities (dB)');
colorbar;

%% ---------- 13. 画图：NARMA10 输入和目标 ----------
figure('Name','NARMA10 输入与目标');
subplot(2,1,1);
plot(u(1:N_steps), 'b', 'LineWidth', 1);
xlabel('Time step');
ylabel('u(k)');
title('NARMA10 输入序列');
grid on;

subplot(2,1,2);
plot(y(1:N_steps+tau_pred), 'k', 'LineWidth', 1);
xlabel('Time step');
ylabel('y(k)');
title('NARMA10 目标序列');
grid on;

end


%% 子函数：单步 RK4IP 求解 LLE
function E_next = rk4ip_lle_step(E, Ein, delta0, params, dt, OM)

t_R   = params.t_R;
L     = params.L;
Beta2 = params.Beta2;
gamma = params.gamma;
alpha = params.alpha;
theta = params.theta;

D_lin = (-(alpha + 1i*delta0)/t_R) + ...
        1i*(L/t_R)*(Beta2/2) .* (OM.^2);

expD_h2      = exp(D_lin * (dt/2));
expD_h       = exp(D_lin * dt);
expMinusD_h2 = exp(-D_lin * (dt/2));
expMinusD_h  = exp(-D_lin * dt);

c_nl  = 1i * gamma * L / t_R;
drive = (sqrt(theta)/t_R) * Ein;

apply_expD = @(X,expD) ifft(ifftshift(fftshift(fft(X)) .* expD));

U0 = E;

% K1
E0 = U0;
N1 = c_nl .* abs(E0).^2 .* E0 + drive;
K1 = N1;

% K2
Utemp = U0 + (dt/2) * K1;
E2    = apply_expD(Utemp, expD_h2);
N2    = c_nl .* abs(E2).^2 .* E2 + drive;
K2    = apply_expD(N2, expMinusD_h2);

% K3
Utemp = U0 + (dt/2) * K2;
E3    = apply_expD(Utemp, expD_h2);
N3    = c_nl .* abs(E3).^2 .* E3 + drive;
K3    = apply_expD(N3, expMinusD_h2);

% K4
Utemp = U0 + dt * K3;
E4    = apply_expD(Utemp, expD_h);
N4    = c_nl .* abs(E4).^2 .* E4 + drive;
K4    = apply_expD(N4, expMinusD_h);

U1 = U0 + (dt/6) * (K1 + 2*K2 + 2*K3 + K4);
E_next = apply_expD(U1, expD_h);

end


%% 子函数：逐列 z-score 标准化
function [Xn, mu_X, sigma_X] = zscore_features(X)
mu_X    = mean(X, 1);
sigma_X = std(X, 0, 1);
sigma_X(sigma_X < 1e-12) = 1e-12;
Xn = (X - mu_X) ./ sigma_X;
end