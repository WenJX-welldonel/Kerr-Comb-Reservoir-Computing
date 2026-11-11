function kerr_comb_application_narma10_soliton_ftmle()

% =========================================================
% Kerr microcomb + NARMA10 + FTMLE
%
% Operating regime:
%   Soliton
%
% Preparation:
%   600-symbol detuning scan
%   400-symbol hold
%
% Data stage:
%   800-symbol washout
%   5360-symbol training
%   800-symbol testing
%
% Total Stage-2:
%   6960 symbols = 6.960 us at 1 GSa/s
%
% Cache strategy:
%
%   1. Try complete Stage-2 cache
%   2. Try old complete result MAT
%   3. Try prepared soliton cache
%   4. Run simulation only if necessary
%
% FTMLE starts immediately after the hold stage.
%
% Only Figure 6 (FTMLE curve) is enabled.
% All other figures are commented out.
% =========================================================

clc;
clear;
close all;

%% =========================================================
% 0. CACHE CONTROLS
% =========================================================

force_reprepare = false;

force_recompute_stage2 = false;

%% =========================================================
% 1. NARMA10 TASK PARAMETERS
% =========================================================

input_seed = 1;

rng(input_seed,'twister');

tau_pred = 1;

N_washout = 800;
N_train = 5360;
N_test = 800;

N_steps = N_washout + N_train + N_test;

N_total = N_steps + tau_pred + 50;

u = 0.5*rand(N_total,1);

y = zeros(N_total,1);

for k = 10:N_total-1

    term1 = 0.3*y(k);

    term2 = 0.05*y(k)*sum(y(k-9:k));

    term3 = 1.5*u(k-9)*u(k);

    y(k+1) = term1 + term2 + term3 + 0.1;

end

discard = 50;

u = u(discard+1:end);

y = y(discard+1:end);

if length(u) < N_steps + tau_pred
    error('NARMA10 sequence is too short.');
end

%% =========================================================
% 2. INPUT -> PUMP POWER
% =========================================================

u_drive = 4*u - 1;

u_drive_used = u_drive(1:N_steps);

Pav = 0.325;

Delta = 0.05;

Pin_series = Pav + (Delta/2)*u_drive_used;

Pin_series = max(Pin_series,0);

Ein_series = sqrt(Pin_series);

%% =========================================================
% 3. KERR MRR / LLE PARAMETERS
% =========================================================

c0 = 2.998e8;

nn = 2^8;

N = 2*nn;

FSR = 200e9;

t_R = 1/FSR;

dtau = t_R/N;

tau = (-nn:nn-1)*dtau;

D1 = 2*pi*FSR;

k_mode = -nn:nn-1;

OM = D1*k_mode;

lam_center = 1.600e-6;

f_center = c0/lam_center;

f_axis = f_center + OM/(2*pi);

%% =========================================================
% 3.1 SUPER-GAUSSIAN FILTER
% =========================================================

f_width = 4e12;

filter_order = 4;

f_half = f_width/2;

H_filter = exp(-((f_axis-f_center)/f_half).^(2*filter_order));

%% =========================================================
% 3.2 LLE PARAMETERS
% =========================================================

Beta_2 = -4.0e-9*1e-18;

gamma = 1.6e-2;

delta0_target = 4e-3;

neff = 1.444;

L = t_R*c0/neff;

Q_int = 1.0e7;

Qext1 = 1.0e7;

Qext2 = 1.0e7;

Q_tot = 1/(1/Q_int + 1/Qext1 + 1/Qext2);

nu0 = f_center;

t_ph = Q_tot/(2*pi*nu0);

alpha = t_R/(2*t_ph);

t_e1 = Qext1/(2*pi*nu0);

theta = t_R/t_e1;

fprintf('==== Physical parameters ====\n');
fprintf('Q_int = %.2e\n',Q_int);
fprintf('Q_ext1 = %.2e\n',Qext1);
fprintf('Q_ext2 = %.2e\n',Qext2);
fprintf('Q_tot = %.2e\n',Q_tot);
fprintf('FSR = %.3f GHz\n',FSR*1e-9);
fprintf('t_R = %.3f ps\n',t_R*1e12);
fprintf('t_ph = %.3f ns\n',t_ph*1e9);
fprintf('alpha = %.3e\n',alpha);
fprintf('theta = %.3e\n',theta);
fprintf('Ring length L = %.3f mm\n\n',L*1e3);

%% =========================================================
% 3.3 SLOW-TIME INTEGRATION
% =========================================================

symbol_rate = 1.0e9;

T_sym = 1/symbol_rate;

N_substeps_sym = 1000;

dt_slow = T_sym/N_substeps_sym;

fprintf('Data rate = %.1f GSa/s\n',symbol_rate/1e9);
fprintf('T_sym = %.3f ns\n',T_sym*1e9);
fprintf('dt_slow = %.3f ps\n',dt_slow*1e12);
fprintf('RK4IP steps per symbol = %d\n\n',N_substeps_sym);

params.t_R = t_R;
params.L = L;
params.Beta2 = Beta_2;
params.gamma = gamma;
params.alpha = alpha;
params.theta = theta;

%% =========================================================
% 4. PREPARATION PARAMETERS
% =========================================================

initial_noise_seed = 1;

N_scan_sym = 600;

N_hold_sym = 400;

delta_start = -2e-2;

Ein_bias = sqrt(Pav);

%% =========================================================
% 5. FTMLE PARAMETERS
% =========================================================

perturbation_seed = 2026;

perturbation_relative = 1e-8;

%% =========================================================
% 6. RESULT AND CACHE DIRECTORIES
% =========================================================

result_dir = fullfile(pwd,'narma10_ftmle_soliton_seed1');

if ~exist(result_dir,'dir')
    mkdir(result_dir);
end

cache_dir = fullfile(result_dir,'cache');

if ~exist(cache_dir,'dir')
    mkdir(cache_dir);
end

prepared_file = fullfile(cache_dir,'prepared_soliton_seed1.mat');

stage2_cache_file = fullfile(cache_dir,'stage2_soliton_seed1.mat');

result_file = fullfile(result_dir,'narma10_ftmle_soliton_seed1.mat');

%% =========================================================
% 7. CONFIGURATION SIGNATURES
% =========================================================

prepared_cfg.Pav = Pav;
prepared_cfg.FSR = FSR;
prepared_cfg.N = N;
prepared_cfg.Beta_2 = Beta_2;
prepared_cfg.gamma = gamma;
prepared_cfg.delta0_target = delta0_target;
prepared_cfg.neff = neff;
prepared_cfg.Q_int = Q_int;
prepared_cfg.Qext1 = Qext1;
prepared_cfg.Qext2 = Qext2;
prepared_cfg.symbol_rate = symbol_rate;
prepared_cfg.N_substeps_sym = N_substeps_sym;
prepared_cfg.dt_slow = dt_slow;
prepared_cfg.N_scan_sym = N_scan_sym;
prepared_cfg.N_hold_sym = N_hold_sym;
prepared_cfg.delta_start = delta_start;
prepared_cfg.initial_noise_seed = initial_noise_seed;

stage2_cfg = prepared_cfg;

stage2_cfg.input_seed = input_seed;
stage2_cfg.tau_pred = tau_pred;
stage2_cfg.N_washout = N_washout;
stage2_cfg.N_train = N_train;
stage2_cfg.N_test = N_test;
stage2_cfg.N_steps = N_steps;
stage2_cfg.Delta = Delta;
stage2_cfg.perturbation_seed = perturbation_seed;
stage2_cfg.perturbation_relative = perturbation_relative;

%% =========================================================
% 8. FIXED TARGET-DETUNING RK4IP OPERATOR
% =========================================================

rk_cache_target = rk4ip_lle_prepare(delta0_target,params,dt_slow,OM);

%% =========================================================
% 9. TIME AXIS
% =========================================================

time_s = (1:N_steps).'*T_sym;

time_ns = time_s*1e9;

%% =========================================================
% 10. CHECK COMPLETE STAGE-2 CACHE
% =========================================================

use_stage2_cache = false;

stage2_cache_source = '';

if ~force_recompute_stage2 && exist(stage2_cache_file,'file')

    tmp_cfg = load(stage2_cache_file,'stage2_cfg');

    if isfield(tmp_cfg,'stage2_cfg')

        if isequaln(tmp_cfg.stage2_cfg,stage2_cfg)

            use_stage2_cache = true;

            stage2_cache_source = stage2_cache_file;

        end

    end

end

%% =========================================================
% 11. CHECK LEGACY RESULT MAT
% =========================================================

if ~use_stage2_cache && ~force_recompute_stage2 && exist(result_file,'file')

    fprintf('Checking existing legacy result MAT...\n');

    required_names = {'X','E_prepared','time_s','time_ns','log_growth','distance_before_renorm','FTMLE_s_inv','FTMLE_ns_inv','d0'};

    file_info = whos('-file',result_file);

    file_names = {file_info.name};

    has_required_data = all(ismember(required_names,file_names));

    if has_required_data

        meta = load(result_file,'initial_noise_seed','perturbation_seed','perturbation_relative','N_scan_sym','N_hold_sym','N_washout','N_train','N_test','N_steps','Pav','Delta','delta0_target','symbol_rate','T_sym','dt_slow');

        legacy_ok = true;

        if ~isfield(meta,'initial_noise_seed') || meta.initial_noise_seed ~= initial_noise_seed
            legacy_ok = false;
        end

        if ~isfield(meta,'perturbation_seed') || meta.perturbation_seed ~= perturbation_seed
            legacy_ok = false;
        end

        if ~isfield(meta,'perturbation_relative') || meta.perturbation_relative ~= perturbation_relative
            legacy_ok = false;
        end

        if ~isfield(meta,'N_scan_sym') || meta.N_scan_sym ~= N_scan_sym
            legacy_ok = false;
        end

        if ~isfield(meta,'N_hold_sym') || meta.N_hold_sym ~= N_hold_sym
            legacy_ok = false;
        end

        if ~isfield(meta,'N_washout') || meta.N_washout ~= N_washout
            legacy_ok = false;
        end

        if ~isfield(meta,'N_train') || meta.N_train ~= N_train
            legacy_ok = false;
        end

        if ~isfield(meta,'N_test') || meta.N_test ~= N_test
            legacy_ok = false;
        end

        if ~isfield(meta,'N_steps') || meta.N_steps ~= N_steps
            legacy_ok = false;
        end

        if ~isfield(meta,'Pav') || meta.Pav ~= Pav
            legacy_ok = false;
        end

        if ~isfield(meta,'Delta') || meta.Delta ~= Delta
            legacy_ok = false;
        end

        if ~isfield(meta,'delta0_target') || meta.delta0_target ~= delta0_target
            legacy_ok = false;
        end

        if ~isfield(meta,'symbol_rate') || meta.symbol_rate ~= symbol_rate
            legacy_ok = false;
        end

        if ~isfield(meta,'T_sym') || meta.T_sym ~= T_sym
            legacy_ok = false;
        end

        if ~isfield(meta,'dt_slow') || meta.dt_slow ~= dt_slow
            legacy_ok = false;
        end

        if legacy_ok

            use_stage2_cache = true;

            stage2_cache_source = result_file;

            fprintf('Existing result MAT is compatible.\n');
            fprintf('The previous RK4IP result will be reused.\n\n');

        else

            fprintf('Existing result MAT is not compatible with current parameters.\n\n');

        end

    else

        fprintf('Existing result MAT does not contain all required simulation data.\n\n');

    end

end

%% =========================================================
% 12. LOAD COMPLETE STAGE-2 CACHE WHEN AVAILABLE
% =========================================================

if use_stage2_cache

    fprintf('=========================================================\n');
    fprintf('LOADING SAVED STAGE-2 DATA\n');
    fprintf('=========================================================\n');

    fprintf('Source:\n%s\n\n',stage2_cache_source);

    tmp = load(stage2_cache_source,'X','E_prepared','time_s','time_ns','log_growth','distance_before_renorm','FTMLE_s_inv','FTMLE_ns_inv','d0');

    X = tmp.X;

    E_prepared = tmp.E_prepared;

    time_s = tmp.time_s;

    time_ns = tmp.time_ns;

    log_growth = tmp.log_growth;

    distance_before_renorm = tmp.distance_before_renorm;

    FTMLE_s_inv = tmp.FTMLE_s_inv;

    FTMLE_ns_inv = tmp.FTMLE_ns_inv;

    d0 = tmp.d0;

    FTMLE_final_s_inv = FTMLE_s_inv(end);

    FTMLE_final_ns_inv = FTMLE_ns_inv(end);

    fprintf('Loaded %d Stage-2 symbols.\n',length(FTMLE_ns_inv));
    fprintf('Final FTMLE = %+.8e ns^-1\n',FTMLE_final_ns_inv);
    fprintf('No preparation or Stage-2 RK4IP calculation is required.\n');
    fprintf('=========================================================\n\n');

    if strcmp(stage2_cache_source,result_file)

        save(stage2_cache_file,'X','E_prepared','time_s','time_ns','log_growth','distance_before_renorm','FTMLE_s_inv','FTMLE_ns_inv','FTMLE_final_s_inv','FTMLE_final_ns_inv','d0','stage2_cfg','-v7.3');

        fprintf('New Stage-2 cache created from legacy result:\n');
        fprintf('%s\n\n',stage2_cache_file);

    end

else

    %% =====================================================
    % 13. CHECK PREPARED SOLITON CACHE
    % ======================================================

    use_prepared_cache = false;

    if ~force_reprepare && exist(prepared_file,'file')

        tmp_cfg = load(prepared_file,'prepared_cfg');

        if isfield(tmp_cfg,'prepared_cfg')

            if isequaln(tmp_cfg.prepared_cfg,prepared_cfg)

                use_prepared_cache = true;

            end

        end

    end

    %% =====================================================
    % 14. LOAD OR GENERATE PREPARED SOLITON STATE
    % ======================================================

    if use_prepared_cache

        fprintf('=========================================================\n');
        fprintf('LOADING PREPARED SOLITON STATE\n');
        fprintf('=========================================================\n');

        tmp = load(prepared_file,'E_prepared');

        E_prepared = tmp.E_prepared;

        fprintf('Prepared soliton state loaded.\n');
        fprintf('600-symbol scan and 400-symbol hold are skipped.\n');
        fprintf('=========================================================\n\n');

    else

        %% =================================================
        % 14.1 INITIAL INTRACAVITY FIELD
        % ==================================================

        rng(initial_noise_seed,'twister');

        E = 1e-6*(randn(1,N) + 1i*randn(1,N));

        %% =================================================
        % 14.2 STAGE 1: DETUNING SCAN
        % ==================================================

        fprintf('=========================================================\n');
        fprintf('Stage 1: detuning scan\n');
        fprintf('=========================================================\n');
        fprintf('Initial noise seed = %d\n',initial_noise_seed);
        fprintf('delta0: %.4g -> %.4g\n',delta_start,delta0_target);
        fprintf('Scan length = %d symbols\n\n',N_scan_sym);

        for n = 1:N_scan_sym

            frac = (n-1)/max(N_scan_sym-1,1);

            delta0 = delta_start + (delta0_target-delta_start)*frac;

            rk_cache_scan = rk4ip_lle_prepare(delta0,params,dt_slow,OM);

            drive_bias_scan = rk_cache_scan.drive_scale*Ein_bias;

            for kk = 1:N_substeps_sym

                E = rk4ip_lle_step_fast(E,drive_bias_scan,rk_cache_scan);

            end

            if mod(n,100) == 0 || n == 1 || n == N_scan_sym

                fprintf('Scan %4d / %4d | delta0 = %.6e\n',n,N_scan_sym,delta0);

            end

        end

        fprintf('\nScan finished.\n\n');

        %% =================================================
        % 14.3 HOLD STAGE
        % ==================================================

        drive_bias_hold = rk_cache_target.drive_scale*Ein_bias;

        fprintf('=========================================================\n');
        fprintf('Hold stage\n');
        fprintf('=========================================================\n');
        fprintf('Fixed delta0 = %.6e\n',delta0_target);
        fprintf('Fixed pump power = %.6f W\n',Pav);
        fprintf('Hold length = %d symbols\n\n',N_hold_sym);

        for n = 1:N_hold_sym

            for kk = 1:N_substeps_sym

                E = rk4ip_lle_step_fast(E,drive_bias_hold,rk_cache_target);

            end

            if mod(n,100) == 0 || n == 1 || n == N_hold_sym

                fprintf('Hold %4d / %4d\n',n,N_hold_sym);

            end

        end

        fprintf('\nHold finished.\n\n');

        %% =================================================
        % 14.4 SAVE PREPARED SOLITON STATE
        % ==================================================

        E_prepared = E;

        save(prepared_file,'E_prepared','prepared_cfg','-v7.3');

        fprintf('Prepared soliton state saved:\n');
        fprintf('%s\n\n',prepared_file);

    end

    %% =====================================================
    % 15. FTMLE INITIALIZATION
    % ======================================================

    E_ref = E_prepared;

    rs_pert = RandStream('mt19937ar','Seed',perturbation_seed);

    v0_real = randn(rs_pert,1,N);

    v0_imag = randn(rs_pert,1,N);

    v0 = v0_real + 1i*v0_imag;

    v0_norm = norm(v0,2);

    if ~isfinite(v0_norm) || v0_norm <= 1e-15
        error('Invalid perturbation direction.');
    end

    v0 = v0/v0_norm;

    E_ref_norm = norm(E_ref,2);

    if ~isfinite(E_ref_norm) || E_ref_norm <= 1e-15
        error('Invalid prepared soliton field norm.');
    end

    d0 = perturbation_relative*E_ref_norm;

    E_pert = E_ref + d0*v0;

    fprintf('=========================================================\n');
    fprintf('FTMLE initialization\n');
    fprintf('=========================================================\n');
    fprintf('Initial noise seed = %d\n',initial_noise_seed);
    fprintf('Perturbation seed = %d\n',perturbation_seed);
    fprintf('Relative perturbation = %.3e\n',perturbation_relative);
    fprintf('Prepared field norm = %.8e\n',E_ref_norm);
    fprintf('Initial perturbation d0 = %.8e\n\n',d0);

    %% =====================================================
    % 16. ALLOCATE STAGE-2 DATA
    % ======================================================

    X = zeros(N_steps,N);

    time_s = (1:N_steps).'*T_sym;

    time_ns = time_s*1e9;

    distance_before_renorm = zeros(N_steps,1);

    log_growth = zeros(N_steps,1);

    FTMLE_s_inv = zeros(N_steps,1);

    FTMLE_ns_inv = zeros(N_steps,1);

    cumulative_log_growth = 0;

    fprintf('=========================================================\n');
    fprintf('Stage 2: NARMA10 + FTMLE\n');
    fprintf('=========================================================\n');
    fprintf('Washout = %d symbols\n',N_washout);
    fprintf('Training = %d symbols\n',N_train);
    fprintf('Testing = %d symbols\n',N_test);
    fprintf('Total = %d symbols\n',N_steps);
    fprintf('Physical FTMLE window = %.3f us\n\n',N_steps*T_sym*1e6);

    %% =====================================================
    % 17. STAGE 2: NARMA10 + FTMLE
    % ======================================================

    for n = 1:N_steps

        Ein_n = Ein_series(n);

        drive_n = rk_cache_target.drive_scale*Ein_n;

        for kk = 1:N_substeps_sym

            E_ref = rk4ip_lle_step_fast(E_ref,drive_n,rk_cache_target);

            E_pert = rk4ip_lle_step_fast(E_pert,drive_n,rk_cache_target);

        end

        F_spec = fftshift(fft(E_ref));

        F_filt = F_spec;

        % To enable the filter:
        % F_filt = F_spec .* H_filter;

        X(n,:) = abs(F_filt).^2;

        delta_E = E_pert-E_ref;

        d_n = norm(delta_E,2);

        if ~isfinite(d_n) || d_n <= 1e-300
            error('Invalid perturbation distance at symbol %d.',n);
        end

        distance_before_renorm(n) = d_n;

        g_n = log(d_n/d0);

        log_growth(n) = g_n;

        cumulative_log_growth = cumulative_log_growth + g_n;

        current_time_s = n*T_sym;

        FTMLE_s_inv(n) = cumulative_log_growth/current_time_s;

        FTMLE_ns_inv(n) = FTMLE_s_inv(n)*1e-9;

        perturbation_direction = delta_E/d_n;

        E_pert = E_ref + d0*perturbation_direction;

        if n == 1 || mod(n,100) == 0 || n == N_steps

            fprintf('Symbol %4d / %4d | t = %7.1f ns | FTMLE = %+.8e ns^-1\n',n,N_steps,time_ns(n),FTMLE_ns_inv(n));

        end

    end

    fprintf('\nStage 2 finished.\n\n');

    %% =====================================================
    % 18. FINAL FTMLE
    % ======================================================

    FTMLE_final_s_inv = FTMLE_s_inv(end);

    FTMLE_final_ns_inv = FTMLE_ns_inv(end);

    E_ref_final = E_ref;

    E_pert_final = E_pert;

    %% =====================================================
    % 19. SAVE HEAVY STAGE-2 CACHE
    % ======================================================

    save(stage2_cache_file,'X','E_prepared','E_ref_final','E_pert_final','time_s','time_ns','log_growth','distance_before_renorm','FTMLE_s_inv','FTMLE_ns_inv','FTMLE_final_s_inv','FTMLE_final_ns_inv','d0','stage2_cfg','-v7.3');

    fprintf('=========================================================\n');
    fprintf('STAGE-2 CACHE SAVED\n');
    fprintf('=========================================================\n');
    fprintf('%s\n',stage2_cache_file);
    fprintf('Future runs can skip all RK4IP calculations.\n');
    fprintf('=========================================================\n\n');

end

%% =========================================================
% 20. PREPARED SOLITON STATE
% =========================================================

F_after_prepare = fftshift(fft(E_prepared));

E_plot = fftshift(E_prepared);

I_t_prepare = abs(E_plot).^2;

I_f_prepare = abs(F_after_prepare).^2;

%% =========================================================
% 21. FTMLE SUMMARY VALUES
% =========================================================

washout_end_index = N_washout;

training_end_index = N_washout + N_train;

FTMLE_washout_end_ns_inv = FTMLE_ns_inv(washout_end_index);

FTMLE_training_end_ns_inv = FTMLE_ns_inv(training_end_index);

FTMLE_final_s_inv = FTMLE_s_inv(end);

FTMLE_final_ns_inv = FTMLE_ns_inv(end);

washout_end_ns = N_washout*T_sym*1e9;

training_end_ns = (N_washout + N_train)*T_sym*1e9;

total_end_ns = N_steps*T_sym*1e9;

fprintf('=========================================================\n');
fprintf('FTMLE RESULTS\n');
fprintf('=========================================================\n');
fprintf('FTMLE at washout end (%.0f ns) = %+.8e ns^-1\n',washout_end_ns,FTMLE_washout_end_ns_inv);
fprintf('FTMLE at training end (%.0f ns) = %+.8e ns^-1\n',training_end_ns,FTMLE_training_end_ns_inv);
fprintf('Final FTMLE (%.0f ns) = %+.8e ns^-1\n',total_end_ns,FTMLE_final_ns_inv);
fprintf('Final FTMLE = %+.8e s^-1\n',FTMLE_final_s_inv);
fprintf('=========================================================\n\n');

%% =========================================================
% 22. CONSTRUCT TRAINING AND TESTING DATA
% =========================================================

X_start = N_washout + 1;

X_end = N_washout + N_train + N_test;

X_eff = X(X_start:X_end,:);

y_start = N_washout + 1 + tau_pred;

y_end = N_washout + N_train + N_test + tau_pred;

y_eff = y(y_start:y_end);

X_train = X_eff(1:N_train,:);

y_train = y_eff(1:N_train);

X_test_start = N_train + 1;

X_test_end = N_train + N_test;

X_test = X_eff(X_test_start:X_test_end,:);

y_test = y_eff(X_test_start:X_test_end);

%% =========================================================
% 23. FEATURE NORMALIZATION
% =========================================================

[X_train_n,mu_X,sigma_X] = zscore_features(X_train);

X_test_n = (X_test-mu_X)./sigma_X;

bias_train = ones(N_train,1);

bias_test = ones(N_test,1);

Xtr = [X_train_n,bias_train];

Xte = [X_test_n,bias_test];

%% =========================================================
% 24. RIDGE REGRESSION LAMBDA SCAN
% =========================================================

lambdas = logspace(-9,2,200);

nmse_scan = zeros(size(lambdas));

p = size(Xtr,2);

XtX = Xtr.'*Xtr;

Xty = Xtr.'*y_train;

Ireg = eye(p);

Ireg(end,end) = 0;

var_y_train = var(y_train,1);

for ii = 1:length(lambdas)

    lam = lambdas(ii);

    ridge_matrix = XtX + lam*Ireg;

    beta_tmp = ridge_matrix\Xty;

    y_train_pred_tmp = Xtr*beta_tmp;

    error_tmp = y_train_pred_tmp-y_train;

    nmse_scan(ii) = mean(error_tmp.^2)/var_y_train;

end

[best_nmse,best_idx] = min(nmse_scan);

lambda_opt = lambdas(best_idx);

fprintf('=========================================================\n');
fprintf('RIDGE REGRESSION\n');
fprintf('=========================================================\n');
fprintf('Optimal lambda = %.8e\n',lambda_opt);
fprintf('Training NMSE during scan = %.8e\n\n',best_nmse);

%% =========================================================
% FIGURE 1
% LAMBDA SCAN
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','Lambda vs Train NMSE');

semilogx(lambdas,nmse_scan,'-o','LineWidth',1.2,'MarkerSize',6);

grid on;

hold on;

semilogx(lambda_opt,best_nmse,'s','MarkerSize',8,'LineWidth',1.5);

ylambda = ylim;

plot([lambda_opt lambda_opt],ylambda,'k--','LineWidth',1.0);

ylim(ylambda);

xlabel('\lambda');

ylabel('Train NMSE');

title('Ridge: Train NMSE vs \lambda');

legend('Scan NMSE','Best point','Best \lambda','Location','best');
%}

%% =========================================================
% 26. FINAL RIDGE REGRESSION
% =========================================================

ridge_matrix = XtX + lambda_opt*Ireg;

beta = ridge_matrix\Xty;

y_train_pred = Xtr*beta;

y_test_pred = Xte*beta;

train_error = y_train_pred-y_train;

test_error = y_test_pred-y_test;

nmse_train = mean(train_error.^2)/var(y_train,1);

nmse_test = mean(test_error.^2)/var(y_test,1);

fprintf('Training NMSE = %.8e\n',nmse_train);

fprintf('Testing NMSE = %.8e\n',nmse_test);

fprintf('=========================================================\n\n');

%% =========================================================
% FIGURE 2
% PREPARED SOLITON STATE
% COMMENTED OUT
% =========================================================
%{
lam_axis = c0./f_axis*1e9;

figure('Color','w','Name','Prepared soliton state');

subplot(2,1,1);

plot(tau*1e12,10*log10(I_t_prepare+1e-12),'b','LineWidth',1);

xlabel('Fast time \tau (ps)');

ylabel('Intensity (dB)');

title('Prepared soliton: intracavity temporal intensity');

grid on;

subplot(2,1,2);

plot(lam_axis,10*log10(I_f_prepare+1e-12),'r','LineWidth',1);

xlabel('Wavelength (nm)');

ylabel('Power (dB)');

title('Prepared soliton: optical spectrum');

xlim([1500 1700]);

grid on;
%}

%% =========================================================
% FIGURE 3
% NARMA10 PREDICTION
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','NARMA10 one-step prediction');

plot(y_test,'k','LineWidth',1);

hold on;

plot(y_test_pred,'r--','LineWidth',1);

xlabel('Time step');

ylabel('NARMA10 output');

legend('Ground truth','Prediction','Location','best');

title(sprintf('NARMA10 one-step prediction, NMSE = %.3e',nmse_test));

grid on;
%}

%% =========================================================
% FIGURE 4
% RESERVOIR SPECTRAL STATES
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','Reservoir spectral states');

X_plot_dB = 10*log10(X_eff.'+1e-12);

imagesc(X_plot_dB);

axis xy;

xlabel('Time index after washout');

ylabel('Spectral channel index');

title('Microcomb reservoir spectral states');

colorbar;
%}

%% =========================================================
% FIGURE 5
% NARMA10 INPUT AND TARGET
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','NARMA10 input and target');

subplot(2,1,1);

plot(u(1:N_steps),'b','LineWidth',1);

xlabel('Time step');

ylabel('u(k)');

title('NARMA10 input');

grid on;

subplot(2,1,2);

plot(y(1:N_steps+tau_pred),'k','LineWidth',1);

xlabel('Time step');

ylabel('y(k)');

title('NARMA10 target');

grid on;
%}

%% =========================================================
% FIGURE 6
% FTMLE PLOT
% ENABLED
% =========================================================

figure('Color','w','Name','Soliton FTMLE','Position',[100 100 850 560]);

plot(time_ns,FTMLE_ns_inv,'LineWidth',1.5);

hold on;

yline(0,'k--','LineWidth',1.0,'HandleVisibility','off');

xline(washout_end_ns,'--','Washout end','LineWidth',1.0,'HandleVisibility','off');

xline(training_end_ns,'--','Training end','LineWidth',1.0,'HandleVisibility','off');

xlabel('Time (ns)','FontName','Times New Roman','FontSize',16);

ylabel('FTMLE (ns^{-1})','FontName','Times New Roman','FontSize',16);

title('Finite-time maximum Lyapunov exponent','FontName','Times New Roman','FontSize',16,'FontWeight','normal');

xlim([0 total_end_ns]);

grid on;

set(gca,'FontName','Times New Roman');

set(gca,'FontSize',14);

set(gca,'LineWidth',1.0);

set(gca,'Box','on');

set(gca,'TickDir','out');

%% =========================================================
% FIGURE 7
% LOCAL LOGARITHMIC GROWTH
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','Local logarithmic growth');

plot(time_ns,log_growth,'LineWidth',1.0);

hold on;

yline(0,'k--','LineWidth',1.0);

xline(washout_end_ns,'--','Washout end','LineWidth',1.0);

xline(training_end_ns,'--','Training end','LineWidth',1.0);

xlabel('Time (ns)');

ylabel('log(d_n/d_0)');

title('Local logarithmic perturbation growth');

xlim([0 total_end_ns]);

grid on;
%}

%% =========================================================
% FIGURE 8
% PERTURBATION DISTANCE
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','Perturbation distance');

semilogy(time_ns,distance_before_renorm,'LineWidth',1.2);

hold on;

yline(d0,'k--','d_0','LineWidth',1.0);

xline(washout_end_ns,'--','Washout end','LineWidth',1.0);

xline(training_end_ns,'--','Training end','LineWidth',1.0);

xlabel('Time (ns)');

ylabel('Distance before renormalization');

title('Perturbation distance before Benettin renormalization');

xlim([0 total_end_ns]);

grid on;
%}

%% =========================================================
% 34. CUMULATIVE LOGARITHMIC GROWTH DATA
% =========================================================

cumulative_log_growth_curve = cumsum(log_growth);

%% =========================================================
% FIGURE 9
% CUMULATIVE LOGARITHMIC GROWTH
% COMMENTED OUT
% =========================================================
%{
figure('Color','w','Name','Cumulative logarithmic growth');

plot(time_ns,cumulative_log_growth_curve,'LineWidth',1.4);

hold on;

yline(0,'k--','LineWidth',1.0);

xline(washout_end_ns,'--','Washout end','LineWidth',1.0);

xline(training_end_ns,'--','Training end','LineWidth',1.0);

xlabel('Time (ns)');

ylabel('\Sigma log(d_n/d_0)');

title('Cumulative logarithmic perturbation growth');

xlim([0 total_end_ns]);

grid on;
%}

%% =========================================================
% 35. SAVE COMPLETE RESULT MAT
% =========================================================

save(result_file,'initial_noise_seed','input_seed','perturbation_seed','perturbation_relative','d0','u','y','Pin_series','Ein_series','X','X_eff','y_train','y_test','y_train_pred','y_test_pred','nmse_train','nmse_test','lambda_opt','time_s','time_ns','log_growth','cumulative_log_growth_curve','distance_before_renorm','FTMLE_s_inv','FTMLE_ns_inv','FTMLE_washout_end_ns_inv','FTMLE_training_end_ns_inv','FTMLE_final_s_inv','FTMLE_final_ns_inv','N_scan_sym','N_hold_sym','N_washout','N_train','N_test','N_steps','Pav','Delta','delta0_target','symbol_rate','T_sym','dt_slow','E_prepared','prepared_cfg','stage2_cfg','-v7.3');

fprintf('Complete result MAT saved:\n');

fprintf('%s\n\n',result_file);

%% =========================================================
% 36. SAVE FTMLE CURVE CSV
% =========================================================

Time_ns = time_ns;

FTMLE_ns_inverse = FTMLE_ns_inv;

LogGrowth = log_growth;

CumulativeLogGrowth = cumulative_log_growth_curve;

DistanceBeforeRenorm = distance_before_renorm;

FTMLE_table = table(Time_ns,FTMLE_ns_inverse,LogGrowth,CumulativeLogGrowth,DistanceBeforeRenorm);

ftmle_csv = fullfile(result_dir,'FTMLE_soliton_seed1_curve.csv');

writetable(FTMLE_table,ftmle_csv);

%% =========================================================
% 37. SAVE SUMMARY CSV
% =========================================================

NoiseSeed = initial_noise_seed;

ScanSymbols = N_scan_sym;

HoldSymbols = N_hold_sym;

WashoutSymbols = N_washout;

TrainSymbols = N_train;

TestSymbols = N_test;

TrainNMSE = nmse_train;

TestNMSE = nmse_test;

FTMLE_WashoutEnd_ns_inv = FTMLE_washout_end_ns_inv;

FTMLE_TrainingEnd_ns_inv = FTMLE_training_end_ns_inv;

FTMLE_Final_ns_inv = FTMLE_final_ns_inv;

summary_table = table(NoiseSeed,ScanSymbols,HoldSymbols,WashoutSymbols,TrainSymbols,TestSymbols,TrainNMSE,TestNMSE,FTMLE_WashoutEnd_ns_inv,FTMLE_TrainingEnd_ns_inv,FTMLE_Final_ns_inv);

summary_csv = fullfile(result_dir,'soliton_seed1_summary.csv');

writetable(summary_table,summary_csv);

%% =========================================================
% 38. FINAL SUMMARY
% =========================================================

fprintf('============================================================\n');

fprintf(' FINAL SUMMARY\n');

fprintf('============================================================\n');

fprintf('Operating regime         = Soliton\n');

fprintf('Initial noise seed       = %d\n',initial_noise_seed);

fprintf('Scan                     = %d symbols\n',N_scan_sym);

fprintf('Hold                     = %d symbols\n',N_hold_sym);

fprintf('Washout                  = %d symbols\n',N_washout);

fprintf('Training                 = %d symbols\n',N_train);

fprintf('Testing                  = %d symbols\n',N_test);

fprintf('Stage-2 total            = %d symbols\n',N_steps);

fprintf('Stage-2 physical window  = %.3f us\n',N_steps*T_sym*1e6);

fprintf('\n');

fprintf('NARMA10 train NMSE       = %.8e\n',nmse_train);

fprintf('NARMA10 test NMSE        = %.8e\n',nmse_test);

fprintf('\n');

fprintf('FTMLE at %.0f ns          = %+.8e ns^-1\n',washout_end_ns,FTMLE_washout_end_ns_inv);

fprintf('FTMLE at %.0f ns         = %+.8e ns^-1\n',training_end_ns,FTMLE_training_end_ns_inv);

fprintf('FTMLE at %.0f ns         = %+.8e ns^-1\n',total_end_ns,FTMLE_final_ns_inv);

fprintf('\n');

fprintf('Prepared-state cache:\n');

fprintf('%s\n\n',prepared_file);

fprintf('Stage-2 cache:\n');

fprintf('%s\n\n',stage2_cache_file);

fprintf('Complete result MAT:\n');

fprintf('%s\n',result_file);

fprintf('============================================================\n');

end

%% =========================================================
% RK4IP OPERATOR PREPARATION
% =========================================================

function cache = rk4ip_lle_prepare(delta0,params,dt,OM)

t_R = params.t_R;

L = params.L;

Beta2 = params.Beta2;

gamma = params.gamma;

alpha = params.alpha;

theta = params.theta;

loss_term = -(alpha + 1i*delta0)/t_R;

dispersion_term = 1i*(L/t_R)*(Beta2/2).*(OM.^2);

D_lin = loss_term + dispersion_term;

cache.expD_h2 = ifftshift(exp(D_lin*(dt/2)));

cache.expD_h = ifftshift(exp(D_lin*dt));

cache.expMinusD_h2 = ifftshift(exp(-D_lin*(dt/2)));

cache.expMinusD_h = ifftshift(exp(-D_lin*dt));

cache.c_nl = 1i*gamma*L/t_R;

cache.drive_scale = sqrt(theta)/t_R;

cache.dt = dt;

cache.dt_h = dt/2;

cache.dt_6 = dt/6;

end

%% =========================================================
% ONE RK4IP STEP
% =========================================================

function E_next = rk4ip_lle_step_fast(E,drive,cache)

dt = cache.dt;

dt_h = cache.dt_h;

dt_6 = cache.dt_6;

c_nl = cache.c_nl;

expD_h2 = cache.expD_h2;

expD_h = cache.expD_h;

expMinusD_h2 = cache.expMinusD_h2;

expMinusD_h = cache.expMinusD_h;

U0 = E;

%% K1

N1 = c_nl.*abs(U0).^2.*U0 + drive;

K1 = N1;

%% K2

Utemp2 = U0 + dt_h*K1;

E2 = ifft(fft(Utemp2).*expD_h2);

N2 = c_nl.*abs(E2).^2.*E2 + drive;

K2 = ifft(fft(N2).*expMinusD_h2);

%% K3

Utemp3 = U0 + dt_h*K2;

E3 = ifft(fft(Utemp3).*expD_h2);

N3 = c_nl.*abs(E3).^2.*E3 + drive;

K3 = ifft(fft(N3).*expMinusD_h2);

%% K4

Utemp4 = U0 + dt*K3;

E4 = ifft(fft(Utemp4).*expD_h);

N4 = c_nl.*abs(E4).^2.*E4 + drive;

K4 = ifft(fft(N4).*expMinusD_h);

%% RK4 COMBINATION

RK4_sum = K1 + 2*K2 + 2*K3 + K4;

U1 = U0 + dt_6*RK4_sum;

E_next = ifft(fft(U1).*expD_h);

end

%% =========================================================
% FEATURE NORMALIZATION
% =========================================================

function [Xn,mu_X,sigma_X] = zscore_features(X)

mu_X = mean(X,1);

sigma_X = std(X,0,1);

small_sigma = sigma_X < 1e-12;

sigma_X(small_sigma) = 1e-12;

Xn = (X-mu_X)./sigma_X;

end