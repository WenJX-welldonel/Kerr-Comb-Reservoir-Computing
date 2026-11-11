function kerr_comb_ftmle_narma10_seed1_washout800()

% =========================================================
% Kerr-comb FTMLE
%
% Operating regime:
%   Near-critical
%
% Initial intracavity noise seed:
%   seed = 1
%
% Protocol:
%   1. Generate initial small complex noise
%   2. Prepare near-critical state by 320-symbol scan
%   3. Define FTMLE t = 0 at first washout symbol
%   4. Propagate washout -> training -> testing
%   5. Use fixed NARMA10 input
%   6. Use fixed perturbation direction seed
%   7. Renormalize perturbation once per symbol
%
% Actual running window:
%   washout = 800
%   train   = 5360
%   test    = 800
%   total   = 6960 symbols
%
% At 1 GSa/s:
%   total physical window = 6960 ns = 6.96 us
%
% IMPORTANT:
%   Calculation parameters are unchanged.
%   Plotting range is generated automatically from the
%   actual calculated data.
% =========================================================

clc;
close all;

%% =========================================================
% 0. USER CONFIGURATION
% ==========================================================

% Initial intracavity noise seed
initial_noise_seed = 1;

% NARMA10 input seed
input_seed = 1;

% Perturbation-direction seed
perturbation_seed = 2026;

% Relative perturbation
perturbation_relative = 1e-8;

%% =========================================================
% 0.1 Plot settings
% ==========================================================

P.font_name = 'Times New Roman';

P.tick_font_size = 14;

P.label_font_size = 16;

P.title_font_size = 16;

P.legend_font_size = 12;

P.line_width = 1.6;

P.axis_line_width = 1.0;

P.grid_alpha = 0.18;

P.export_resolution = 600;

P.xlabel_text = 'Time (ns)';

P.ylabel_text = 'FTMLE (ns^{-1})';

P.title_text = 'Finite-time maximum Lyapunov exponent';

P.show_zero_line = true;

P.show_stage_boundaries = true;

% X-axis is now determined automatically from actual data.
P.x_tick_step = 1000;

% [] means automatic
P.y_lim = [];

P.y_ticks = [];

P.figure_position = [100 100 850 560];

%% =========================================================
% 0.2 Cache controls
% ==========================================================

% true = redo 320-symbol preparation
force_reprepare = false;

% true = redo FTMLE even if compatible cache exists
force_recompute_ftmle = false;

%% =========================================================
% 1. Output directories
% ==========================================================

result_dir = fullfile(pwd, 'ftmle_results_seed1_washout800');

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

prepared_dir = fullfile(result_dir, 'prepared_states');

if ~exist(prepared_dir, 'dir')
    mkdir(prepared_dir);
end

ftmle_cache_dir = fullfile(result_dir, 'ftmle_cache');

if ~exist(ftmle_cache_dir, 'dir')
    mkdir(ftmle_cache_dir);
end

%% =========================================================
% 2. Reservoir-computing time window
% ==========================================================

tau_pred = 1;

N_washout = 800;

N_train = 5360;

N_test = 800;

N_steps = N_washout + N_train + N_test;

%% =========================================================
% 3. Generate NARMA10 input
% ==========================================================

rs_input = RandStream('mt19937ar', 'Seed', input_seed);

N_total = N_steps + tau_pred + 50;

u_full = 0.5 * rand(rs_input, N_total, 1);

discard = 50;

u = u_full(discard+1:end);

if length(u) < N_steps
    error('Input sequence is too short.');
end

%% =========================================================
% 4. Input -> pump power
% ==========================================================

u_drive = 4*u - 1;

u_drive_used = u_drive(1:N_steps);

Pav = 0.325;

Delta = 0.05;

Pin_series = Pav + (Delta/2) .* u_drive_used;

Pin_series = max(Pin_series, 0);

Ein_series = sqrt(Pin_series);

%% =========================================================
% 5. Kerr-MRR / LLE parameters
% ==========================================================

c0 = 2.998e8;

nn = 2^8;

N = 2*nn;

FSR = 200e9;

t_R = 1/FSR;

D1 = 2*pi*FSR;

k = -nn:nn-1;

OM = D1*k;

lam_center = 1.600e-6;

f_center = c0/lam_center;

Beta_2 = -4.0e-9 * 1e-18;

gamma = 1.6e-2;

delta0_target = 4e-3;

neff = 1.444;

L = t_R*c0/neff;

Q_int = 1.0e7;

Qext1 = 1.0e7;

Qext2 = 1.0e7;

Q_tot = 1 / (1/Q_int + 1/Qext1 + 1/Qext2);

nu0 = f_center;

t_ph = Q_tot / (2*pi*nu0);

alpha = t_R / (2*t_ph);

t_e1 = Qext1 / (2*pi*nu0);

theta = t_R/t_e1;

%% =========================================================
% 6. Slow-time integration
% ==========================================================

symbol_rate = 1.0e9;

T_sym = 1/symbol_rate;

N_substeps_sym = 1000;

dt_slow = T_sym / N_substeps_sym;

params.t_R = t_R;

params.L = L;

params.Beta2 = Beta_2;

params.gamma = gamma;

params.alpha = alpha;

params.theta = theta;

%% =========================================================
% 7. Near-critical preparation protocol
% ==========================================================

N_scan_sym = 320;

delta_start = -2e-2;

Ein_bias = sqrt(Pav);

%% =========================================================
% 8. Configuration signatures
% ==========================================================

prepared_cfg.Pav = Pav;

prepared_cfg.FSR = FSR;

prepared_cfg.N = N;

prepared_cfg.Beta2 = Beta_2;

prepared_cfg.gamma = gamma;

prepared_cfg.delta0_target = delta0_target;

prepared_cfg.Q_int = Q_int;

prepared_cfg.Qext1 = Qext1;

prepared_cfg.Qext2 = Qext2;

prepared_cfg.symbol_rate = symbol_rate;

prepared_cfg.N_substeps_sym = N_substeps_sym;

prepared_cfg.dt_slow = dt_slow;

prepared_cfg.N_scan_sym = N_scan_sym;

prepared_cfg.delta_start = delta_start;

prepared_cfg.initial_noise_seed = initial_noise_seed;

ftmle_cfg = prepared_cfg;

ftmle_cfg.input_seed = input_seed;

ftmle_cfg.perturbation_seed = perturbation_seed;

ftmle_cfg.perturbation_relative = perturbation_relative;

ftmle_cfg.N_washout = N_washout;

ftmle_cfg.N_train = N_train;

ftmle_cfg.N_test = N_test;

ftmle_cfg.N_steps = N_steps;

ftmle_cfg.Delta = Delta;

%% =========================================================
% 9. Print protocol
% ==========================================================

fprintf('\n');

fprintf('=========================================================\n');

fprintf(' FTMLE - SINGLE INITIAL-NOISE REALIZATION\n');

fprintf('=========================================================\n');

fprintf('Initial noise seed       : %d\n', initial_noise_seed);

fprintf('Input seed               : %d\n', input_seed);

fprintf('Perturbation seed        : %d\n', perturbation_seed);

fprintf('Relative perturbation    : %.3e\n', perturbation_relative);

fprintf('\n');

fprintf('Preparation scan         : %d symbols\n', N_scan_sym);

fprintf('Washout                  : %d symbols\n', N_washout);

fprintf('Training                 : %d symbols\n', N_train);

fprintf('Testing                  : %d symbols\n', N_test);

fprintf('Total FTMLE window       : %d symbols\n', N_steps);

fprintf('Physical window          : %.3f us\n', N_steps*T_sym*1e6);

fprintf('=========================================================\n\n');

%% =========================================================
% 10. Fixed RK4IP operator for data stage
% ==========================================================

rk_cache_target = rk4ip_lle_prepare(delta0_target, params, dt_slow, OM);

%% =========================================================
% 11. Time axis
% ==========================================================

time_s = (1:N_steps).' * T_sym;

time_ns = time_s * 1e9;

%% =========================================================
% 12. Cache filenames
% ==========================================================

prepared_file = fullfile(prepared_dir, sprintf('prepared_nearcritical_noise_seed_%d.mat', initial_noise_seed));

ftmle_file = fullfile(ftmle_cache_dir, sprintf('FTMLE_noise_seed_%d.mat', initial_noise_seed));

%% =========================================================
% 13. Check completed FTMLE cache
% ==========================================================

use_ftmle_cache = false;

if ~force_recompute_ftmle && exist(ftmle_file, 'file')

    tmp_cfg = load(ftmle_file, 'ftmle_cfg');

    if isfield(tmp_cfg, 'ftmle_cfg') && isequaln(tmp_cfg.ftmle_cfg, ftmle_cfg)
        use_ftmle_cache = true;
    end

end

%% =========================================================
% 14. Load FTMLE result or calculate
% ==========================================================

if use_ftmle_cache

    fprintf('Compatible FTMLE cache found.\n');

    fprintf('Loading result; no RK4IP calculation needed.\n');

    tmp = load(ftmle_file, 'FTMLE_s_inv', 'FTMLE_ns_inv', 'FTMLE_final_s_inv', 'FTMLE_final_ns_inv', 'log_growth', 'distance_before_renorm', 'd0');

    FTMLE_s_inv = tmp.FTMLE_s_inv;

    FTMLE_ns_inv = tmp.FTMLE_ns_inv;

    FTMLE_final_s_inv = tmp.FTMLE_final_s_inv;

    FTMLE_final_ns_inv = tmp.FTMLE_final_ns_inv;

    log_growth = tmp.log_growth;

    distance_before_renorm = tmp.distance_before_renorm;

    d0 = tmp.d0;

    fprintf('Loaded final FTMLE = %+.8e ns^-1\n', FTMLE_final_ns_inv);

else

    %% =====================================================
    % 14.1 Check prepared-state cache
    % ======================================================

    use_prepared_cache = false;

    if ~force_reprepare && exist(prepared_file, 'file')

        tmp_cfg = load(prepared_file, 'prepared_cfg');

        if isfield(tmp_cfg, 'prepared_cfg') && isequaln(tmp_cfg.prepared_cfg, prepared_cfg)
            use_prepared_cache = true;
        end

    end

    %% =====================================================
    % 14.2 Load or prepare near-critical state
    % ======================================================

    if use_prepared_cache

        fprintf('Loading prepared near-critical state...\n');

        tmp = load(prepared_file, 'E_prepared');

        E_prepared = tmp.E_prepared;

        fprintf('Prepared state loaded.\n');

    else

        fprintf('Preparing near-critical state from noise seed %d...\n', initial_noise_seed);

        rng(initial_noise_seed, 'twister');

        E = 1e-6 * (randn(1,N) + 1i*randn(1,N));

        for n = 1:N_scan_sym

            frac = (n-1) / max(N_scan_sym-1, 1);

            delta0 = delta_start + (delta0_target - delta_start) * frac;

            rk_cache = rk4ip_lle_prepare(delta0, params, dt_slow, OM);

            drive_bias = rk_cache.drive_scale * Ein_bias;

            for kk = 1:N_substeps_sym
                E = rk4ip_lle_step_fast(E, drive_bias, rk_cache);
            end

        end

        E_prepared = E;

        save(prepared_file, 'E_prepared', 'prepared_cfg', '-v7.3');

        fprintf('Prepared state saved:\n%s\n', prepared_file);

    end

    %% =====================================================
    % 14.3 Reference trajectory at FTMLE t = 0
    % ======================================================

    E_ref = E_prepared;

    %% =====================================================
    % 14.4 Perturbation direction
    % ======================================================

    rs_pert = RandStream('mt19937ar', 'Seed', perturbation_seed);

    v0 = randn(rs_pert, 1, N) + 1i * randn(rs_pert, 1, N);

    v0 = v0 / norm(v0, 2);

    E_ref_norm = norm(E_ref, 2);

    if ~isfinite(E_ref_norm) || E_ref_norm <= 1e-15
        error('Invalid prepared field norm for noise seed %d.', initial_noise_seed);
    end

    d0 = perturbation_relative * E_ref_norm;

    E_pert = E_ref + d0*v0;

    fprintf('Prepared field norm      = %.6e\n', E_ref_norm);

    fprintf('Initial perturbation d0  = %.6e\n', d0);

    %% =====================================================
    % 14.5 Allocate FTMLE arrays
    % ======================================================

    log_growth = zeros(N_steps, 1);

    distance_before_renorm = zeros(N_steps, 1);

    FTMLE_s_inv = zeros(N_steps, 1);

    FTMLE_ns_inv = zeros(N_steps, 1);

    cumulative_log_growth = 0;

    %% =====================================================
    % 14.6 FTMLE propagation
    % ======================================================

    fprintf('Calculating FTMLE...\n');

    for n = 1:N_steps

        Ein_n = Ein_series(n);

        drive_n = rk_cache_target.drive_scale * Ein_n;

        % Propagate reference and perturbed trajectories
        for kk = 1:N_substeps_sym

            E_ref = rk4ip_lle_step_fast(E_ref, drive_n, rk_cache_target);

            E_pert = rk4ip_lle_step_fast(E_pert, drive_n, rk_cache_target);

        end

        % Perturbation distance
        delta_E = E_pert - E_ref;

        d_n = norm(delta_E, 2);

        if ~isfinite(d_n) || d_n <= 1e-300
            error('Invalid perturbation distance at seed %d, symbol %d.', initial_noise_seed, n);
        end

        distance_before_renorm(n) = d_n;

        % Local logarithmic expansion
        g_n = log(d_n / d0);

        log_growth(n) = g_n;

        cumulative_log_growth = cumulative_log_growth + g_n;

        % Cumulative FTMLE
        FTMLE_s_inv(n) = cumulative_log_growth / time_s(n);

        FTMLE_ns_inv(n) = FTMLE_s_inv(n) * 1e-9;

        % Benettin renormalization
        E_pert = E_ref + d0 * (delta_E / d_n);

        % Progress
        if mod(n,100) == 0 || n == 1 || n == N_steps

            fprintf('  noise seed %d | %4d / %4d | t = %7.1f ns | FTMLE = %+.6e ns^-1\n', initial_noise_seed, n, N_steps, time_ns(n), FTMLE_ns_inv(n));

        end

    end

    %% =====================================================
    % 14.7 Final FTMLE
    % ======================================================

    FTMLE_final_s_inv = FTMLE_s_inv(end);

    FTMLE_final_ns_inv = FTMLE_ns_inv(end);

    %% =====================================================
    % 14.8 Save FTMLE result
    % ======================================================

    noise_seed = initial_noise_seed;

    save(ftmle_file, 'noise_seed', 'input_seed', 'perturbation_seed', 'perturbation_relative', 'd0', 'time_s', 'time_ns', 'FTMLE_s_inv', 'FTMLE_ns_inv', 'FTMLE_final_s_inv', 'FTMLE_final_ns_inv', 'log_growth', 'distance_before_renorm', 'ftmle_cfg', '-v7.3');

    fprintf('FTMLE result saved:\n%s\n', ftmle_file);

end

%% =========================================================
% 15. Print final result
% ==========================================================

fprintf('\n');

fprintf('=========================================================\n');

fprintf(' SINGLE-SEED FTMLE RESULT\n');

fprintf('=========================================================\n');

fprintf('Noise seed %d : %+.8e ns^-1\n', initial_noise_seed, FTMLE_final_ns_inv);

fprintf('\n');

fprintf('FINAL FTMLE = %+.8e ns^-1\n', FTMLE_final_ns_inv);

fprintf('=========================================================\n\n');

%% =========================================================
% 16. Save summary MAT
% ==========================================================

summary_file = fullfile(result_dir, 'FTMLE_seed1_washout800_summary.mat');

save(summary_file, 'initial_noise_seed', 'input_seed', 'perturbation_seed', 'perturbation_relative', 'time_s', 'time_ns', 'FTMLE_s_inv', 'FTMLE_ns_inv', 'FTMLE_final_s_inv', 'FTMLE_final_ns_inv', 'log_growth', 'distance_before_renorm', 'd0', 'N_scan_sym', 'N_washout', 'N_train', 'N_test', 'N_steps', 'symbol_rate', 'T_sym', 'dt_slow', 'P', '-v7.3');

fprintf('Summary MAT saved:\n%s\n\n', summary_file);

%% =========================================================
% 17. Save final-value CSV
% ==========================================================

NoiseSeed = initial_noise_seed;

Final_FTMLE_ns_inv = FTMLE_final_ns_inv;

Final_FTMLE_s_inv = FTMLE_final_s_inv;

InitialPerturbation = d0;

T_final = table(NoiseSeed, Final_FTMLE_ns_inv, Final_FTMLE_s_inv, InitialPerturbation);

final_csv = fullfile(result_dir, 'FTMLE_seed1_washout800_final_value.csv');

writetable(T_final, final_csv);

%% =========================================================
% 18. Save time-dependent FTMLE CSV
% ==========================================================

Time_ns = time_ns;

FTMLE_seed1_ns_inv = FTMLE_ns_inv;

LogGrowth = log_growth;

DistanceBeforeRenorm = distance_before_renorm;

T_curve = table(Time_ns, FTMLE_seed1_ns_inv, LogGrowth, DistanceBeforeRenorm);

curve_csv = fullfile(result_dir, 'FTMLE_seed1_washout800_time_curve.csv');

writetable(T_curve, curve_csv);

%% =========================================================
% 19. Plot
% ==========================================================

fig = figure('Color', 'w', 'Name', 'FTMLE - seed 1, washout 800', 'Position', P.figure_position);

ax = axes(fig);

hold(ax, 'on');

colors = get(groot, 'defaultAxesColorOrder');

h_ftmle = plot(ax, time_ns, FTMLE_ns_inv, '-', 'Color', colors(1,:), 'LineWidth', P.line_width);

%% =========================================================
% 19.1 Zero line
% ==========================================================

if P.show_zero_line

    yline(ax, 0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

end

%% =========================================================
% 19.2 Stage boundaries
% ==========================================================

washout_end_ns = N_washout * T_sym * 1e9;

training_end_ns = (N_washout + N_train) * T_sym * 1e9;

total_end_ns = N_steps * T_sym * 1e9;

if P.show_stage_boundaries

    xline(ax, washout_end_ns, '--', 'Washout end', 'FontName', P.font_name, 'FontSize', P.tick_font_size, 'HandleVisibility', 'off');

    xline(ax, training_end_ns, '--', 'Training end', 'FontName', P.font_name, 'FontSize', P.tick_font_size, 'HandleVisibility', 'off');

end

%% =========================================================
% 19.3 Labels
% ==========================================================

xlabel(ax, P.xlabel_text, 'FontName', P.font_name, 'FontSize', P.label_font_size);

ylabel(ax, P.ylabel_text, 'FontName', P.font_name, 'FontSize', P.label_font_size);

title(ax, P.title_text, 'FontName', P.font_name, 'FontSize', P.title_font_size, 'FontWeight', 'normal');

%% =========================================================
% 19.4 Legend
% ==========================================================

lgd = legend(ax, h_ftmle, 'Noise seed 1', 'Location', 'best');

set(lgd, 'FontName', P.font_name, 'FontSize', P.legend_font_size);

%% =========================================================
% 19.5 Axis style
% ==========================================================

set(ax, 'FontName', P.font_name, 'FontSize', P.tick_font_size, 'LineWidth', P.axis_line_width, 'Box', 'on', 'TickDir', 'out');

grid(ax, 'on');

if isprop(ax, 'GridAlpha')
    ax.GridAlpha = P.grid_alpha;
end

%% =========================================================
% 19.6 X-axis range based on actual calculated data
% ==========================================================

xlim(ax, [0 total_end_ns]);

x_tick_values = 0:P.x_tick_step:total_end_ns;

if isempty(x_tick_values)
    x_tick_values = [0 total_end_ns];
elseif x_tick_values(end) < total_end_ns
    x_tick_values = [x_tick_values total_end_ns];
end

xticks(ax, x_tick_values);

%% =========================================================
% 19.7 Y-axis range
% ==========================================================

if ~isempty(P.y_lim)
    ylim(ax, P.y_lim);
end

if ~isempty(P.y_ticks)
    yticks(ax, P.y_ticks);
end

%% =========================================================
% 20. Save figure
% ==========================================================

png_file = fullfile(result_dir, 'FTMLE_seed1_washout800.png');

pdf_file = fullfile(result_dir, 'FTMLE_seed1_washout800.pdf');

fig_file = fullfile(result_dir, 'FTMLE_seed1_washout800.fig');

exportgraphics(fig, png_file, 'Resolution', P.export_resolution);

exportgraphics(fig, pdf_file, 'ContentType', 'vector');

savefig(fig, fig_file);

%% =========================================================
% 21. Finish
% ==========================================================

fprintf('\n');

fprintf('=========================================================\n');

fprintf(' SINGLE-SEED FTMLE RESULT SAVED\n');

fprintf('=========================================================\n');

fprintf('Noise seed        : %d\n', initial_noise_seed);

fprintf('Preparation scan  : %d ns\n', N_scan_sym);

fprintf('Washout end       : %.0f ns\n', washout_end_ns);

fprintf('Training end      : %.0f ns\n', training_end_ns);

fprintf('FTMLE end         : %.0f ns\n', total_end_ns);

fprintf('Final FTMLE       : %+.8e ns^-1\n\n', FTMLE_final_ns_inv);

fprintf('Summary MAT:\n%s\n\n', summary_file);

fprintf('Prepared-state cache:\n%s\n\n', prepared_dir);

fprintf('FTMLE cache:\n%s\n\n', ftmle_cache_dir);

fprintf('Figure:\n%s\n', png_file);

fprintf('=========================================================\n');

end


%% =========================================================
% RK4IP operator preparation
% ==========================================================

function cache = rk4ip_lle_prepare(delta0, params, dt, OM)

t_R = params.t_R;

L = params.L;

Beta2 = params.Beta2;

gamma = params.gamma;

alpha = params.alpha;

theta = params.theta;

D_lin = -(alpha + 1i*delta0)/t_R + 1i*(L/t_R)*(Beta2/2).*(OM.^2);

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
% One RK4IP step
% ==========================================================

function E_next = rk4ip_lle_step_fast(E, drive, cache)

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

N1 = c_nl .* abs(U0).^2 .* U0 + drive;

K1 = N1;

%% K2

Utemp = U0 + dt_h*K1;

E2 = ifft(fft(Utemp) .* expD_h2);

N2 = c_nl .* abs(E2).^2 .* E2 + drive;

K2 = ifft(fft(N2) .* expMinusD_h2);

%% K3

Utemp = U0 + dt_h*K2;

E3 = ifft(fft(Utemp) .* expD_h2);

N3 = c_nl .* abs(E3).^2 .* E3 + drive;

K3 = ifft(fft(N3) .* expMinusD_h2);

%% K4

Utemp = U0 + dt*K3;

E4 = ifft(fft(Utemp) .* expD_h);

N4 = c_nl .* abs(E4).^2 .* E4 + drive;

K4 = ifft(fft(N4) .* expMinusD_h);

%% RK4 combination

U1 = U0 + dt_6*(K1 + 2*K2 + 2*K3 + K4);

E_next = ifft(fft(U1) .* expD_h);

end