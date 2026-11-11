function kerr_comb_mc_iid_soliton_multirealization_800_5360_800()
% =========================================================
% Kerr-comb linear Memory Capacity (MC)
% SOLITON regime
%
% IID input:
% u(n) ~ U[-1,1]
%
% Data protocol:
% Washout          = 800
% Training history = 200
% Training         = 5360
% Testing history  = 200
% Testing          = 800
% D_max            = 200
%
% Total driven symbols = 7360
%
% Soliton preparation:
% 600-symbol detuning scan
% 400-symbol hold
%
% Parallel computing:
% parfor over independent random seeds
%
% Checkpoint:
% Save every 200 symbols
% Resume automatically after interruption
%
% Cache:
% 1. Prepared soliton state
% 2. Per-seed reservoir states X
% 3. Per-seed MC metrics
% 4. Aggregate results
%
% MATLAB R2020a compatible
% =========================================================

clc;
close all;

%% =========================================================
% 1. USER CONFIGURATION
% =========================================================

P.font_name = 'Times New Roman';

P.tick_font_size = 14;
P.label_font_size = 16;
P.title_font_size = 16;
P.legend_font_size = 14;

P.line_width = 1.8;
P.axis_line_width = 1.0;
P.marker_size = 7;

P.export_resolution = 600;

P.delay_xlabel = 'Delay d';
P.delay_ylabel = 'Memory capacity, MC_d';
P.delay_title = 'Soliton regime: delay-dependent memory capacity';

P.cumulative_xlabel = 'Delay d';
P.cumulative_ylabel = 'Cumulative MC';
P.cumulative_title = 'Soliton regime: cumulative memory capacity';

P.total_xlabel = 'Random realization';
P.total_ylabel = 'Total memory capacity';
P.total_title = 'Soliton regime: total memory capacity';

%% =========================================================
% 2. CACHE CONTROL
% =========================================================

% Normally keep all three options false.

% true: rerun complete reservoir evolution
force_resimulate = false;

% true: recompute Ridge regression and MC
force_recompute_mc = false;

% true: regenerate prepared soliton
force_reprepare_soliton = false;

%% =========================================================
% 3. PARALLEL AND CHECKPOINT SETTINGS
% =========================================================

% Enable parallel computing
use_parfor = true;

% Requested worker count
requested_workers = 4;

% Checkpoint interval
checkpoint_symbols = 200;

% Remove redundant chunks after saving complete X
remove_chunks_after_state_saved = true;

%% =========================================================
% 4. MC EVALUATION PROTOCOL
% =========================================================

N_washout = 800;

D_max = 200;

N_train_history = D_max;

N_train = 5360;

N_test_history = D_max;

N_test = 800;

N_steps = sum([N_washout, N_train_history, N_train, N_test_history, N_test]);

% Ridge regularization
% Modified from 1e-9 to 1e-8
lambda_mc = 1e-8;

% Number of independent realizations
N_realizations = 10;

% Same random seeds as near-critical regime
seed_list = 2026 + (0:N_realizations-1);

% Fixed intracavity initial-noise seed
initial_field_seed = 1;

%% =========================================================
% 5. SAFETY CHECKS
% =========================================================

if N_train_history < D_max
    error('N_train_history must be >= D_max.');
end

if N_test_history < D_max
    error('N_test_history must be >= D_max.');
end

if checkpoint_symbols < 1 || checkpoint_symbols ~= floor(checkpoint_symbols)
    error('checkpoint_symbols must be a positive integer.');
end

%% =========================================================
% 6. PUMP PARAMETERS
% =========================================================

Pav = 0.325;

Delta = 0.05;

%% =========================================================
% 7. KERR-MRR AND LLE PARAMETERS
% =========================================================

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

%% =========================================================
% 8. SYMBOL RATE AND RK4IP TIME STEP
% =========================================================

symbol_rate = 1.0e9;

T_sym = 1/symbol_rate;

N_substeps_sym = 1000;

dt_slow = T_sym/N_substeps_sym;

params.t_R = t_R;
params.L = L;
params.Beta2 = Beta_2;
params.gamma = gamma;
params.alpha = alpha;
params.theta = theta;

%% =========================================================
% 9. SOLITON PREPARATION
% =========================================================

N_scan_sym = 600;

N_hold_sym = 400;

delta_start = -2e-2;

Ein_bias = sqrt(Pav);

%% =========================================================
% 10. OUTPUT FOLDERS
% =========================================================

run_tag = 'iid_soliton_800_5360_800_D200';

result_dir = fullfile(pwd, 'compare_results_mc_iid_800_5360_800_D200');

cache_dir = fullfile(result_dir, 'soliton_cache');

checkpoint_dir = fullfile(cache_dir, 'checkpoints');

state_cache_dir = fullfile(cache_dir, 'reservoir_states');

metric_cache_dir = fullfile(cache_dir, 'mc_metrics');

figure_dir = fullfile(result_dir, 'figures');

make_dir(result_dir);
make_dir(cache_dir);
make_dir(checkpoint_dir);
make_dir(state_cache_dir);
make_dir(metric_cache_dir);
make_dir(figure_dir);

%% =========================================================
% 11. PREPARED-STATE CONFIGURATION
% =========================================================

prepared_cfg = struct();

prepared_cfg.regime = 'soliton';
prepared_cfg.initial_field_seed = initial_field_seed;
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
prepared_cfg.N_hold_sym = N_hold_sym;
prepared_cfg.delta_start = delta_start;

% Additional physical parameters for stricter cache validation
prepared_cfg.lam_center = lam_center;
prepared_cfg.neff = neff;

%% =========================================================
% 12. RESERVOIR-STATE CONFIGURATION
% =========================================================

state_cfg = prepared_cfg;

state_cfg.Delta = Delta;
state_cfg.N_steps = N_steps;
state_cfg.N_washout = N_washout;
state_cfg.N_train_history = N_train_history;
state_cfg.N_train = N_train;
state_cfg.N_test_history = N_test_history;
state_cfg.N_test = N_test;
state_cfg.D_max = D_max;
state_cfg.checkpoint_symbols = checkpoint_symbols;

%% =========================================================
% 13. MC METRIC CONFIGURATION
% =========================================================

metric_cfg = struct();

metric_cfg.D_max = D_max;
metric_cfg.lambda_mc = lambda_mc;
metric_cfg.N_washout = N_washout;
metric_cfg.N_train_history = N_train_history;
metric_cfg.N_train = N_train;
metric_cfg.N_test_history = N_test_history;
metric_cfg.N_test = N_test;

%% =========================================================
% 14. TRAINING AND TESTING INDICES
% =========================================================

idx_train_start = N_washout + N_train_history + 1;

idx_train_state = idx_train_start:(idx_train_start+N_train-1);

idx_test_history_start = idx_train_state(end)+1;

idx_test_history = idx_test_history_start:(idx_test_history_start+N_test_history-1);

idx_test_start = idx_test_history(end)+1;

idx_test_state = idx_test_start:(idx_test_start+N_test-1);

delays = 1:D_max;

train_target_idx = bsxfun(@minus, idx_train_state(:), delays);

test_target_idx = bsxfun(@minus, idx_test_state(:), delays);

%% =========================================================
% 15. INDEX VALIDATION
% =========================================================

earliest_train_target = min(train_target_idx(:));

earliest_test_target = min(test_target_idx(:));

if earliest_train_target < N_washout+1
    error('Training delayed target entered washout.');
end

if earliest_test_target < idx_train_state(end)+1
    error('Testing delayed target entered training.');
end

if idx_test_state(end) ~= N_steps
    error('Incorrect training/testing indexing.');
end

%% =========================================================
% 16. PRINT SIMULATION PROTOCOL
% =========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' SOLITON MC - MULTIPLE IID REALIZATIONS\n');
fprintf('=========================================================\n');

fprintf('Input distribution        : IID U[-1,1]\n');
fprintf('Random realizations       : %d\n', N_realizations);
fprintf('Seed range                : %d - %d\n', seed_list(1), seed_list(end));

fprintf('\n');

fprintf('Washout                   : %d\n', N_washout);
fprintf('Training history          : %d\n', N_train_history);
fprintf('Training                  : %d\n', N_train);
fprintf('Testing history           : %d\n', N_test_history);
fprintf('Testing                   : %d\n', N_test);
fprintf('Driven symbols per run    : %d\n', N_steps);

fprintf('\n');

fprintf('D_max                     : %d\n', D_max);
fprintf('lambda_MC                 : %.3e\n', lambda_mc);
fprintf('Spectral features         : %d\n', N);
fprintf('Symbol rate               : %.3f GSa/s\n', symbol_rate*1e-9);
fprintf('RK4IP slow-time step      : %.3f ps\n', dt_slow*1e12);
fprintf('Checkpoint interval       : %d symbols\n', checkpoint_symbols);

fprintf('\n');

fprintf('Soliton scan              : %d symbols\n', N_scan_sym);
fprintf('Soliton hold              : %d symbols\n', N_hold_sym);

fprintf('\n');

fprintf('Washout indices           : 1 : %d\n', N_washout);
fprintf('Train history indices     : %d : %d\n', N_washout+1, idx_train_start-1);
fprintf('Training indices          : %d : %d\n', idx_train_state(1), idx_train_state(end));
fprintf('Test history indices      : %d : %d\n', idx_test_history(1), idx_test_history(end));
fprintf('Testing indices           : %d : %d\n', idx_test_state(1), idx_test_state(end));

fprintf('\n');

fprintf('Output directory:\n%s\n', result_dir);
fprintf('=========================================================\n');

%% =========================================================
% 17. PREPARED SOLITON CACHE
% =========================================================

prepared_file = fullfile(cache_dir, 'prepared_soliton_state.mat');

legacy_prepared_file = fullfile(pwd, 'compare_results_mc_iid', 'soliton_cache', 'prepared_soliton_state.mat');

E_prepared = [];

if ~force_reprepare_soliton

    E_prepared = load_prepared_if_valid(prepared_file, prepared_cfg, N, false);

    if isempty(E_prepared)

        % Reuse the original soliton preparation when compatible.
        % Original code did not include lam_center and neff
        % in prepared_cfg, so legacy compatibility is checked.

        E_prepared = load_prepared_if_valid(legacy_prepared_file, prepared_cfg, N, true);

        if ~isempty(E_prepared)

            prepared_data = struct();
            prepared_data.E_prepared = E_prepared;
            prepared_data.prepared_cfg = prepared_cfg;

            atomic_save_mat(prepared_file, prepared_data, true);

            fprintf('Compatible original soliton preparation reused.\n');

        end

    end

end

%% =========================================================
% 18. PREPARE SOLITON IF NO VALID CACHE
% =========================================================

if isempty(E_prepared)

    fprintf('\n');
    fprintf('No compatible prepared soliton cache found.\n');
    fprintf('Preparing soliton state...\n');

    % -----------------------------------------------------
    % Initial intracavity field
    % -----------------------------------------------------

    rng(initial_field_seed, 'twister');

    E = 1e-6*(randn(1,N)+1i*randn(1,N));

    % -----------------------------------------------------
    % Stage 1: detuning scan
    % -----------------------------------------------------

    fprintf('Stage 1: detuning scan (%d symbols)\n', N_scan_sym);

    for n = 1:N_scan_sym

        frac = (n-1)/max(N_scan_sym-1,1);

        delta0 = delta_start+(delta0_target-delta_start)*frac;

        op_scan = prepare_rk4ip_operator(delta0, params, dt_slow, OM);

        for kk = 1:N_substeps_sym

            E = rk4ip_step_precomputed(E, Ein_bias, op_scan);

        end

        if mod(n,20)==0 || n==N_scan_sym

            fprintf('Soliton scan: %d / %d symbols\n', n, N_scan_sym);

        end

    end

    fprintf('Detuning scan completed.\n');

    % -----------------------------------------------------
    % Stage 2: soliton hold
    % -----------------------------------------------------

    op_target = prepare_rk4ip_operator(delta0_target, params, dt_slow, OM);

    fprintf('Stage 2: soliton hold (%d symbols)\n', N_hold_sym);

    for n = 1:N_hold_sym

        for kk = 1:N_substeps_sym

            E = rk4ip_step_precomputed(E, Ein_bias, op_target);

        end

        if mod(n,20)==0 || n==N_hold_sym

            fprintf('Soliton hold: %d / %d symbols\n', n, N_hold_sym);

        end

    end

    fprintf('Soliton hold completed.\n');

    E_prepared = E;

    % -----------------------------------------------------
    % Save prepared soliton state
    % -----------------------------------------------------

    prepared_data = struct();

    prepared_data.E_prepared = E_prepared;
    prepared_data.prepared_cfg = prepared_cfg;

    atomic_save_mat(prepared_file, prepared_data, true);

    fprintf('Prepared soliton state saved:\n%s\n', prepared_file);

else

    fprintf('Compatible prepared soliton state loaded.\n');

end

%% =========================================================
% 19. FIXED RK4IP OPERATOR FOR DATA STAGE
% =========================================================

op_target = prepare_rk4ip_operator(delta0_target, params, dt_slow, OM);

%% =========================================================
% 20. PARALLEL POOL
% =========================================================

parallel_ready = false;

if use_parfor

    try

        pool = gcp('nocreate');

        if isempty(pool)

            pool = parpool('local', requested_workers);

        end

        parallel_ready = true;

        fprintf('Parallel pool active: %d workers.\n', pool.NumWorkers);

    catch ME

        fprintf('Parallel pool unavailable: %s\n', ME.message);
        fprintf('Continuing with serial execution.\n');

    end

end

%% =========================================================
% 21. BUILD PER-SEED COMPUTATION CONFIGURATION
% =========================================================

C = struct();

C.N_steps = N_steps;
C.N = N;
C.N_substeps_sym = N_substeps_sym;

C.Pav = Pav;
C.Delta = Delta;

C.state_cfg = state_cfg;
C.metric_cfg = metric_cfg;

C.idx_train_state = idx_train_state;
C.idx_test_state = idx_test_state;

C.train_target_idx = train_target_idx;
C.test_target_idx = test_target_idx;

C.lambda_mc = lambda_mc;
C.D_max = D_max;

C.checkpoint_symbols = checkpoint_symbols;

C.checkpoint_dir = checkpoint_dir;
C.state_cache_dir = state_cache_dir;
C.metric_cache_dir = metric_cache_dir;

C.force_resimulate = force_resimulate;
C.force_recompute_mc = force_recompute_mc;
C.remove_chunks_after_state_saved = remove_chunks_after_state_saved;

%% =========================================================
% 22. MULTIPLE RANDOM REALIZATIONS
% =========================================================

result_cell = cell(N_realizations,1);

if parallel_ready

    parfor r = 1:N_realizations

        result_cell{r} = run_one_seed(seed_list(r), C, E_prepared, op_target);

    end

else

    for r = 1:N_realizations

        result_cell{r} = run_one_seed(seed_list(r), C, E_prepared, op_target);

    end

end

%% =========================================================
% 23. COLLECT RESULTS
% =========================================================

MC_d_all = zeros(D_max,N_realizations);

MC_cumulative_all = zeros(D_max,N_realizations);

MC_total_all = zeros(N_realizations,1);

peak_MC_all = zeros(N_realizations,1);

peak_delay_all = zeros(N_realizations,1);

for r = 1:N_realizations

    metric = result_cell{r};

    MC_d_all(:,r) = metric.MC_d(:);

    MC_cumulative_all(:,r) = metric.MC_cumulative(:);

    MC_total_all(r) = metric.MC_total;

    peak_MC_all(r) = metric.peak_MC;

    peak_delay_all(r) = metric.peak_delay;

end

%% =========================================================
% 24. STATISTICS
% =========================================================

MC_d_mean = mean(MC_d_all,2);

MC_d_std = std(MC_d_all,0,2);

MC_cumulative_mean = mean(MC_cumulative_all,2);

MC_cumulative_std = std(MC_cumulative_all,0,2);

MC_total_mean = mean(MC_total_all);

MC_total_std = std(MC_total_all,0);

%% =========================================================
% 25. PRINT FINAL RESULTS
% =========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' SOLITON MEMORY CAPACITY STATISTICS\n');
fprintf('=========================================================\n');

for r = 1:N_realizations

    fprintf('Run %2d | seed = %d | MC_total = %.8f\n', r, seed_list(r), MC_total_all(r));

end

fprintf('\n');

fprintf('MC_total mean = %.8f\n', MC_total_mean);
fprintf('MC_total std  = %.8f\n', MC_total_std);

fprintf('\n');

fprintf('FINAL RESULT:\n');
fprintf('MC_total = %.4f +/- %.4f\n', MC_total_mean, MC_total_std);

fprintf('\n');

fprintf('Mean cumulative MC at d = 50  : %.8f\n', MC_cumulative_mean(min(50,D_max)));
fprintf('Mean cumulative MC at d = 100 : %.8f\n', MC_cumulative_mean(min(100,D_max)));
fprintf('Mean cumulative MC at d = 200 : %.8f\n', MC_cumulative_mean(min(200,D_max)));

fprintf('=========================================================\n');

%% =========================================================
% 26. SAVE AGGREGATE DATA
% =========================================================

summary_mat = fullfile(result_dir,[run_tag '_summary.mat']);

summary_data = struct();

summary_data.run_tag = run_tag;
summary_data.seed_list = seed_list;
summary_data.N_realizations = N_realizations;

summary_data.symbol_rate = symbol_rate;

summary_data.N_washout = N_washout;
summary_data.N_train_history = N_train_history;
summary_data.N_train = N_train;
summary_data.N_test_history = N_test_history;
summary_data.N_test = N_test;

summary_data.N_steps = N_steps;
summary_data.N = N;
summary_data.D_max = D_max;
summary_data.lambda_mc = lambda_mc;

summary_data.idx_train_state = idx_train_state;
summary_data.idx_test_history = idx_test_history;
summary_data.idx_test_state = idx_test_state;

summary_data.MC_d_all = MC_d_all;
summary_data.MC_cumulative_all = MC_cumulative_all;
summary_data.MC_total_all = MC_total_all;

summary_data.peak_MC_all = peak_MC_all;
summary_data.peak_delay_all = peak_delay_all;

summary_data.MC_d_mean = MC_d_mean;
summary_data.MC_d_std = MC_d_std;

summary_data.MC_cumulative_mean = MC_cumulative_mean;
summary_data.MC_cumulative_std = MC_cumulative_std;

summary_data.MC_total_mean = MC_total_mean;
summary_data.MC_total_std = MC_total_std;

summary_data.prepared_cfg = prepared_cfg;
summary_data.state_cfg = state_cfg;
summary_data.metric_cfg = metric_cfg;

summary_data.P = P;

atomic_save_mat(summary_mat, summary_data, true);

%% =========================================================
% 27. SAVE PER-REALIZATION CSV
% =========================================================

Realization = (1:N_realizations).';

Seed = seed_list(:);

T_runs = table(Realization, Seed, MC_total_all, peak_MC_all, peak_delay_all, 'VariableNames', {'Realization','Seed','MC_Total','Peak_MC_d','Peak_Delay'});

runs_csv = fullfile(result_dir,[run_tag '_runs.csv']);

writetable(T_runs,runs_csv);

%% =========================================================
% 28. SAVE DELAY STATISTICS CSV
% =========================================================

Delay = (1:D_max).';

T_delay = table(Delay, MC_d_mean, MC_d_std, MC_cumulative_mean, MC_cumulative_std, 'VariableNames', {'Delay','MC_Mean','MC_Std','Cumulative_MC_Mean','Cumulative_MC_Std'});

delay_csv = fullfile(result_dir,[run_tag '_delay_statistics.csv']);

writetable(T_delay,delay_csv);

%% =========================================================
% 29. FIGURE 1: DELAY-DEPENDENT MC
% =========================================================

x = (1:D_max).';

fig1_base = fullfile(figure_dir,[run_tag '_MC_delay_mean_std']);

plot_mc_band(x, MC_d_mean, MC_d_std, P.delay_xlabel, P.delay_ylabel, P.delay_title, P, fig1_base);

%% =========================================================
% 30. FIGURE 2: CUMULATIVE MC
% =========================================================

fig2_base = fullfile(figure_dir,[run_tag '_MC_cumulative_mean_std']);

plot_mc_band(x, MC_cumulative_mean, MC_cumulative_std, P.cumulative_xlabel, P.cumulative_ylabel, P.cumulative_title, P, fig2_base);

%% =========================================================
% 31. FIGURE 3: MC REPRODUCIBILITY
% =========================================================

fig3 = figure('Name','Soliton MC total reproducibility','Color','w');

hold on;

plot(1:N_realizations, MC_total_all, 'o-', 'MarkerSize', P.marker_size, 'LineWidth', 1.2);

yline(MC_total_mean, '-', sprintf('Mean = %.3f', MC_total_mean), 'LineWidth', 1.5);

yline(MC_total_mean+MC_total_std, '--', '+1 SD', 'LineWidth', 1);

yline(MC_total_mean-MC_total_std, '--', '-1 SD', 'LineWidth', 1);

xlabel(P.total_xlabel, 'FontName', P.font_name, 'FontSize', P.label_font_size);

ylabel(P.total_ylabel, 'FontName', P.font_name, 'FontSize', P.label_font_size);

title(sprintf('%s: %.3f \\pm %.3f', P.total_title, MC_total_mean, MC_total_std), 'FontName', P.font_name, 'FontSize', P.title_font_size, 'FontWeight', 'normal');

set(gca, 'FontName', P.font_name, 'FontSize', P.tick_font_size, 'LineWidth', P.axis_line_width);

xlim([0.5 N_realizations+0.5]);

xticks(1:N_realizations);

grid on;
box on;

fig3_base = fullfile(figure_dir,[run_tag '_MC_total_realizations']);

exportgraphics(fig3,[fig3_base '.png'],'Resolution',P.export_resolution);

savefig(fig3,[fig3_base '.fig']);

%% =========================================================
% 32. FINAL INFORMATION
% =========================================================

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' ALL SOLITON MC RESULTS SAVED\n');
fprintf('=========================================================\n');

fprintf('Summary MAT:\n%s\n\n',summary_mat);
fprintf('Per-run CSV:\n%s\n\n',runs_csv);
fprintf('Delay statistics CSV:\n%s\n\n',delay_csv);

fprintf('Figures:\n%s\n\n',figure_dir);
fprintf('Reservoir-state cache:\n%s\n\n',state_cache_dir);
fprintf('Checkpoint cache:\n%s\n\n',checkpoint_dir);
fprintf('MC metric cache:\n%s\n',metric_cache_dir);

fprintf('=========================================================\n');

end

%% =========================================================
% SUBFUNCTION 1
% RUN ONE RANDOM REALIZATION
% =========================================================

function metric = run_one_seed(seed, C, E_prepared, op_target)

state_file = fullfile(C.state_cache_dir,sprintf('soliton_state_seed_%d.mat',seed));

metric_file = fullfile(C.metric_cache_dir,sprintf('soliton_mc_seed_%d.mat',seed));

seed_checkpoint_dir = fullfile(C.checkpoint_dir,sprintf('seed_%d',seed));

make_dir(seed_checkpoint_dir);

%% ---------------------------------------------------------
% 1. Check completed MC cache
% ---------------------------------------------------------

if ~C.force_resimulate && ~C.force_recompute_mc

    metric = load_metric_if_valid(metric_file,C.state_cfg,C.metric_cfg,C.D_max);

    if ~isempty(metric)

        fprintf('Seed %d: completed MC cache loaded.\n',seed);

        return;

    end

end

%% ---------------------------------------------------------
% 2. Check complete reservoir-state cache
% ---------------------------------------------------------

X = [];
u_mc = [];

if ~C.force_resimulate

    [X,u_mc] = load_state_if_valid(state_file,C.state_cfg,C.N_steps,C.N);

end

%% ---------------------------------------------------------
% 3. Generate or resume reservoir states
% ---------------------------------------------------------

if isempty(X)

    % -----------------------------------------------------
    % Generate deterministic IID input
    % -----------------------------------------------------

    rs_mc = RandStream('mt19937ar','Seed',seed);

    u_mc = 2*rand(rs_mc,C.N_steps,1)-1;

    Pin_series = max(C.Pav+(C.Delta/2)*u_mc,0);

    Ein_series = sqrt(Pin_series);

    % -----------------------------------------------------
    % Initialize reservoir
    % -----------------------------------------------------

    X = zeros(C.N_steps,C.N);

    E = E_prepared;

    n_completed = 0;

    num_chunks = ceil(C.N_steps/C.checkpoint_symbols);

    %% -----------------------------------------------------
    % 4. Recover completed checkpoints
    % -----------------------------------------------------

    if ~C.force_resimulate

        for block = 1:num_chunks

            chunk_file = fullfile(seed_checkpoint_dir,sprintf('chunk_%04d.mat',block));

            start_n = (block-1)*C.checkpoint_symbols+1;

            end_n = min(block*C.checkpoint_symbols,C.N_steps);

            [ok,chunk] = load_chunk_if_valid(chunk_file,seed,start_n,end_n,C.state_cfg,C.N);

            if ~ok

                delete_chunks_from(seed_checkpoint_dir,block,num_chunks);

                break;

            end

            X(start_n:end_n,:) = chunk.X_chunk;

            E = chunk.E_end;

            n_completed = end_n;

        end

    else

        delete_chunks_from(seed_checkpoint_dir,1,num_chunks);

    end

    fprintf('Seed %d: starting from symbol %d / %d.\n',seed,n_completed+1,C.N_steps);

    %% -----------------------------------------------------
    % 5. Continue RK4IP evolution
    % -----------------------------------------------------

    first_new_block = ceil(n_completed/C.checkpoint_symbols)+1;

    for block = first_new_block:num_chunks

        start_n = (block-1)*C.checkpoint_symbols+1;

        end_n = min(block*C.checkpoint_symbols,C.N_steps);

        X_chunk = zeros(end_n-start_n+1,C.N);

        for n = start_n:end_n

            Ein_n = Ein_series(n);

            for kk = 1:C.N_substeps_sym

                E = rk4ip_step_precomputed(E,Ein_n,op_target);

            end

            F_spec = fftshift(fft(E));

            X_chunk(n-start_n+1,:) = abs(F_spec).^2;

        end

        %% -------------------------------------------------
        % 6. Check numerical validity
        % -------------------------------------------------

        if ~all(isfinite(X_chunk(:))) || ~all(isfinite(E(:)))

            error('Nonfinite reservoir state in seed %d, chunk %d.',seed,block);

        end

        X(start_n:end_n,:) = X_chunk;

        %% -------------------------------------------------
        % 7. Save checkpoint
        % -------------------------------------------------

        chunk = struct();

        chunk.seed = seed;
        chunk.start_n = start_n;
        chunk.end_n = end_n;

        chunk.state_cfg = C.state_cfg;

        chunk.X_chunk = X_chunk;
        chunk.E_end = E;

        chunk_file = fullfile(seed_checkpoint_dir,sprintf('chunk_%04d.mat',block));

        checkpoint_data = struct();

        checkpoint_data.chunk_data = chunk;

        atomic_save_mat(chunk_file,checkpoint_data,false);

        fprintf('Seed %d: checkpoint %d / %d saved (%d / %d symbols).\n',seed,block,num_chunks,end_n,C.N_steps);

    end

    %% -----------------------------------------------------
    % 8. Save complete reservoir states
    % -----------------------------------------------------

    full_state = struct();

    full_state.X = X;
    full_state.u_mc = u_mc;
    full_state.seed = seed;
    full_state.state_cfg = C.state_cfg;

    atomic_save_mat(state_file,full_state,true);

    fprintf('Seed %d: complete reservoir states saved.\n',seed);

    %% -----------------------------------------------------
    % 9. Remove redundant checkpoint files
    % -----------------------------------------------------

    if C.remove_chunks_after_state_saved

        delete_chunks_from(seed_checkpoint_dir,1,num_chunks);

    end

else

    fprintf('Seed %d: reusing complete reservoir states.\n',seed);

end

%% =========================================================
% 10. TRAINING / TESTING STATES
% =========================================================

X_train = X(C.idx_train_state,:);

X_test = X(C.idx_test_state,:);

%% =========================================================
% 11. FEATURE NORMALIZATION
% =========================================================

[X_train_n,mu_X,sigma_X] = zscore_features(X_train);

X_test_n = bsxfun(@rdivide,bsxfun(@minus,X_test,mu_X),sigma_X);

%% =========================================================
% 12. BIAS
% =========================================================

Xtr = [X_train_n,ones(numel(C.idx_train_state),1)];

Xte = [X_test_n,ones(numel(C.idx_test_state),1)];

%% =========================================================
% 13. RIDGE REGULARIZATION
% =========================================================

Ireg = eye(size(Xtr,2));

% Bias is not regularized
Ireg(end,end) = 0;

%% =========================================================
% 14. DELAYED TARGETS
% =========================================================

Y_train = u_mc(C.train_target_idx);

Y_test = u_mc(C.test_target_idx);

%% =========================================================
% 15. RIDGE REGRESSION
% =========================================================

A = Xtr.'*Xtr+C.lambda_mc*Ireg;

B = A\(Xtr.'*Y_train);

Y_pred = Xte*B;

%% =========================================================
% 16. DELAY-DEPENDENT MC
% =========================================================

Yt = bsxfun(@minus,Y_test,mean(Y_test,1));

Yp = bsxfun(@minus,Y_pred,mean(Y_pred,1));

numerator = sum(Yt.*Yp,1).^2;

denominator = sum(Yt.^2,1).*sum(Yp.^2,1);

MC_d = zeros(C.D_max,1);

valid = denominator>0;

MC_d(valid) = numerator(valid)./denominator(valid);

%% =========================================================
% 17. TOTAL AND CUMULATIVE MC
% =========================================================

MC_cumulative = cumsum(MC_d);

MC_total = sum(MC_d);

[peak_MC,peak_delay] = max(MC_d);

%% =========================================================
% 18. SAVE MC METRIC
% =========================================================

metric = struct();

metric.seed = seed;
metric.MC_d = MC_d;
metric.MC_cumulative = MC_cumulative;
metric.MC_total = MC_total;

metric.peak_MC = peak_MC;
metric.peak_delay = peak_delay;

metric.metric_cfg = C.metric_cfg;
metric.state_cfg = C.state_cfg;

atomic_save_mat(metric_file,metric,false);

fprintf('Seed %d finished. MC_total = %.8f\n',seed,MC_total);

end

%% =========================================================
% SUBFUNCTION 2
% LOAD PREPARED SOLITON CACHE
% =========================================================

function E_prepared = load_prepared_if_valid(file_path,prepared_cfg,N,allow_legacy)

E_prepared = [];

if exist(file_path,'file') ~= 2
    return;
end

try

    tmp = load(file_path,'prepared_cfg','E_prepared');

    if ~isfield(tmp,'prepared_cfg') || ~isfield(tmp,'E_prepared')
        return;
    end

    compatible = isequaln(tmp.prepared_cfg,prepared_cfg);

    % -----------------------------------------------------
    % Compatibility with original preparation cache
    % -----------------------------------------------------

    if ~compatible && allow_legacy

        legacy_cfg = rmfield(prepared_cfg,{'lam_center','neff'});

        % The original program used these two fixed values.
        if prepared_cfg.lam_center==1.600e-6 && prepared_cfg.neff==1.444

            compatible = isequaln(tmp.prepared_cfg,legacy_cfg);

        end

    end

    if ~compatible
        return;
    end

    if numel(tmp.E_prepared) ~= N
        return;
    end

    if ~all(isfinite(tmp.E_prepared(:)))
        return;
    end

    E_prepared = reshape(tmp.E_prepared,1,[]);

catch

    E_prepared = [];

end

end

%% =========================================================
% SUBFUNCTION 3
% LOAD MC METRIC CACHE
% =========================================================

function metric = load_metric_if_valid(file_path,state_cfg,metric_cfg,D_max)

metric = [];

if exist(file_path,'file') ~= 2
    return;
end

try

    tmp = load(file_path);

    required = {'state_cfg','metric_cfg','MC_d','MC_cumulative','MC_total','peak_MC','peak_delay'};

    if ~all(isfield(tmp,required))
        return;
    end

    if ~isequaln(tmp.state_cfg,state_cfg)
        return;
    end

    if ~isequaln(tmp.metric_cfg,metric_cfg)
        return;
    end

    if numel(tmp.MC_d) ~= D_max
        return;
    end

    if numel(tmp.MC_cumulative) ~= D_max
        return;
    end

    if ~isfinite(tmp.MC_total)
        return;
    end

    if ~all(isfinite(tmp.MC_d(:)))
        return;
    end

    if ~all(isfinite(tmp.MC_cumulative(:)))
        return;
    end

    metric = tmp;

catch

    metric = [];

end

end

%% =========================================================
% SUBFUNCTION 4
% LOAD COMPLETE RESERVOIR STATES
% =========================================================

function [X,u_mc] = load_state_if_valid(file_path,state_cfg,N_steps,N)

X = [];
u_mc = [];

if exist(file_path,'file') ~= 2
    return;
end

try

    cfg = load(file_path,'state_cfg');

    if ~isfield(cfg,'state_cfg')
        return;
    end

    if ~isequaln(cfg.state_cfg,state_cfg)
        return;
    end

    tmp = load(file_path,'X','u_mc');

    if ~isfield(tmp,'X') || ~isfield(tmp,'u_mc')
        return;
    end

    if ~isequal(size(tmp.X),[N_steps,N])
        return;
    end

    if numel(tmp.u_mc) ~= N_steps
        return;
    end

    if ~all(isfinite(tmp.X(:)))
        return;
    end

    if ~all(isfinite(tmp.u_mc(:)))
        return;
    end

    X = tmp.X;

    u_mc = tmp.u_mc(:);

catch

    X = [];
    u_mc = [];

end

end

%% =========================================================
% SUBFUNCTION 5
% LOAD CHECKPOINT
% =========================================================

function [ok,chunk] = load_chunk_if_valid(file_path,seed,start_n,end_n,state_cfg,N)

ok = false;

chunk = [];

if exist(file_path,'file') ~= 2
    return;
end

try

    tmp = load(file_path,'chunk_data');

    if ~isfield(tmp,'chunk_data')
        return;
    end

    chunk = tmp.chunk_data;

    required = {'seed','start_n','end_n','state_cfg','X_chunk','E_end'};

    if ~all(isfield(chunk,required))
        return;
    end

    if chunk.seed ~= seed
        return;
    end

    if chunk.start_n ~= start_n
        return;
    end

    if chunk.end_n ~= end_n
        return;
    end

    if ~isequaln(chunk.state_cfg,state_cfg)
        return;
    end

    if ~isequal(size(chunk.X_chunk),[end_n-start_n+1,N])
        return;
    end

    if numel(chunk.E_end) ~= N
        return;
    end

    if ~all(isfinite(chunk.X_chunk(:)))
        return;
    end

    if ~all(isfinite(chunk.E_end(:)))
        return;
    end

    chunk.E_end = reshape(chunk.E_end,1,[]);

    ok = true;

catch

    ok = false;
    chunk = [];

end

end

%% =========================================================
% SUBFUNCTION 6
% DELETE REDUNDANT CHECKPOINTS
% =========================================================

function delete_chunks_from(folder,first_block,last_block)

for block = first_block:last_block

    file_path = fullfile(folder,sprintf('chunk_%04d.mat',block));

    if exist(file_path,'file') == 2

        delete(file_path);

    end

end

end

%% =========================================================
% SUBFUNCTION 7
% SAFE MAT SAVE
% =========================================================

function atomic_save_mat(file_path,data,use_v73)

% Save to a temporary file first.
% A partially written temporary file is ignored on restart.

temp_path = [file_path '.writing.mat'];

if use_v73

    save(temp_path,'-struct','data','-v7.3');

else

    save(temp_path,'-struct','data');

end

[ok,message] = movefile(temp_path,file_path,'f');

if ~ok

    error('Unable to finalize cache %s: %s',file_path,message);

end

end

%% =========================================================
% SUBFUNCTION 8
% CREATE DIRECTORY
% =========================================================

function make_dir(folder)

if exist(folder,'dir') ~= 7

    [ok,msg] = mkdir(folder);

    if ~ok

        error('Cannot create directory %s: %s',folder,msg);

    end

end

end

%% =========================================================
% SUBFUNCTION 9
% PLOT MC MEAN AND STANDARD DEVIATION
% =========================================================

function plot_mc_band(x,y_mean,y_std,x_name,y_name,plot_title,P,file_base)

fig = figure('Name',plot_title,'Color','w');

hold on;

colors = get(gca,'ColorOrder');

c1 = colors(1,:);

upper = y_mean+y_std;

lower = max(y_mean-y_std,0);

% ---------------------------------------------------------
% Standard deviation band
% ---------------------------------------------------------

h_band = fill([x;flipud(x)],[upper;flipud(lower)],c1,'FaceAlpha',0.18,'EdgeColor','none');

% ---------------------------------------------------------
% Mean curve
% ---------------------------------------------------------

h_mean = plot(x,y_mean,'-','Color',c1,'LineWidth',P.line_width);

xlabel(x_name,'FontName',P.font_name,'FontSize',P.label_font_size);

ylabel(y_name,'FontName',P.font_name,'FontSize',P.label_font_size);

title(plot_title,'FontName',P.font_name,'FontSize',P.title_font_size,'FontWeight','normal');

lgd = legend([h_mean,h_band],{'Mean','Mean \pm 1 SD'},'Location','best');

set(lgd,'FontName',P.font_name,'FontSize',P.legend_font_size);

set(gca,'FontName',P.font_name,'FontSize',P.tick_font_size,'LineWidth',P.axis_line_width);

xlim([1 numel(x)]);

ymax = max(upper);

if ~isfinite(ymax) || ymax<=0

    ymax = 1;

end

ylim([0 1.05*ymax]);

grid on;
box on;

exportgraphics(fig,[file_base '.png'],'Resolution',P.export_resolution);

savefig(fig,[file_base '.fig']);

end

%% =========================================================
% SUBFUNCTION 10
% PRECOMPUTE RK4IP OPERATOR
% =========================================================

function op = prepare_rk4ip_operator(delta0,params,dt,OM)

t_R = params.t_R;

L = params.L;

Beta2 = params.Beta2;

gamma = params.gamma;

alpha = params.alpha;

theta = params.theta;

D_lin = (-(alpha+1i*delta0)/t_R)+1i*(L/t_R)*(Beta2/2).*(OM.^2);

op.expD_h2 = exp(D_lin*(dt/2));

op.expD_h = exp(D_lin*dt);

op.expMinusD_h2 = exp(-D_lin*(dt/2));

op.expMinusD_h = exp(-D_lin*dt);

op.dt = dt;

op.c_nl = 1i*gamma*L/t_R;

op.drive_scale = sqrt(theta)/t_R;

end

%% =========================================================
% SUBFUNCTION 11
% ONE RK4IP STEP
% =========================================================

function E_next = rk4ip_step_precomputed(E,Ein,op)

dt = op.dt;

c_nl = op.c_nl;

drive = op.drive_scale*Ein;

expD_h2 = op.expD_h2;

expD_h = op.expD_h;

expMinusD_h2 = op.expMinusD_h2;

expMinusD_h = op.expMinusD_h;

apply_expD = @(X,expD) ifft(ifftshift(fftshift(fft(X)).*expD));

U0 = E;

% ---------------------------------------------------------
% K1
% ---------------------------------------------------------

N1 = c_nl.*abs(U0).^2.*U0+drive;

K1 = N1;

% ---------------------------------------------------------
% K2
% ---------------------------------------------------------

Utemp = U0+(dt/2)*K1;

E2 = apply_expD(Utemp,expD_h2);

N2 = c_nl.*abs(E2).^2.*E2+drive;

K2 = apply_expD(N2,expMinusD_h2);

% ---------------------------------------------------------
% K3
% ---------------------------------------------------------

Utemp = U0+(dt/2)*K2;

E3 = apply_expD(Utemp,expD_h2);

N3 = c_nl.*abs(E3).^2.*E3+drive;

K3 = apply_expD(N3,expMinusD_h2);

% ---------------------------------------------------------
% K4
% ---------------------------------------------------------

Utemp = U0+dt*K3;

E4 = apply_expD(Utemp,expD_h);

N4 = c_nl.*abs(E4).^2.*E4+drive;

K4 = apply_expD(N4,expMinusD_h);

% ---------------------------------------------------------
% RK4 update
% ---------------------------------------------------------

U1 = U0+(dt/6)*(K1+2*K2+2*K3+K4);

E_next = apply_expD(U1,expD_h);

end

%% =========================================================
% SUBFUNCTION 12
% FEATURE-WISE Z-SCORE NORMALIZATION
% =========================================================

function [Xn,mu_X,sigma_X] = zscore_features(X)

mu_X = mean(X,1);

sigma_X = std(X,0,1);

sigma_X(sigma_X<1e-12) = 1e-12;

Xn = bsxfun(@rdivide,bsxfun(@minus,X,mu_X),sigma_X);

end