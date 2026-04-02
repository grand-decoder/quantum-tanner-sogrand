% This script runs Monte Carlo simulation to evaluate the performance of 
% SOGRAND decoding of quantum Tanner codes. The simulation parameters can
% be specified below under "simulation parameters".
% Run setup_env.m before running a simulation.

% Simulation results are stored in the folder 
% data/outputs/{yy_MM_dd}/{HH_mm_ss}_osd{osd_post_processing}_corr{correlation_aware_decoding}_{file_name_part}
% containing:
% - log.txt: simulation log file.
% - result.mat:
%   - stores the simulation parameters.
%   - p_list: array containing the p_physical values that were simulated.
%   - bler_deg_list: Logical block error rate (BLER).
%   - bler_deg_per_iteration_list: if calc_performance_per_iter=true, this
%     stores bler_deg_list for each decoding iteration.
%   - bler_list, bler_per_iteration_list: same as bler_deg_list, bler_deg_per_iteration_list
%     but count a decoding as an error if not the exact error pattern is
%     recovered whereas bler_deg_list and bler_deg_per_iteration_list take
%     degeneracy into account.
%   - num_guesses_total_list, num_guesses_per_iteration_list: Number of
%     total guesses that SOGRAND carried out for all component codewords to
%     decode.

% [2] P. Yuan, M. Médard, K. Galligan, and K. R. Duffy, "Soft-output (SO) GRAND and long, low rate codes to outperform 5 LDPCs," IEEE Trans. Wireless Commun., vol. 24, no. 4, pp. 3386–3399, Apr. 2025.
% [3] R. K. Radebold, S. D. Bartlett, and A. C. Doherty, "Explicit instances of quantum Tanner codes," Nov. 14, 2025, arXiv: arXiv:2508.05095. doi: 10.48550/arXiv.2508.05095.
% [8] K. R. Duffy, W. An, and M. Médard, "Ordered reliability bits guessing random additive noise decoding," IEEE Trans. Signal Proc., vol. 70, pp. 4528–4542, Aug. 2022, doi: 10.1109/TSP.2022.3203251.

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters (change parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of the quantum Tanner code. Currently, the quantum tanner codes from
% [3] are supported, which can be downloaded from https://github.com/RebKatRad/qTanner or via setup_data.m.
code_name = "[[250_10_15]]"; % [[250_10_15]] or [[200_10_10]]

% Maximal number of block errors before simulation terminates.
max_num_deg_errors = 10;
% Maximal number of codewords before the simulation terminates.
max_num_simulations = 100;

% Physical error rates p_physical are simulated at p_count log-uniformly
% spaced values over [10^p_start, 10^p_end].
p_start = -3;
p_end = -0.8;
p_count = 15;

% If true, the decoding performance is evaluates for each decoding
% iteration.
calc_performance_per_iter = true;

% Runs the Monte Carlo simulation in parallel on multiple cores.
use_parallel = false;

% If use_parallel = true, simulation_batch_time (in seconds) sets how long
% each worker core simulates before results are accumulated and checked
% against the termination criterion. Smaller values enable earlier
% termination but increase synchronization overhead; larger values reduce
% idle time but may overshoot max_num_deg_errors.
simulation_batch_time = 20;

% SOGRAND list size
decoder_property.list_size = 4;

% SOGRAND terminates list decoding early if the estimated list error 
% probability is below this threshold (see [2] for details).
decoder_property.sogrand_thr = 10^-5;

% Maximal number of decoding iterations
decoder_property.max_iter = 20;

% Scaling factor with which extrinsic messages are scaled during message
% passing.
decoder_property.alpha = 1.0;

% Maximal number of queries that SOGRAND carries out before stopping list
% decoding. If max_query=inf, querying only stops of the list_size is
% reached or the estimated list error probability is below 
% decoder_property.sogrand_thr.
decoder_property.max_query=inf;

% Intercept (integer >= 0) for ORBGRAND noise generator. 
% If IC = -1, the intercept is automatically determined based on the
% current LLRs (see [8] for details).
decoder_property.IC = -1;

% If true, correlation between X and Z error from the depolarizing channel
% are considered.
decoder_property.use_x_z_correlation = true;

% If true, OSD with combination strategy using the LDPCv2 package is used 
% as a post-processing step.
decoder_property.use_osd_decoder = true;
% OSD order for post-processing.
decoder_property.osd_order = 9;

% file prefix to name simulation run.
file_name_part = "test";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation (no changes necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_list = logspace(p_start, p_end, p_count);

file_name = code_name + "_osd" + string(double(decoder_property.use_osd_decoder)) + "_corr" + string(double(decoder_property.use_x_z_correlation)) + "_" + file_name_part;

% Get the root folder relative to this script's location
script_dir = fileparts(mfilename('fullpath'));
root_folder = fullfile(script_dir, '..');

if code_name=="[[250_10_15]]"
    lrcc_path = fullfile(root_folder, "data/codes/[[250_10_15]]/LRCC");
    code_path = fullfile(root_folder, "data/codes/classical_codes/delta_5_for_[[250_10_15]]");
elseif code_name=="[[200_10_10]]"
    lrcc_path = fullfile(root_folder, "data/codes/[[200_10_10]]/LRCC");
    code_path = fullfile(root_folder, "data/codes/classical_codes/delta_5_for_[[200_10_10]]");
else
    error("code not found");
end

% Code construction
quantum_tanner_code=quantum_tanner_code_construction(lrcc_path, code_path);

result_folder_path = create_result_folder(file_name);
log_message(result_folder_path, "Started");

result.git_hash = get_git_hash();
log_message(result_folder_path, "Git hash: " + result.git_hash);

slurm_id = getenv('SLURM_JOB_ID');
if isempty(slurm_id)
    slurm_id = 'N/A_Local_Session';
end
log_message(result_folder_path, "Slurm Job ID: " + slurm_id);

num_decoded_sum_list = zeros(length(p_list),1);
num_error_sum_list = zeros(length(p_list),1);
num_degenereated_error_sum_list = zeros(length(p_list),1);
num_error_per_iteration_sum_list = zeros(length(p_list), decoder_property.max_iter);
num_error_deg_per_iteration_sum_list = zeros(length(p_list), decoder_property.max_iter);
num_guesses_total_sum_list = zeros(length(p_list), 1);
num_guesses_per_iteration_sum_list  = zeros(length(p_list), decoder_property.max_iter);

if use_parallel
    poolObj = gcp('nocreate');
    if isempty(poolObj)
        poolObj = parpool(); % Starts default pool if none exists
    end
    num_workers = poolObj.NumWorkers; % Number of available cores/workers
    result.num_workers = num_workers;
    log_message(result_folder_path, "num_workers: " + string(num_workers));

    parallel_quantum_tanner_code = parallel.pool.Constant(quantum_tanner_code);
end



for ii=1:length(p_list)
    log_message(result_folder_path, "Running simulation: p_idx: " + string(ii) + "/" + string(length(p_list)));
    disp(ii);
    start_time = tic;
    % Noise variance

    p_physical = p_list(ii);
    if ~use_parallel
        osd_decoder_x = []; 
        osd_decoder_z = [];
        if decoder_property.use_osd_decoder
            % x and z are swapped here on purpose!
            osd_decoder_x = OsdDecoderWrapper(quantum_tanner_code.H_global_z, p_physical, decoder_property.osd_order);
            osd_decoder_z = OsdDecoderWrapper(quantum_tanner_code.H_global_x, p_physical, decoder_property.osd_order);
        end
    end
    
    num_decoded_sum = 0;
    num_error_sum = 0;
    num_degenereated_error_sum = 0;
    num_error_per_iteration_sum = zeros(decoder_property.max_iter, 1);
    num_error_deg_per_iteration_sum = zeros(decoder_property.max_iter, 1);
    num_guesses_total_sum = 0;
    num_guesses_per_iteration_sum = zeros(decoder_property.max_iter, 1);
    
    while num_degenereated_error_sum < max_num_deg_errors && num_decoded_sum < max_num_simulations
        if use_parallel
            parfor batch_idx = 1:num_workers
                % the python osd decoder can not be serialized. hence,
                % generate it here inside every time. it does not seem to
                % cause any relevant overhead.
                osd_decoder_x = []; 
                osd_decoder_z = [];
                if decoder_property.use_osd_decoder
                    % x and z are swapped here on purpose!
                    osd_decoder_x = OsdDecoderWrapper(parallel_quantum_tanner_code.Value.H_global_z, p_physical, decoder_property.osd_order);
                    osd_decoder_z = OsdDecoderWrapper(parallel_quantum_tanner_code.Value.H_global_x, p_physical, decoder_property.osd_order);
                end

                [batch_decoded, batch_error, batch_deg_error, batch_error_iter, ...
                 batch_deg_iter, batch_guesses_total, batch_guesses_iter] = run_decoding_batch( ...
                 simulation_batch_time, p_physical, parallel_quantum_tanner_code.Value, decoder_property, max_num_deg_errors, ...
                 max_num_simulations, osd_decoder_x, osd_decoder_z, calc_performance_per_iter);
        
                num_decoded_sum                 = num_decoded_sum + batch_decoded;
                num_error_sum                   = num_error_sum + batch_error;
                num_degenereated_error_sum      = num_degenereated_error_sum + batch_deg_error;
                num_error_per_iteration_sum     = num_error_per_iteration_sum + batch_error_iter;
                num_error_deg_per_iteration_sum = num_error_deg_per_iteration_sum + batch_deg_iter;
                num_guesses_total_sum           = num_guesses_total_sum + batch_guesses_total;
                num_guesses_per_iteration_sum   = num_guesses_per_iteration_sum + batch_guesses_iter;
            end
        else
            [batch_decoded, batch_error, batch_deg_error, batch_error_iter, ...
             batch_deg_iter, batch_guesses_total, batch_guesses_iter] = run_decoding( ...
             p_physical, quantum_tanner_code, decoder_property, osd_decoder_x, ...
             osd_decoder_z, calc_performance_per_iter);
    
            num_decoded_sum                 = num_decoded_sum + batch_decoded;
            num_error_sum                   = num_error_sum + batch_error;
            num_degenereated_error_sum      = num_degenereated_error_sum + batch_deg_error;
            num_error_per_iteration_sum     = num_error_per_iteration_sum + batch_error_iter;
            num_error_deg_per_iteration_sum = num_error_deg_per_iteration_sum + batch_deg_iter;
            num_guesses_total_sum           = num_guesses_total_sum + batch_guesses_total;
            num_guesses_per_iteration_sum   = num_guesses_per_iteration_sum + batch_guesses_iter;
        end
    end

    num_decoded_sum_list(ii)                 = num_decoded_sum;
    num_error_sum_list(ii)                   = num_error_sum;
    num_degenereated_error_sum_list(ii)      = num_degenereated_error_sum;
    num_error_per_iteration_sum_list(ii, :)  = num_error_per_iteration_sum;
    num_error_deg_per_iteration_sum_list(ii, :) = num_error_deg_per_iteration_sum;
    num_guesses_total_sum_list(ii)           = num_guesses_total_sum;
    num_guesses_per_iteration_sum_list(ii, :) = num_guesses_per_iteration_sum;
    
    log_message(result_folder_path, "p_idx: " + string(ii) + " runtime: " + string(toc(start_time)) + " s");
    
    if num_decoded_sum > 0
        bler = num_error_sum / num_decoded_sum;
        bler_deg = num_degenereated_error_sum / num_decoded_sum;
    else
        bler = -1;
        bler_deg = -1;
    end
    sum_details = sprintf('Sums -> Decoded: %d, Errors: %d, DegenErrors: %d | BLER: %.4e, BLER_deg: %.4e', ...
            num_decoded_sum, num_error_sum, num_degenereated_error_sum, bler, bler_deg);
    log_message(result_folder_path, sum_details);
end

log_message(result_folder_path, "Storing Results");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
result.p_list = p_list;
result.p_start = p_start;
result.p_end = p_end;
result.p_count = p_count;

result.decoder_property = decoder_property;
result.lrcc_path = lrcc_path;
result.code_path = code_path;
result.code_name = code_name;

result.max_num_deg_errors = max_num_deg_errors;
result.max_num_simulations = max_num_simulations;
result.simulation_batch_time = simulation_batch_time;
result.calc_performance_per_iter = calc_performance_per_iter;

% Raw data
result.num_decoded_sum_list = num_decoded_sum_list;
result.num_error_sum_list = num_error_sum_list;
result.num_degenereated_error_sum_list = num_degenereated_error_sum_list;
result.num_error_per_iteration_sum_list = num_error_per_iteration_sum_list;
result.num_error_deg_per_iteration_sum_list = num_error_deg_per_iteration_sum_list;

% Normalized Performance data
bler_list = num_error_sum_list ./ num_decoded_sum_list;
bler_deg_list = num_degenereated_error_sum_list ./ num_decoded_sum_list;

bler_per_iteration_list = num_error_per_iteration_sum_list ./ num_decoded_sum_list;
bler_deg_per_iteration_list = num_error_deg_per_iteration_sum_list ./ num_decoded_sum_list;

result.bler_list = bler_list;
result.bler_deg_list = bler_deg_list;

result.bler_per_iteration_list = bler_per_iteration_list;
result.bler_deg_per_iteration_list = bler_deg_per_iteration_list;

result.num_guesses_total_list = num_guesses_total_sum_list  ./ num_decoded_sum_list;
result.num_guesses_per_iteration_list = num_guesses_per_iteration_sum_list  ./ num_decoded_sum_list;

result_file_path = fullfile(result_folder_path, "result.mat");
save(result_file_path, 'result');

log_message(result_folder_path, "Finished");

function [num_decoded, num_error_sum, num_degenereated_error_sum, num_error_per_iteration_sum, ...
    num_error_deg_per_iteration_sum, num_guesses_total_sum, num_guesses_per_iteration_sum]=run_decoding_batch( ...
    simulation_batch_time, p_physical, quantum_tanner_code, decoder_property, max_num_deg_errors, max_num_simulations, osd_decoder_x, osd_decoder_z, calc_performance_per_iter)
% Runs repeated decoding trials for a fixed duration (used in parallel mode).
% Calls run_decoding in a loop until simulation_batch_time seconds have elapsed
% or the per-worker termination criteria (max_num_deg_errors, max_num_simulations)
% are reached.
    num_decoded = 0;
    num_error_sum = 0;
    num_degenereated_error_sum = 0;
    num_error_per_iteration_sum = zeros(decoder_property.max_iter, 1);
    num_error_deg_per_iteration_sum = zeros(decoder_property.max_iter, 1);
    num_guesses_total_sum = 0;
    num_guesses_per_iteration_sum = zeros(decoder_property.max_iter, 1);

    tStart = tic;

    while (toc(tStart) < simulation_batch_time) && (num_degenereated_error_sum < max_num_deg_errors) && (num_decoded < max_num_simulations)
        [~, num_error, num_degenereated_error, num_error_per_iteration, num_error_deg_per_iteration, ...
            num_guesses_total, num_guesses_iteration_list]=run_decoding(p_physical, quantum_tanner_code, decoder_property, osd_decoder_x, osd_decoder_z, calc_performance_per_iter);
        
        num_decoded = num_decoded + 1;
        num_error_sum = num_error_sum + num_error;
        num_degenereated_error_sum = num_degenereated_error_sum + num_degenereated_error;
        num_error_per_iteration_sum = num_error_per_iteration_sum + num_error_per_iteration;
        num_error_deg_per_iteration_sum = num_error_deg_per_iteration_sum + num_error_deg_per_iteration;
        num_guesses_total_sum = num_guesses_total_sum + num_guesses_total;
        num_guesses_per_iteration_sum = num_guesses_per_iteration_sum + num_guesses_iteration_list;
    end
end

function [num_decoded, num_error, num_degenereated_error, num_error_per_iteration, num_error_deg_per_iteration, ...
    num_guesses_total, num_guesses_iteration_list]=run_decoding(p_physical, quantum_tanner_code, decoder_property, osd_decoder_x, osd_decoder_z, calc_performance_per_iter)
% Runs a single decoding trial: samples a depolarizing error, computes the
% syndrome, decodes with the selected decoder (correlated or independent X/Z),
% optionally applies OSD post-processing, and returns error and guess counts.

    num_decoded = 1;

    % Depolarization channel
    % [0, X, Z, Y] = [1, 2, 3, 4]
    correct_error_gf4 = randsample([1, 2, 3, 4], quantum_tanner_code.num_qubits, true, [1-p_physical, p_physical/3, p_physical/3, p_physical/3])';
    correct_error_x = (correct_error_gf4 == 2) | (correct_error_gf4 == 4);
    correct_error_z = (correct_error_gf4 == 3) | (correct_error_gf4 == 4);
    correct_error = [correct_error_x; correct_error_z];
    
    syndrome = mod(quantum_tanner_code.H_global * swap_z_x(correct_error), 2);
    syndrome_z = syndrome(1:size(quantum_tanner_code.H_global_z, 1));
    syndrome_x = syndrome((size(quantum_tanner_code.H_global_z, 1) + 1):end);
    
    if decoder_property.use_x_z_correlation
        [error_hat_x, error_hat_z, num_error_per_iteration, num_error_deg_per_iteration, ...
            l_app_x, l_app_z, num_guesses_total, num_guesses_iteration_list] = decode_quantum_tanner_corr( ...
            p_physical, quantum_tanner_code, decoder_property, syndrome_x, syndrome_z, correct_error, osd_decoder_x, osd_decoder_z, calc_performance_per_iter);
    else
        [error_hat_x, error_hat_z, num_error_per_iteration, num_error_deg_per_iteration, ...
            l_app_x, l_app_z, num_guesses_total, num_guesses_iteration_list] = decode_quantum_tanner( ...
            p_physical, quantum_tanner_code, decoder_property, syndrome_x, syndrome_z, correct_error, osd_decoder_x, osd_decoder_z, calc_performance_per_iter);
    end
    if decoder_property.use_osd_decoder
        if ~isequal(mod(quantum_tanner_code.H_global_z * error_hat_x, 2), syndrome_x)
            error_hat_x = osd_decoder_x.decode(syndrome_x, l_app_x);
        end
        if ~isequal(mod(quantum_tanner_code.H_global_x * error_hat_z, 2), syndrome_z)
            error_hat_z = osd_decoder_z.decode(syndrome_z, l_app_z);
        end
    end
   
    error_hat = [error_hat_x; error_hat_z];
    
    error_total = mod(error_hat + correct_error, 2);
    
    num_degenereated_error = double(~all(mod(quantum_tanner_code.H_dual_symplectic * swap_z_x(error_total), 2) == 0, "all"));
    num_error = double(~isequal(error_hat, correct_error));
end