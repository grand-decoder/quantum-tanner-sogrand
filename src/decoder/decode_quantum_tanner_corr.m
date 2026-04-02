% Decode quantum tanner code. decodes x and z errors combined.

function [error_hat_x, error_hat_z, normal_error_list, deg_error_list, l_app_x, l_app_z, n_guess_total, n_guess_iteration_list] = decode_quantum_tanner_corr( ...
    p_physical, quantum_tanner_code, decoder_property, syndrome_global_x, syndrome_global_z, correct_error, osd_decoder_x, osd_decoder_z, calc_performance_per_iter)
% Iterative belief propagation decoder for quantum Tanner codes with correlation aware message updates.
% While [1] describes the decoder in probability domain, this implementation operates in the log-domain for numerical stability and efficiency.
%
% [1] L. Rapp, M. Médard, E. Tang, and K. R. Duffy, 
% "Efficient Soft-Output Guessing for Enhanced Quantum Tanner Code Decoding," 
% arXiv preprint arXiv:2603.18318, 2026.
%
% Inputs:
%   p_physical               - Physical error probability p (depolarizing channel:
%                              P(I)=1-p, P(X)=P(Z)=P(Y)=p/3)
%   quantum_tanner_code      - Struct with code parameters (parity-check matrices,
%                              qubit mappings, component codes, etc.)
%   decoder_property         - Struct with decoder settings (max_iter, list_size,
%                              sogrand_thr, alpha, IC, use_osd_decoder)
%   syndrome_global_x        - Measured syndrome of X errors (Z stabilizers) (column vector)
%   syndrome_global_z        - Measured syndrome of Z errors (X stabilizers) (column vector)
%   correct_error            - True error vector [x-part; z-part], used for
%                              performance evaluation only
%   osd_decoder_x            - OSD decoder object for X errors (used if
%                              decoder_property.use_osd_decoder is set)
%   osd_decoder_z            - OSD decoder object for Z errors
%   calc_performance_per_iter - If true, compute error statistics at each iteration
%
% Outputs:
%   error_hat_x              - Estimated X error pattern (binary column vector)
%   error_hat_z              - Estimated Z error pattern (binary column vector)
%   normal_error_list        - Per-iteration flag: 1 if a logical error occurred
%   deg_error_list           - Per-iteration flag: 1 if a non-degenerate logical
%                              error occurred
%   l_app_x                  - A-posteriori LLRs for X errors
%   l_app_z                  - A-posteriori LLRs for Z errors
%   n_guess_total            - Total number of SOGRAND guesses across all iterations
%   n_guess_iteration_list   - Per-iteration SOGRAND guess counts

    assert(decoder_property.alpha == 1, 'At the moment only alpha=1 is supported.');
    
    % --- Pre-compute Channel Log-Likelihoods ---
    % P(I) = 1-p, P(X)=P(Z)=P(Y)=p/3
    p_ch_vec = [
        (1-p_physical)  % I
        p_physical/3    % X
        p_physical/3    % Z
        p_physical/3    % Y
    ];
    l_ch_vec = log(p_ch_vec); 
    % Reshape for broadcasting
    l_ch_vec = reshape(l_ch_vec, 1, 1, 4);
    
    % --- Initialize Messages ---
    p_ch_x_or_z = 2 * p_physical / 3;
    init_llr = log((1-p_ch_x_or_z) ./ p_ch_x_or_z);
    
    l_a_x_v2c = zeros(quantum_tanner_code.num_qubits, 2);
    l_a_x_v2c(:) = init_llr;
    
    l_a_z_v2c = zeros(quantum_tanner_code.num_qubits, 2);
    l_a_z_v2c(:) = init_llr;
    
    normal_error_list = zeros(decoder_property.max_iter, 1);
    deg_error_list = zeros(decoder_property.max_iter, 1);
    
    n_guess_total=0;
    n_guess_iteration_list = zeros(decoder_property.max_iter, 1);
    
    for iteration=1:decoder_property.max_iter
        % Constraint node update
        [l_e_x_c2v, n_guess_iter_x] = constraint_node_update(l_a_x_v2c, quantum_tanner_code.num_qubits, ...
            quantum_tanner_code.mapping_z, quantum_tanner_code.H_comp_z, syndrome_global_x, quantum_tanner_code.num_comp_codes, decoder_property);

        [l_e_z_c2v, n_guess_iter_z] = constraint_node_update(l_a_z_v2c, quantum_tanner_code.num_qubits, ...
            quantum_tanner_code.mapping_x, quantum_tanner_code.H_comp_x, syndrome_global_z, quantum_tanner_code.num_comp_codes, decoder_property);
        
        n_guess_total = n_guess_total + n_guess_iter_x + n_guess_iter_z;
        n_guess_iteration_list(iteration) = n_guess_iteration_list(iteration) + n_guess_iter_x + n_guess_iter_z;
    
        % Convert incoming messages into symbol probabilities (here implemented in log domain) by assuming
        % that a message coming from X knows nothing about Z, hence, prob
        % for Z is 1/2 and vise versa.
        l_e_x_c2v_sym = zeros(quantum_tanner_code.num_qubits, 2, 4); %0, X, Z, Y
        l_e_x_c2v_sym(:, :, 2) = -l_e_x_c2v; % X
        l_e_x_c2v_sym(:, :, 4) = -l_e_x_c2v; % Y
    
        l_e_z_c2v_sym = zeros(quantum_tanner_code.num_qubits, 2, 4);
        l_e_z_c2v_sym(:, :, 3) = -l_e_z_c2v; % Z
        l_e_z_c2v_sym(:, :, 4) = -l_e_z_c2v; % Y
        
        l_total_sum = l_ch_vec + l_e_x_c2v_sym(:, 1, :) + l_e_x_c2v_sym(:, 2, :) + l_e_z_c2v_sym(:, 1, :) + l_e_z_c2v_sym(:, 2, :);
        
        l_a_x_v2c_sym = zeros(quantum_tanner_code.num_qubits, 2, 4);
        l_a_x_v2c_sym(:, 1, :) = l_total_sum - l_e_x_c2v_sym(:, 1, :);
        l_a_x_v2c_sym(:, 2, :) = l_total_sum - l_e_x_c2v_sym(:, 2, :);
    
        l_a_z_v2c_sym = zeros(quantum_tanner_code.num_qubits, 2, 4);
        l_a_z_v2c_sym(:, 1, :) = l_total_sum - l_e_z_c2v_sym(:, 1, :);
        l_a_z_v2c_sym(:, 2, :) = l_total_sum - l_e_z_c2v_sym(:, 2, :);

        num_x = maxstar(l_a_x_v2c_sym(:, :, 1), l_a_x_v2c_sym(:, :, 3));
        den_x = maxstar(l_a_x_v2c_sym(:, :, 2), l_a_x_v2c_sym(:, :, 4));
        l_a_x_v2c = num_x - den_x;
    
        num_z = maxstar(l_a_z_v2c_sym(:, :, 1), l_a_z_v2c_sym(:, :, 2));
        den_z = maxstar(l_a_z_v2c_sym(:, :, 3), l_a_z_v2c_sym(:, :, 4));
        l_a_z_v2c = num_z - den_z;
        
        l_app_vec = squeeze(l_total_sum);
        l_app_x_num = maxstar(l_app_vec(:, 1), l_app_vec(:, 3));
        l_app_x_den = maxstar(l_app_vec(:, 2), l_app_vec(:, 4));
        l_app_x = l_app_x_num - l_app_x_den;
        
        l_app_z_num = maxstar(l_app_vec(:, 1), l_app_vec(:, 2));
        l_app_z_den = maxstar(l_app_vec(:, 3), l_app_vec(:, 4));
        l_app_z = l_app_z_num - l_app_z_den;
    
        error_hat_z = l_app_z < 0;
        error_hat_x = l_app_x < 0;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is only for error counting and does not affect decoding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if calc_performance_per_iter
            error_hat_x_only_counting = error_hat_x;
            error_hat_z_only_counting = error_hat_z;

            if decoder_property.use_osd_decoder
                if ~isequal(mod(quantum_tanner_code.H_global_z * error_hat_x_only_counting, 2), syndrome_global_x)
                    error_hat_x_only_counting = osd_decoder_x.decode(syndrome_global_x, l_app_x);
                end
                if ~isequal(mod(quantum_tanner_code.H_global_x * error_hat_z_only_counting, 2), syndrome_global_z)
                    error_hat_z_only_counting = osd_decoder_z.decode(syndrome_global_z, l_app_z);
                end
            end

            error_hat_only_counting = [error_hat_x_only_counting; error_hat_z_only_counting];

            error_total_only_counting = mod(error_hat_only_counting + correct_error, 2);
            is_normal_error = ~isequal(error_hat_only_counting, correct_error);
            is_deg_error = ~all(mod(quantum_tanner_code.H_dual_symplectic * swap_z_x(error_total_only_counting), 2) == 0, "all");
        else
            is_normal_error = false;
            is_deg_error = false;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is only for error counting and does not affect decoding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Termination Condition
        current_syndrome_x = mod(quantum_tanner_code.H_global_z * error_hat_x, 2);
        current_syndrome_z = mod(quantum_tanner_code.H_global_x * error_hat_z, 2);
    
        if isequal(current_syndrome_x, syndrome_global_x) && isequal(current_syndrome_z, syndrome_global_z)
            normal_error_list(iteration:end) = is_normal_error;
            deg_error_list(iteration:end) = is_deg_error;
            break;
        else
            normal_error_list(iteration) = is_normal_error;
            deg_error_list(iteration) = is_deg_error;
        end
    end
end

function [le_new, n_guess_iter] = constraint_node_update(l_a_v2c, num_qubits, mapping, H_comp, syndrome_global, ...
    num_component_codes, decoder_property)
% Constraint node update step for one iteration of belief propagation.
% In contrast to decode_quantum_tanner, here we only carry out the constraint node update while the variable node
% is integrated into the message update above.
% For each component code, the incoming a-priori V2C LLRs are passed to the
% SOGRAND component decoder, and the resulting extrinsic C2V LLRs are written
% back at the corresponding qubit positions.
%
% Inputs:
%   l_a_v2c             - A-priori variable-to-constraint LLR matrix (num_qubits x 2)
%   num_qubits          - Total number of qubits in the code
%   mapping             - (num_component_codes x n_comp x 2) array mapping each
%                         component code view to global qubit indices (dim 3 = 1)
%                         and label directions (dim 3 = 2)
%   H_comp              - Parity-check matrix of the component code
%   syndrome_global     - Concatenated syndromes for all component codes
%   num_component_codes - Number of component codes
%   decoder_property    - Struct with decoder settings (list_size, sogrand_thr, IC)
%
% Outputs:
%   le_new       - Extrinsic constraint node to variable node LLR matrix (num_qubits x 2)
%   n_guess_iter - Total number of SOGRAND guesses across all component codes

    m_comp = size(H_comp, 1);    

    n_guess_iter = 0;
    
    le_new = zeros(num_qubits, 2);
    for comp_code_idx=1:num_component_codes
        view_qubit_global_indices = col(mapping(comp_code_idx, :, 1));

        view_labels = col(mapping(comp_code_idx, :, 2));
        l_e_linear_indices = sub2ind(size(l_a_v2c), view_qubit_global_indices, view_labels);

        l_in = l_a_v2c(l_e_linear_indices);

        view_syndrome = syndrome_global((m_comp * (comp_code_idx - 1) + 1):(m_comp * comp_code_idx));
        
        [le_output, n_guess] = component_decoding(l_in, decoder_property.list_size, decoder_property.sogrand_thr, H_comp, view_syndrome, decoder_property.IC);
        n_guess_iter = n_guess_iter + n_guess;

        le_new(l_e_linear_indices) = le_output;
    end
end

function y = maxstar(a, b)
% Computes the Jacobi logarithm (max* operator): max*(a,b) = log(exp(a) + exp(b)).
%
% Inputs:
%   a, b - Scalar or array log-domain values (same size)
%
% Output:
%   y    - Element-wise max*(a, b), same size as inputs

    % max*(a,b) = max(a,b) + log(1 + exp(-abs(a-b)))
    % Handles -Inf correctly
    y = max(a, b) + log(1 + exp(-abs(a - b)));
end
