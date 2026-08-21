% Decode quantum tanner code. decodes x and z errors combined.

function [error_hat_x, error_hat_z, normal_error_list, deg_error_list, l_app_x, l_app_z, n_guess_total, n_guess_iteration_list] = decode_quantum_tanner_relay_x_z( ...
    p_physical, quantum_tanner_code, decoder_property, syndrome_global_x, syndrome_global_z, correct_error_x, correct_error_z, calc_performance_per_iter)
% Iterative belief propagation decoder for quantum Tanner codes with relay principle [6] message updates.
%
% Inputs:
%   p_physical               - Physical error probability p (depolarizing channel:
%                              P(I)=1-p, P(X)=P(Z)=P(Y)=p/3)
%   quantum_tanner_code      - Struct with code parameters (parity-check matrices,
%                              qubit mappings, component codes, etc.)
%   decoder_property         - Struct with decoder settings (max_iter, list_size,
%                              sogrand_thr, alpha, IC, use_osd_decoder, and relay parameters)
%   syndrome_global_x        - Measured syndrome of X errors (Z stabilizers) (column vector)
%   syndrome_global_z        - Measured syndrome of Z errors (X stabilizers) (column vector)
%   correct_error_x          - True error vector [x-part], used for performance evaluation only
%   correct_error_z          - True error vector [z-part], used for performance evaluation only
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

    [error_hat_x, normal_error_list_x, deg_error_list_x, l_app_x, n_guess_total_x, n_guess_iteration_list_x] = decode_quantum_tanner_relay( ...
        p_physical, quantum_tanner_code.num_qubits,  quantum_tanner_code.num_comp_codes, quantum_tanner_code.mapping_z, quantum_tanner_code.H_comp_z, ...
        quantum_tanner_code.H_global_z, quantum_tanner_code.H_dual_symplectic_x, decoder_property, syndrome_global_x, correct_error_x, calc_performance_per_iter);

    [error_hat_z, normal_error_list_z, deg_error_list_z, l_app_z, n_guess_total_z, n_guess_iteration_list_z] = decode_quantum_tanner_relay( ...
        p_physical, quantum_tanner_code.num_qubits, quantum_tanner_code.num_comp_codes, quantum_tanner_code.mapping_x, quantum_tanner_code.H_comp_x, ...
        quantum_tanner_code.H_global_x, quantum_tanner_code.H_dual_symplectic_z, decoder_property, syndrome_global_z, correct_error_z, calc_performance_per_iter);

    n_guess_iteration_list = n_guess_iteration_list_x + n_guess_iteration_list_z;
    n_guess_total = n_guess_total_x + n_guess_total_z;

    % decoding is in error if either x or z part is in error. "or" corresponds here to "max".
    normal_error_list = max(normal_error_list_x, normal_error_list_z);
    deg_error_list = max(deg_error_list_x, deg_error_list_z);
end

function [error_hat, normal_error_list_total, deg_error_list_total, l_app, n_guess_total, n_guess_iteration_list_total] = decode_quantum_tanner_relay( ...
    p_physical, num_qubits, num_comp_codes, mapping, H_comp, H_global, H_dual_symplectic, decoder_property, syndrome_global, correct_error, calc_performance_per_iter)

    p_tilde = 2 * p_physical / 3;
    l_ch = zeros(num_qubits, 1);
    l_ch(:) = log((1-p_tilde) ./ p_tilde);

    total_num_iter = decoder_property.rel_bp_pre_iter + (decoder_property.rel_bp_num_sets - 1) * decoder_property.rel_bp_set_max_iter;
    normal_error_list_total = [];
    deg_error_list_total = [];
    n_guess_iteration_list_total = [];

    l_app_previous_leg = l_ch;
    for leg_idx = 1:decoder_property.rel_bp_num_sets
        first_set = leg_idx == 1 && ~decoder_property.rel_bp_skip_first_leg;

        [dec_successful, error_hat, normal_error_iteration_list, deg_error_iteration_list, l_app, n_guess_total, n_guess_iteration_list] = relay_leg( ...
            first_set, l_ch, num_qubits, num_comp_codes, mapping, H_comp, H_global, H_dual_symplectic, decoder_property, syndrome_global, ...
            correct_error, calc_performance_per_iter, l_app_previous_leg);
        
        normal_error_list_total = [normal_error_list_total; normal_error_iteration_list];
        deg_error_list_total = [deg_error_list_total; deg_error_iteration_list];
        n_guess_iteration_list_total = [n_guess_iteration_list_total; n_guess_iteration_list];

        l_app_previous_leg = l_app;

        if dec_successful
            break;
        end
    end
    
    if numel(normal_error_list_total) < total_num_iter
        % In case we terminate early because of decoding success, the
        % remaining entries are filled with the last decoding result.
        entries_to_add = total_num_iter - numel(normal_error_list_total);
        normal_error_list_total = [normal_error_list_total; normal_error_list_total(end) * ones(entries_to_add, 1)];
        deg_error_list_total = [deg_error_list_total; deg_error_list_total(end) * ones(entries_to_add, 1)];
        n_guess_iteration_list_total = [n_guess_iteration_list_total; zeros(entries_to_add, 1)];
    end
end

function [dec_successful, error_hat, normal_error_iteration_list, deg_error_iteration_list, l_app, n_guess_total, n_guess_iteration_list] = relay_leg( ...
    first_set, l_bias_init, num_qubits, num_comp_codes, mapping, H_comp, H_global, H_dual_symplectic, decoder_property, syndrome_global, correct_error, ...
    calc_performance_per_iter, l_app_previous_leg)
    
    dec_successful = false;

    l_e_current = zeros(num_qubits, 2);
    
    % Select/sample memory strengths
    if first_set
        memory_strength_list = decoder_property.rel_bp_gamma0 * ones(num_qubits, 1);
        max_iter = decoder_property.rel_bp_pre_iter;
    else
        memory_strength_list = decoder_property.rel_bp_gamma_dist_interval_start + ...
            (decoder_property.rel_bp_gamma_dist_interval_end - decoder_property.rel_bp_gamma_dist_interval_start) * rand(num_qubits, 1);
        max_iter = decoder_property.rel_bp_set_max_iter;
    end
    
    normal_error_iteration_list = zeros(max_iter, 1);
    deg_error_iteration_list = zeros(max_iter, 1);
    
    n_guess_total = 0;
    n_guess_iteration_list = zeros(max_iter, 1);
    
    l_app_previous_iter = l_app_previous_leg;
    for iteration=1:max_iter
        % DMem-BP
        l_bias_iteration = (1 - memory_strength_list) .* l_bias_init + memory_strength_list .* l_app_previous_iter;
        [l_e_current, n_guess_iter, l_app, error_hat] = belief_propagation_update(l_e_current, num_qubits, ...
            mapping, H_comp, syndrome_global, num_comp_codes, decoder_property, l_bias_iteration);
        l_app_previous_iter = l_app;

        n_guess_total = n_guess_total + n_guess_iter;
        n_guess_iteration_list(iteration) = n_guess_iteration_list(iteration) + n_guess_iter;
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is only for error counting and does not affect decoding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if calc_performance_per_iter
            error_hat_only_counting = error_hat;
            error_total_only_counting = mod(error_hat_only_counting + correct_error, 2);
            is_normal_error = ~isequal(error_hat_only_counting, correct_error);
            is_deg_error = ~all(mod(H_dual_symplectic * error_total_only_counting, 2) == 0, "all");
        else
            is_normal_error = false;
            is_deg_error = false;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is only for error counting and does not affect decoding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Termination Condition
        current_syndrome = mod(H_global * error_hat, 2);
        if isequal(current_syndrome, syndrome_global)
            normal_error_iteration_list(iteration:end) = is_normal_error;
            deg_error_iteration_list(iteration:end) = is_deg_error;
            dec_successful = true;
            break;
        else
            normal_error_iteration_list(iteration) = is_normal_error;
            deg_error_iteration_list(iteration) = is_deg_error;
        end
    end
end

function [le_new, n_guess_iter, l_app, error_hat] = belief_propagation_update(l_e_current, num_qubits, mapping, H_comp, syndrome_global, ...
    num_component_codes, decoder_property, l_bias)
% Applies one belief propagation update including constraint node and variable node update for all component codes and qubits.
% For each component code, incoming extrinsic LLRs and channel LLRs are combined,
% then passed to the SOGRAND component decoder. The resulting output LLRs are
% written back at the corresponding qubit positions.
%
% Inputs:
%   l_e_current         - Current extrinsic LLR matrix (num_qubits x 2) containing the extrinsic LLRs that each variable node 
%                         sends to its two connected component codes
%   num_qubits          - Total number of qubits in the code
%   mapping             - (num_component_codes x n_comp x 2) array mapping each
%                         component code view to global qubit indices (see quantum_tanner_code_construction.m for details)
%   H_comp              - Parity-check matrix of the component code.
%   syndrome_global     - Concatenated syndromes for all component codes
%   num_component_codes - Number of component codes
%   decoder_property    - Struct with decoder settings (alpha, list_size,
%                         sogrand_thr, IC)
%   l_bias              - Bias LLR vector (corresponds the channel LLR 
%                         vector in normal BP but is here iteration 
%                         dependent [6]) (num_qubits x 1)
%
% Outputs:
%   le_new       - Updated extrinsic LLR matrix containing the new extrinsic LLRs that each variable node sends to its two connected component codes 
%                  in the next iteration. Same format as l_e_current.
%   n_guess_iter - Total number of SOGRAND guesses across all component codes
%   l_app        - A-posteriori LLRs
%   error_hat    - Estimated error pattern
    m_comp = size(H_comp, 1);
    n_guess_iter = 0;
    
    le_new = zeros(num_qubits, 2);
    for comp_code_idx=1:num_component_codes
        view_qubit_global_indices = col(mapping(comp_code_idx, :, 1));
        view_labels = col(mapping(comp_code_idx, :, 2));

        % Obtain the incoming extrinsic LLRs at the current component code by indexing the extrinsic LLRs that the respective other component code sents
        % via view_inverted_labels.
        view_inverted_labels = 3 - view_labels;
        input_l_e_linear_indices = sub2ind(size(l_e_current), view_qubit_global_indices, view_inverted_labels);

        output_l_e_linear_indices = sub2ind(size(l_e_current), view_qubit_global_indices, view_labels);
        
        le_input = decoder_property.alpha * l_e_current(input_l_e_linear_indices);
        l_in = le_input + l_bias(view_qubit_global_indices);

        % Extract syndrome corresponding to current component code.
        view_syndrome = syndrome_global((m_comp * (comp_code_idx - 1) + 1):(m_comp * comp_code_idx));
        
        [le_output, n_guess] = component_decoding(l_in, decoder_property.list_size, decoder_property.sogrand_thr, H_comp, view_syndrome, decoder_property.IC);
        n_guess_iter = n_guess_iter + n_guess;

        le_new(output_l_e_linear_indices) = le_output;
    end

    l_app = sum(le_new, 2) + l_bias;
    error_hat = l_app < 0;
end

