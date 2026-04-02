% Decode quantum tanner code. decodes x and z errors combined.

function [error_hat_x, error_hat_z, normal_error_iteration_list, deg_error_iteration_list, l_app_x, l_app_z, n_guess_total, n_guess_iteration_list] = decode_quantum_tanner( ...
    p_physical, quantum_tanner_code, decoder_property, syndrome_global_x, syndrome_global_z, correct_error, osd_decoder_x, osd_decoder_z, calc_performance_per_iter)
% Iterative belief propagation decoder for quantum Tanner codes.
% X and Z errors are decoded in parallel using their respective component codes.
% Decoding stops early if both syndromes are satisfied. Optionally applies an
% OSD post-processor and tracks per-iteration error statistics.
%
% Inputs:
%   p_physical               - Physical error probability
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
%   error_hat_x                - Estimated X error pattern (binary column vector)
%   error_hat_z                - Estimated Z error pattern (binary column vector)
%   normal_error_iteration_list - Per-iteration flag: 1 if a logical error occurred
%   deg_error_iteration_list    - Per-iteration flag: 1 if a non-degenerate logical
%                                 error occurred
%   l_app_x                    - A-posteriori LLRs for X errors
%   l_app_z                    - A-posteriori LLRs for Z errors
%   n_guess_total              - Total number of SOGRAND guesses across all iterations
%   n_guess_iteration_list     - Per-iteration SOGRAND guess counts

    % l_e_current[qubit_idx, :] = [l_e1, l_e2] 
    % stores the extrinsic LLRs that variable node qubit_idx sends to its two connected component codes.
    l_e_current_x = zeros(quantum_tanner_code.num_qubits, 2);
    l_e_current_z = zeros(quantum_tanner_code.num_qubits, 2);
    
    % all-zero is the most likely error pattern.
    p_tilde = 2 * p_physical / 3;
    l_ch = zeros(quantum_tanner_code.num_qubits, 1);
    l_ch(:) = log((1-p_tilde) ./ p_tilde);
    
    % Set to 1 if syndrome conditions for X or Z errors are satisfied, respectively. Decoding stops when both are 1.
    stop_x = 0;
    stop_z = 0;
    
    normal_error_iteration_list = zeros(decoder_property.max_iter, 1);
    deg_error_iteration_list = zeros(decoder_property.max_iter, 1);
    
    n_guess_total=0;
    n_guess_iteration_list = zeros(decoder_property.max_iter, 1);
    
    for iteration=1:decoder_property.max_iter
        if ~stop_x
            [le_new_x, n_guess_iter_x] = belief_propagation_update(l_e_current_x, quantum_tanner_code.num_qubits, ...
                quantum_tanner_code.mapping_z, quantum_tanner_code.H_comp_z, syndrome_global_x, ...
                quantum_tanner_code.num_comp_codes, decoder_property, l_ch);
            l_e_current_x = le_new_x;
            l_app_x = sum(l_e_current_x, 2) + l_ch;
            error_hat_x = l_app_x < 0;
            n_guess_total = n_guess_total + n_guess_iter_x;
            n_guess_iteration_list(iteration) = n_guess_iteration_list(iteration) + n_guess_iter_x;
        end
        if ~stop_z
            [le_new_z, n_guess_iter_z] = belief_propagation_update(l_e_current_z, quantum_tanner_code.num_qubits, ...
                quantum_tanner_code.mapping_x, quantum_tanner_code.H_comp_x, syndrome_global_z, ...
                quantum_tanner_code.num_comp_codes, decoder_property, l_ch);
            l_e_current_z = le_new_z;
            l_app_z = sum(l_e_current_z, 2) + l_ch;
            error_hat_z = l_app_z < 0;
            n_guess_total = n_guess_total + n_guess_iter_z;
            n_guess_iteration_list(iteration) = n_guess_iteration_list(iteration) + n_guess_iter_z;
        end
    
        % Termination Condition
        current_syndrome_x = mod(quantum_tanner_code.H_global_z * error_hat_x, 2);
        current_syndrome_z = mod(quantum_tanner_code.H_global_x * error_hat_z, 2);
        if isequal(current_syndrome_x, syndrome_global_x) 
            stop_x = 1;
        end
        if isequal(current_syndrome_z, syndrome_global_z)
            stop_z = 1;
        end
        
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

            % count errors.
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
    
        if stop_x && stop_z
            normal_error_iteration_list(iteration:end) = is_normal_error;
            deg_error_iteration_list(iteration:end) = is_deg_error;
            break;
        else
            normal_error_iteration_list(iteration) = is_normal_error;
            deg_error_iteration_list(iteration) = is_deg_error;
        end
    end


end

function [le_new, n_guess_iter] = belief_propagation_update(l_e_current, num_qubits, mapping, H_comp, syndrome_global, ...
    num_component_codes, decoder_property, l_ch)
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
%   l_ch                - Channel LLR vector (num_qubits x 1)
%
% Outputs:
%   le_new       - Updated extrinsic LLR matrix containing the new extrinsic LLRs that each variable node sends to its two connected component codes 
%                  in the next iteration. Same format as l_e_current.
%   n_guess_iter - Total number of SOGRAND guesses across all component codes

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
        l_in = le_input + l_ch(view_qubit_global_indices);

        % Extract syndrome corresponding to current component code.
        view_syndrome = syndrome_global((m_comp * (comp_code_idx - 1) + 1):(m_comp * comp_code_idx));
        
        [le_output, n_guess] = component_decoding(l_in, decoder_property.list_size, decoder_property.sogrand_thr, H_comp, view_syndrome, decoder_property.IC);
        n_guess_iter = n_guess_iter + n_guess;

        le_new(output_l_e_linear_indices) = le_output;
    end
end

