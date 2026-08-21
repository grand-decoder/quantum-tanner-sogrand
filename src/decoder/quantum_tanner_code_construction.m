function quantum_tanner_code=quantum_tanner_code_construction(lrcc_path, code_path)
% QUANTUM_TANNER_CODE_CONSTRUCTION  Constructs a quantum Tanner code from input files
% from https://github.com/RebKatRad/qTanner, designed in [3].
%
% Inputs:
%   lrcc_path - Path to the LRCC file defining the left-right Cayley complex
%               structure, including qubit mappings and component code indices.
%   code_path - Path to the file containing the local component codes
%               (pa, ga, pb, gb) used to build the parity-check matrices.
%
% Output:
%   quantum_tanner_code - Struct containing the quantum Tanner code data:
%     .delta                    - Size of the delta times delta
%                                 component product code (see [3])
%     .num_qubits               - Total number of physical qubits
%     .num_comp_codes           - Number of component codes
%     .H_comp_x/H_comp_z        - Component code of the Tanner code of the X/Z stabilizers
%     .H_global_x/H_global_z    - Global parity check matrix of the X/Z stabilizers
%     .H_global                 - Full symplectic parity-check matrix [H_x | H_z]
%     .H_dual_symplectic        - Dual symplectic matrix for degenerate error detection
%     .mapping_z/mapping_x      - Qubit mapping for the message passing
%                                 (see further explanation below in import_mapping)

% [3] R. K. Radebold, S. D. Bartlett, and A. C. Doherty, "Explicit instances of quantum Tanner codes," Nov. 14, 2025, arXiv: arXiv:2508.05095.


    [num_comp_codes, num_qubits, delta, mapping_z, mapping_x] = import_mapping(lrcc_path);
    [pa, ga, pb, gb] = import_codes(code_path, delta);
    [H_global_z, H_global_x, ~, ~] = calc_global_H_G(lrcc_path, code_path, delta);
    
    H_global = blkdiag(H_global_x, H_global_z);

    tmp = g2h(H_global);
    H_dual_symplectic = [tmp(:, (num_qubits + 1):end), tmp(:, 1:num_qubits)];

    H_dual_symplectic_x = g2h(H_global_x);
    H_dual_symplectic_z = g2h(H_global_z);
    
    H_comp_z = kron(ga, gb);
    H_comp_x = kron(pa, pb);

    quantum_tanner_code.delta = delta;
    quantum_tanner_code.mapping_z = mapping_z;
    quantum_tanner_code.mapping_x = mapping_x;
    quantum_tanner_code.num_comp_codes = num_comp_codes;
    quantum_tanner_code.num_qubits = num_qubits;
    
    quantum_tanner_code.H_global = H_global;
    quantum_tanner_code.H_global_x = H_global_x;
    quantum_tanner_code.H_global_z = H_global_z;
    quantum_tanner_code.H_dual_symplectic = H_dual_symplectic;
    quantum_tanner_code.H_dual_symplectic_x = H_dual_symplectic_x;
    quantum_tanner_code.H_dual_symplectic_z = H_dual_symplectic_z;
    quantum_tanner_code.H_comp_z = H_comp_z;
    quantum_tanner_code.H_comp_x = H_comp_x;
end

function [num_comp_codes, num_qubits, delta, mapping_z, mapping_x] = import_mapping(filename)
% IMPORT_MAPPING  Parses an LRCC file and builds qubit-to-component-code mappings.
%
% Input:
%   filename - Path to the LRCC file
%
% Outputs:
%   num_comp_codes - Number of component codes
%   num_qubits     - Total number of physical qubits
%   delta          - Size component product code
%   mapping_z/x    - [num_comp_codes, delta^2, 2] arrays mapping local qubit
%                    indices to [global_qubit_idx, qubit_label] for Z/X stabilizers
%
% Explanation mapping_z/x:
% mapping_z[component_code_idx, qubit_local_idx, :]
% = [qubit_global_idx, qubit_label]
%
% For each qubit with local index "qubit_local_idx" in component code
% "component_code_idx", this mapping stores which physical qubit
% "qubit_global_idx" it corresponds to. Since in both the X- and Z-
% stabilizer Tanner code, each qubit is connected to exactly two
% component codes, "qubit_label" indicates whether "component_code_idx"
% is the first or second component code that "qubit_global_idx" is
% connected to. This enables message passing between both component
% codes during decoding.


    raw_data = import_LRCC_file_raw(filename);
    
    num_comp_codes = max(raw_data(:, 2)) / 2;
    num_qubits = max(raw_data(:, 8));
    delta = max(raw_data(:, 7));
    
    % [x_z_label, component_code_idx, qubit_local_idx, :]
    % = [qubit_global_idx, qubit_label]
    % each qubit_global_idx is part of two component codes. The first
    % component code sees qubit_label = 1, the second one sees qubit_label = 2
    mapping = zeros(2, num_comp_codes, delta^2, 2);
    
    qubit_seen = zeros(2, num_qubits);
    for row_idx = 1:size(raw_data, 1)
        current_row = raw_data(row_idx, :);
        x_z_label = current_row(1);
        if x_z_label == 1
            comp_code_idx = current_row(2) - num_comp_codes;
        else
            comp_code_idx = current_row(2);
        end
    
        a = current_row(6);
        b = current_row(7);
    
        local_qubit_idx = (a-1) * delta + b;
        global_qubit_idx = current_row(8);
        if qubit_seen(x_z_label+1, global_qubit_idx) == 1
            qubit_label = 2;
        else
            qubit_label = 1;
            qubit_seen(x_z_label+1, global_qubit_idx) = 1;
        end
    
        mapping(x_z_label+1, comp_code_idx, local_qubit_idx, :) = [global_qubit_idx, qubit_label];
    end
    
    mapping_z = squeeze(mapping(1, :, :, :));
    mapping_x = squeeze(mapping(2, :, :, :));

end

function Q_matrix = import_LRCC_file_raw(filename)
% IMPORT_LRCC_FILE_RAW  Imports an LRCC file from https://github.com/RebKatRad/qTanner in table form.
%
% Input:
%   filename - Path to the LRCC text file
%
% Output:
%   Q_matrix - Matrix where each row encodes the connection between a qubit
%              and a component code:
%              [x_z_label, comp_code_idx1-4, a, b, qubit_idx]
%              - x_z_label:       stabilizer type (0 = Z, 1 = X)
%              - comp_code_idx1:  index of the component code this qubit is connected to
%              - comp_code_idx2-4 (unused in construction): indices of the
%                                 remaining component codes of the square
%              - a, b:            row and column index of the qubit's position
%                                 within the component product code
%              - qubit_idx:       global qubit index
    
    % Read the file as a single string to locate the marker
    try
        fileContent = fileread(filename);
    catch ME
        error('Could not read file: %s', ME.message);
    end

    % Find the position of the header "Q_red"
    headerTag = 'Q_red';
    tagIndex = strfind(fileContent, headerTag);

    if isempty(tagIndex)
        error('The marker "Q_red" was not found in the file.');
    end

    % Count how many newline characters exist before the "Q_red" tag
    % to determine the number of lines to skip.
    linesToSkip = sum(fileContent(1:tagIndex) == newline) + 1;

    % Use readmatrix to import the numeric data.
    % It automatically handles commas as delimiters and ignores trailing commas.
    opts = detectImportOptions(filename, 'FileType', 'text');
    opts.DataLines = [linesToSkip + 1, Inf];
    opts.Delimiter = ',';
    
    % [x_z_label, comp_code_idx1, comp_code_idx2, comp_code_idx3,
    % comp_code_idx4, a, b, qubit_idx]
    Q_matrix = readmatrix(filename, opts);
end

function [pa, ga, pb, gb] = import_codes(path, delta)
% IMPORT_CODES  Loads the four component code matrices from a code file.
%
% Inputs:
%   path  - Path to the component codes file
%   delta - Size component product code
%
% Outputs:
%   pa, ga, pb, gb - Parity-check and generator matrices for the A and B codes

    txt = readlines(path);

    row_indices = [1, 4, 7, 10];
    matrices = cell(1, 4);

    for i = 1:length(row_indices)
        line_str = txt(row_indices(i));
        row_data = str2num(line_str);
        matrices{i} = reshape(row_data, [], delta);
    end
    
    pa = matrices{1};
    ga = matrices{2};
    pb = matrices{3};
    gb = matrices{4};
end

function [H_z, H_x, G_z, G_x]=calc_global_H_G(file_name_LRCC, file_name_code, delta)
% CALC_GLOBAL_H_G  Computes global parity-check and generator matrices for X and Z stabilizers.
%
% Inputs:
%   file_name_LRCC - Path to the LRCC file
%   file_name_code - Path to the component codes file
%   delta          - Size component product code
%
% Outputs:
%   H_z/H_x - Global parity-check matrices for Z/X stabilizers
%   G_z/G_x - Corresponding generator matrices

    [pa, ga, pb, gb] = import_codes(file_name_code, delta);
    raw_data = import_LRCC_file_raw(file_name_LRCC);
    
    c0 = kron(ga, gb);
    c1 = kron(pa, pb);
    
    H_z = construct_global_H_matrix(0, c0, raw_data, delta);
    H_x = construct_global_H_matrix(1, c1, raw_data, delta);
    
    G_z = g2h(H_z);
    G_x = g2h(H_x);
end

function H=construct_global_H_matrix(x_z_label, code, raw_data, delta)
% CONSTRUCT_GLOBAL_H_MATRIX  Assembles the global parity-check matrix for one stabilizer type.
% Inputs:
%   x_z_label - 0 for Z stabilizers, 1 for X stabilizers
%   code      - Component code matrix
%   raw_data  - Raw LRCC data from import_LRCC_file_raw
%   delta     - Local code dimension
%
% Output:
%   H - Global parity-check matrix of size (num_comp_codes * rows_per_code x num_qubits)

    num_comp_codes = floor(max(raw_data(:, 2)) / 2); 
    num_qubits = max(raw_data(:, end));

    dic_array_rows = cell(num_comp_codes + 1, 1);
    
    for row = 1:size(raw_data, 1)
        array_row = raw_data(row, :);
    
        if array_row(1) == x_z_label
            idx = array_row(2);
            if x_z_label == 1
                idx = idx - num_comp_codes;
            end

            % Append the tuple (a, b, qubit_idx) to the cell
            new_entry = [array_row(end-2), array_row(end-1), array_row(end)];
            dic_array_rows{idx} = [dic_array_rows{idx}; new_entry];
        end
    end
    
    rows_per_code = size(code, 1);
    H = zeros(num_comp_codes * rows_per_code, num_qubits);
    
    for comp_code_idx = 1:num_comp_codes
        current_entries = dic_array_rows{comp_code_idx};
        
        if isempty(current_entries)
            continue;
        end
        
        for i = 1:size(current_entries, 1)
            a = current_entries(i, 1);
            b = current_entries(i, 2);
            qubit_idx = current_entries(i, 3);
            
            row_range = ((comp_code_idx - 1) * rows_per_code + 1) : (comp_code_idx * rows_per_code);
            
            code_col = (a - 1) * delta + b;
            H(row_range, qubit_idx) = code(:, code_col);
        end
    end
end