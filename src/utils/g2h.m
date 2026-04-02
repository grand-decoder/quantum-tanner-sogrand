function H=g2h(G)
% G2H  Converts a generator matrix to its parity-check matrix over GF(2).
%
% Performs Gaussian elimination over GF(2) to bring G into reduced row
% echelon form, tracking all column swaps. The parity-check matrix is
% then computed via gen2par, and column swaps are reversed to restore
% the original qubit ordering.
%
% Input:
%   G - Generator matrix over GF(2)
%
% Output:
%   H - Corresponding parity-check matrix over GF(2)

    H = G;
    [rows, cols] = size(H);
    colSwap = [];
    
    for i = 1:min(rows, cols)
        % --- Step 1: Find a Pivot in the entire remaining block ---
        if H(i, i) == 0
            % Extract the sub-matrix from (i,i) to the end
            subMatrix = H(i:rows, i:cols);
            
            % Find the FIRST non-zero entry in this sub-matrix
            % 'find' returns linear indices, we need (r, c) relative to subMatrix
            [r_rel, c_rel] = find(subMatrix, 1);
            
            if isempty(r_rel)
                % If the entire sub-matrix is zeros, we are done (rank deficient)
                continue; 
            end
            
            % Convert relative indices to absolute indices
            pivot_row = i + r_rel - 1;
            pivot_col = i + c_rel - 1;
            
            % --- Step 2: Swap to bring pivot to (i,i) ---
            
            % Swap Rows if needed
            if pivot_row ~= i
                H([i, pivot_row], :) = H([pivot_row, i], :);
            end
            
            % Swap Columns if needed
            if pivot_col ~= i
                H(:, [i, pivot_col]) = H(:, [pivot_col, i]);
                colSwap = [colSwap; i pivot_col]; % Record the swap
            end
        end
        
        % --- Step 3: Elimination (Gaussian) ---
        % Now we are guaranteed H(i,i) == 1 (unless matrix is singular)
        for r = 1:rows
            if r ~= i && H(r, i) == 1
                H(r, :) = xor(H(r, :), H(i, :));
            end
        end
    end
    
    H(all(H == 0, 2), :) = [];

    H = gen2par(H);
    % swap columns
    nSwap = size(colSwap, 1);
    for s=nSwap:-1:1
        ord = colSwap(s,:);
        nOrd = ord(end:-1:1);
        H(:,ord) = H(:,nOrd);
    end
end
