function v = swap_z_x(v)
% SWAP_Z_X  Swaps the X and Z halves of a symplectic vector.
%
% Symplectic vectors are stored as [Z | X]. This function swaps the two
% halves to convert between [Z | X] and [X | Z] representations.
%
% Input:
%   v - 1D symplectic vector of even length
%
% Output:
%   v - Vector with first and second halves swapped

    % Check if input is a 1D vector (row or column)
    if ~isvector(v)
        error('Input must be a 1D vector.');
    end

    n = numel(v);

    % Check if the number of entries is even
    if mod(n, 2) ~= 0
        error('Input must have an even number of entries.');
    end

    % Swap the first and second half
    % Indexing automatically preserves the shape (row vs column)
    mid = n / 2;
    v = v([mid+1:end, 1:mid]);
end