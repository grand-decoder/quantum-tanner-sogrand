function out = col(in)
% Reshapes any input array into a column vector.
%
% Arguments:
%   in  - Input array of any shape or size
%   out - (output) Column vector containing all elements of in

    out = in(:);
end