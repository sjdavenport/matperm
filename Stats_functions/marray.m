function [ outarray ] = marray( M, array )
% MARRAY(M, array) 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% M:  a matrix of dimension m by n
% array: a D-dimensional array with final dimension n
%--------------------------------------------------------------------------
% OUTPUT
% outarray: a D-dimensional array 
%--------------------------------------------------------------------------
% EXAMPLES
% M = [1,1;0,1];
% s_array = ones([3,3,2]);
% marray(M, s_array)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_M = size(M);
s_array = size(array);
D = length(s_array);
outarray = zeros([s_array(1:end-1), s_M(1)]);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
index  = repmat( {':'}, 1, D - 1 );
for j = 1:s_M(1)
    for k = 1:s_array(end)
        outarray( index{:}, j ) = outarray( index{:}, j ) +...
                 M( j, k )*array( index{:}, k );
    end
end

end

