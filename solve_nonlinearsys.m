function [ x, T ] = solve_nonlinearsys( F, J, x0, tol, max_iter )
%SOLVE_NONLINEARSYS
%
%   [ x, T ] = solve_nonlinearsys( F, J, x0, tol, max_iter )
%
%   Approximates a solution to a nonlinear system of equations using 
%   Newton's method.
%
%PARAMETERS:
%   F           The system to approximate the solution to.
%
%   J           The Jacobian matrix for the system.
%
%   x0          The initial guess.
%
%   tol         Tolerance.
%
%   max_iter    Maximum number of iterations to be performed.
%
%RETURNS:
%   x           The approximate solution.
%
%   T           Table of approximations (approximate solution at each 
%               iteration).
%
%EXAMPLE:
%   Suppose we want the solution to the following non-linear system
%   of equations:
%       f_1(x,y,z) = x^3 - 2y - 2
%       f_2(x,y,z) = x^3 - 5z^2 + 7
%       f_3(x,y,z) = yz^2 - 1
%
%   The following code gives the approximate solution. 
%-------------------------------------------------------------------------
%   f = @(X) [ X(1)^3-2*X(2)-2  ; ...
%              X(1)^3-5*X(3)^2+7; ...
%              X(2)*X(3)^2-1 ];
%   j = @(X) [ 3*X(1)^2,      -2,            0 ; ...
%              3*X(1)^2,       0,     -10*X(3) ; ...
%                     0,  X(3)^2,  2*X(2)*X(3) ];
%   sol = solve_nonlinearsys(f,j, [1,1,1], 5.0e-5, 10)
%-------------------------------------------------------------------------
% sol =
% 
%    1.442249570335223
%    0.500000000014800
%    1.414213562375909
%
%   Thus, x = 1.44225, y = 0.50000, z = 1.41421.
%
%AUTHOR:    Kenneth Poppie
%DATE:      Oct 28, 2016

% Make sure x0 is in column form.
if size(x0,1) < size(x0,2) 
    x_init = x0';
else
    x_init = x0;
end

% Initialize the return values.
x = x_init;
T = [x_init'];
iter = 1;

while iter < max_iter 
    v = J(x)\-F(x);
    x = x + v;
    T = [T; x'];
    
    % Check tolerance
    if max(abs(v)) < tol 
        
        % tolerance met
        break;
    end
    
    iter = iter + 1;
end
end

