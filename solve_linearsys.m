function [ x, X, iter ] = solve_linearsys(A, b, x0, tol, max_iter, method)
%SOLVE_LINEARSYS
%
%   [x, X, iter] = solve_linearsys(A, b, x0, tol, max_iter, method)
%
%   Approximates the solution to the system Ax=b using the Jacobi 
%   iterative method or the Gauss-Seidel iterative method.
%
%PARAMETERS:
%   A           n by n matrix.
%
%   b           n by 1 column vector.
%
%   x0          n by 1 column vector (the initial guess).
%
%   tol         The tolerance (a zero tolerance will cause the routine
%                   to perform the maximum number of iterations).
%
%   max_iter    The maximum number of iterations the routine will 
%                   perform before stopping.
%
%   method      'jacobi' for the Jacobi method
%               'gauss'  for the Gauss-Seidel method
%
%RETURNS:
%   x           Approximate solution as column vector.
%
%   X           A matrix whose columns are the best solution at each
%                   iteration of the routine.
%
%   iter        The number of iterations the routine actually performed.
%
%AUTHOR:    Kenneth Poppie
%DATE:      10/20/2016
%MODIFIED:  11/6/2016 by Kenneth Poppie
%           Removed the need for the MATLAB backslash operator.
%

% Determine the method to be used.
if lower(method(1)) == 'j' % Jacobi method!!!
    iter_method = 'j';
elseif lower(method(1)) == 'g' % Gauss-Seidel method!!!
    iter_method = 'g';
else 
    error('Unknown iterative method: %s\n%s', method, ...
        'Type ''help solve_linearsys'' for valid methods.');  
end

% Split A into N and P
if iter_method == 'j'
    % Jacobi split
    N = diag(A);
    P = diag(N) - A;
else 
    % Gauss-Seidel split
    N = tril(A);
    P = N - A;
end
   
% Initialize working variables
X = [];
new_x = [];
old_x = x0;

% Begin iterations
for iter = 1:max_iter
    % new_x = N\(P*old_x + b);
    if iter_method == 'j'
        % Jacobi
        new_x = (P*old_x + b)./N;
    else
        % Gauss-Seidel
        new_x = forward_sub(N,(P*old_x+b));
    end

    % Save estimate
    X = [X new_x];

    % Check tolerance
    if norm(old_x - new_x,'inf')<=tol
        break
    end

    % Next iteration
    old_x = new_x;
end
x = new_x;

end % function

% Helper function.
% Solves Lx=b using forward substituion.
%
function [x] = forward_sub(L,b)
    n = max(size(b,1),size(b,2));
    x=zeros(n,1);
    for j=1:n
        if (L(j,j)==0) error('Matrix is singular!'); end;
        x(j)=b(j)/L(j,j);
        b(j+1:n)=b(j+1:n)-L(j+1:n,j)*x(j);
    end
end
