function [A,T,SSqDev] = nl_regression(x, y, model, parts, A, tol, max_iter)
%NL_REGRESSION
%   
%   [A,T,SSqDev] = nl_regression(x, y, model, parts, A, tol, max_iter)
%
%   Nonlinear regression analysis.
%   Attempts to compute the model paramters that give the best fit 
%   to the given data.  Uses the Gauss-Newton algorithm.
%
%PARAMETERS:
%   x           The x-coordinates of the data.
%
%   y           The y-coordinates of the data.
%
%   model       The nonlinear model f(x,A).
%
%   parts       A cell array of the partial derivatives of the model with
%               respect to each of the unknown parameters in order.
%
%   A           The initial guess for the parameters A(1),A(2),...,A(n).
%
%   tol         The tolerance. If zero, then the maximum number of 
%               iterations is performed.
%
%   max_iter    The maximum number of iterations that will be done if 
%               tolerance is not met.
%
%RETURNS:
%   A           The final approximation of the parameters that give the 
%               best fit to the data set.
%
%   T           A table of the approximations at each iteration. Use when
%               checking for convergence.
%
%   SSqDev      A table of the sum of the square of the deviations for 
%               each iteration.
%
%
%EXAMPLE:
%--------------------------------------------------------------------
%   model = @(x,A) A(1)*(1 - exp(-A(2)*x));
%   partials = { ...
%       @(x,A) 1-exp(-A(2)*x), ...      % with respect to A(1)
%       @(x,A) A(1)*x*exp(-A(2)*x) };   % with respect to A(2)
% 
%   xData = [0.25, 0.75, 1.25, 1.75, 2.25];
%   yData = [0.28, 0.57, 0.68, 0.74, 0.79];
% 
%   A = [1, 1]; % initial guess for A(1) and A(2).
%   tolerance = 1e-5;
%   max_iter = 50;
% 
%   [A,T,sqdev] = nl_regression(xData,yData,...
%                               model,partials,A,tolerance,max_iter)
%--------------------------------------------------------------------
% Thus, A(1) and A(2) are the values for the model's parameters
% that give the best fit to the data.
%
%NOTES:
%   This method is NOT guaranteed to converge even for a "good" initial 
%   guess.
%
%AUTHOR:    Kenneth Poppie
%DATE:      Nov. 16, 2016

% More sure A is a column vector.
if size(A,1) < size(A,2) 
    A = A';
end

% Check sizes
nX = length(x); 
nY = length(y); 
nA = length(A); 
nP = length(parts);
if nX ~= nY 
    error('x and y coordinates not matching up.');
end
if nA ~= nP
    error('A and partial derviatives not matching up.');
end

% Initialize
J = zeros(nX,nP);
b = zeros(nX,1);
SSqDev = [];
T = A';

% Begin iterations
iter = 1;
while iter <= max_iter 
    % Build Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 1;
    while ii <= nX 
        jj = 1;
        while jj <= nP
            df = parts{jj};
            J(ii,jj) = df(x(ii),A);
            jj = jj + 1;
        end
        ii = ii + 1;
    end
    
    % Build b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 1;
    while ii <= nX
        b(ii,1) = y(ii) - model(x(ii),A);
        ii = ii + 1;
    end
    SSqDev = [ SSqDev; sum(b(:,1).^2) ];
    
    % Approximation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oldA = A;
    A = A + (J'*J)\(J'*b);
    T = [T ; A'];
    
    % Check tolerance %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tol_met = 1; % assume true
    ii = 1;
    while ii < nA
        if abs((A(ii) - oldA(ii))/A(ii)) >= tol
            tol_met = 0;
            break;
        end
        ii = ii + 1;
    end
    if tol_met 
        % stop iterations
        break;
    end
    
    % Next iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter = iter + 1;
end

% Compute last deviation %%%%%%%%%%%%%%%%%%%%%%%%
ii = 1;
while ii <= nX
    b(ii,1) = y(ii) - model(x(ii),A);
    ii = ii + 1;
end
SSqDev = [ SSqDev; sum(b(:,1).^2) ];

end

