function [ p, p_n ] = findroot_secant( f, tol, max_iter, p_0, p_1 )
%FINDROOT_SECANT
%
%   [p, p_n] = findroot_secant(f, tol, max_iter, p_0, p_1)
%
%   Uses the secant method to approximate the root of the given 
%   function f with initial guesses p_0 and p_1.
%
% PARAMETERS:
%   f           The function f(x) to approximate the root for.
%               This can be either a function handle, 
%               string representation of the function, or a vector
%               representation of a polynomial.
%
%   tol         The tolerance.  If tolerance is zero or less, the 
%               maximum number of iterations will be performed.
%
%   max_iter    The maximum number of iterations the method will do.
%
%   p_0         First initial guess (lower bound). 
%
%   p_1         Second initial guess (upper bound).
%
% RETURNS:
%   p           The approximate root for f(x).
%
%   p_n         A vector containing the root approximations at each 
%               iteration.
%
%AUTHOR: Kenneth Poppie
%DATE:   September 14, 2016  

% Check number of iterations
if max_iter < 1
    error('The number of iterations, max_iter, must be at least 1');
end

% Initialize the return values
p = 0;
p_n = [p_0; p_1];

% Begin secant method
p_n1 = p_1;
p_n2 = p_0;
for n = 1:1:max_iter
    if isa(f,'function_handle')
        f_pn1 = f(p_n1);
        f_pn2 = f(p_n2);
    elseif isa(f,'char')
        x = p_n1;
        f_pn1 = eval(f);
        x = p_n2;
        f_pn2 = eval(f);
    elseif isa(f,'double') && isvector(f)
        x = p_n1;
        f_pn1 = polyval(f,p_n1);
        f_pn2 = polyval(f,p_n2);
    else
        error('f parameter must be a %s', ...
            'function handle, string, or vector.');
    end
        
    p = p_n1 - f_pn1*(p_n1 - p_n2)/(f_pn1 - f_pn2);
    if isnan(p)
        % we got 0/0
        warning('NaN computed. Aborting algorithm.');
        p = p_n1;
        break;
    end
    
    % Add approx. root to vector.
    p_n = [p_n; p];
    
    % Check tolerance.
    if abs(p - p_n1) < tol 
        break; % tolerance met! Done.
    end
    
    % Next iteration
    p_n2 = p_n1;
    p_n1 = p;
end
end % function

