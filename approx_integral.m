function area = approx_integral( f, a, b, rule, n )
%APPROX_INTEGRAL
%
%   area = approx_integral( f, a, b, rule, n )
%
%   Approximates the area under the curve defined by the function f 
%   over the interval [a,b].
%
%PARAMTERS:
%   f       The function to approximate the integral of. 
%           This can be a function handle or a string representation 
%           of the function.
%   
%   a       The x-coordinate of the left boundry of the area.
%
%   b       The x-coordinate of the right boundry of the area.
%
%   rule    Determines the method to use in computing the approximation.
%           The choices are:
%               'trap'      for trapezoidal rule,
%               'simp1/3'   for Simpson's 1/3 rule, 
%               'simp3/8'   for Simpson's 3/8 rule,
%               'mid'       for midpoint rule. 
% 
%   n       (optional) If n is provided and n > 2, the composite 
%           version of the trapezoidal rule or Simpson's 1/3 rule
%           is performed.  For the other rules, this parameter is 
%           ignored.
%
%RETURNS:
%   area    The approximate area under the curve for the interval [a,b].
%
%AUTHOR:    Kenneth Poppie
%DATE:      Dec. 4, 2016

% Check a and b. 
if a > b 
    error('a must be less than b.');
end

% Make sure f is a function or a string.
if ~isa(f,'function_handle') && ~isa(f,'char')
    error('f must be a function handle or string to evalulate.');
end

% See if we are using the composite version (if it applies).
using_composite = 0;
if (nargin > 4) 
    if n > 2
        using_composite = 1;
    else 
        error('Invalid number of subintervals.');
    end
end

% Approximate area.
area = 0;
if strcmpi(rule,'trap')
    
    % trapezoidal rule %%%%%%%%%
    if ~using_composite
        h = b - a;
        if isa(f,'function_handle')
            w1 = f(a);
            w2 = f(b);
        else
            x = a;
            w1 = eval(f);
            x = b;
            w2 = eval(f); 
        end
        area = area + (h/2)*(w1 + w2);
        
    else % using composite rule
        h = (b-a)/n;
        if isa(f,'function_handle')
            % Eval endpoints
            w_0 = f(a);
            w_n = f(b);
            
            % Eval midpoints
            w_mid = 0;
            for ii = 1:(n-1)
                w_mid = w_mid + 2*f(a+ii*h);
            end
            
        else % string
            % Eval endpoints
            x = a; w_0 = eval(f);
            x = b; w_n = eval(f);
            
            % Eval midpoints
            w_mid = 0;
            for ii = 1:(n-1)
                x = a+ii*h;
                w_mid = w_mid + 2*eval(f);
            end
        end
                    
        % Compute area
        area = (h/2)*(w_0 + w_mid + w_n);
    end

elseif strcmpi(rule,'simp1/3')
    
    % Simpson's 1/3 rule %%%%%%%%%%%
    
    if ~using_composite
        h = (b-a)/2;
        if isa(f,'function_handle')
            w1 = f(a);
            w2 = f((a+b)/2);
            w3 = f(b);
        else
            x = a;
            w1 = eval(f);
            x = (a+b)/2;
            w2 = eval(f);
            x = b;
            w3 = eval(f);
        end
        area = (h/3)*(w1 + 4*w2 + w3);
        
    else % using composite rule
        h = (b-a)/n;
        if isa(f,'function_handle')
            % Eval endpoints
            w_0 = f(a);  w_n = f(b);
            
            % Eval midpoints
            w_mid = 0;
            for ii = 1:(n-1)
                x = a + ii*h;
                if mod(ii,2) == 1 % odd
                    coef = 4;
                else
                    coef = 2;
                end
                w_mid = w_mid + coef*f(x);
            end

        else % string
            % Eval endpoints
            x = a;
            w_0 = eval(f);
            x = b; 
            w_n = eval(f);
            
            % Eval midpoints
            w_mid = 0;
            for ii = 1:(n-1)
                if mod(ii,2) == 1 % odd
                    coef = 4;
                else
                    coef = 2;
                end
                x = a + ii*h;
                w_mid = w_mid + coef*eval(f);
            end
        end
        
        % Compute area
        area = (h/3)*(w_0 + w_mid + w_n);
        
    end % if 
    
elseif strcmpi(rule,'simp3/8')
    
    % Simpson's 3/8 rule %%%%%%%%%%%
    
    h = (b-a)/3;
    if isa(f,'function_handle')
        w1 = f(a);
        w2 = f((2*a+b)/3);
        w3 = f((a+2*b)/3);
        w4 = f(b);
    else
        x = a;
        w1 = eval(f);
        x = (2*a+b)/3;
        w2 = eval(f);
        x = (a+2*b)/3;
        w3 = eval(f);
        x = b;
        w4 = eval(f);
    end
    area = 3*(h/8)*(w1 + 3*w2 + 3*w3 + w4);
    
elseif strcmpi(rule,'mid')
    
    % Midpoint rule %%%%%%%%%%
    
    h = b - a;
    if isa(f,'function_handle')
        area = h*f((a+b)/2);
    else
        x = (a+b)/2;
        area = h*eval(f);
    end
    
else
    error('Not a valid rule. Type ''help approx_integral'' %s', ...
        'to see rule choices.');
      
end % if

end % function

