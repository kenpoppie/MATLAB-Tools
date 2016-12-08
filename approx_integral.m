function area = approx_integral( f, a, b, rule, n )
%APPROX_INTEGRAL
%
%   area = approx_integral( f, a, b, rule, n )
%
%   Approximates the area under the curve defined by the function f 
%   over the interval [a,b].
%
%PARAMTERS:
%   f       The function (as a string) to approximate the integral of. 
%   
%   a       The x-coordinate of the left boundry of the area.
%
%   b       The x-coordinate of the right boundry of the area.
%
%   rule    Determines the method to use in computing the approximation.
%           The choices are:
%               'trap'      for composite trapezoidal rule,
%               'simp'      for composite Simpson's 1/3 rule, 
%               'mid'       for composite midpoint rule. 
% 
%   n       The number of subintervals to use. Must be 1 or greater.
%
%RETURNS:
%   area    The approximate area under the curve for the interval [a,b].
%
%AUTHOR:    Kenneth Poppie
%DATE:      Dec. 4, 2016

% Make sure f is a function or a string.
if ~isa(f,'function_handle') && ~isa(f,'char')
    error('f must be a function handle or string to evalulate.');
end

% Check n
if n < 1
    error('Number of subinverals must be 1 or greater.');
end

% Approximate area.
area = 0;
if strcmpi(rule,'trap')
    
    % Composite trapezoidal rule %%%%%%%%%
    
    h = (b-a)/n;
    x = a:h:b;
    x = eval(vectorize(f));
    area = (h/2)*(x(1)+x(n+1)+2*sum(x(2:n)));

elseif strcmpi(rule,'simp')
    
    % Composite Simpson's rule %%%%%%%%%
    
    h=(b-a)/(2*n);
    x = a:h:b;
    y = eval(vectorize(f));
    t = [zeros(1,n-1); ones(1,n-1)];
    t = t(:)';
    t = 2*(ones(size(t))+t);
    t = [1 4 t 1]; %array of coefficients for Simpson's Rule
    area = (h/3)*(t*y');
    
elseif strcmpi(rule,'mid')
    
    % Composite midpoint rule %%%%%%%%%%%%%%
    
    h=(b-a)/(2*n);
    xp = a:h:b;
    x = xp(2:2:length(xp));
    area = 2*h*sum(eval(vectorize(f)));

else
    error('Not a valid rule. Type ''help approx_integral'' %s', ...
        'to see rule choices.');
      
end % if

end % function

