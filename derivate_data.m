function dy = derivate_data( x, y )
%DERIVATE_DATA
%
%   dy = derivate_data( x, y )
%
%   Given a set of data points, this function approximates the
%   first derivative (slope of the tangent line) at each point.
%
%PARAMTERS:
%   x       The x-coordinates of the data set.
%
%   y       The y-coordinates of the data set.
%
%RETURNS:
%   dy      A vector of approximate derivatives.
%
%NOTES:
%   The x-coordinates must be equispaced (same distance 
%   from one another). For example x = [1.0, 1.2, 1.4, 1.6] is 
%   equispaced, but x = [1, 1.3, 1.9, 2.1] is not.
%
%AUTHOR: Kenneth Poppie
%DATE:   September 14, 2016

% Make sure we have at least two data points.
xn = length(x);
yn = length(y);
if xn ~= yn 
    error('Not equal number of data points.');
end
if xn < 2
    error('Not enough data points. Minimum is two.');
end

% Make sure the data is equispaced.
h0 = x(2)-x(1);
for ii = 2:xn-1
    h1 = x(ii+1)-x(ii);
    if h0 ~= h1 
        error('Data is not equispaced.');
    end
    h0 = h1;
end

% Initialize the answer.
dy = zeros(1,yn);

% If there are only two data points, then we use the first order 
% formulas for approximating the derivative at each point.
% h0 is the stepsize.
if xn == 2
    % Forward difference approx.
    dy(1) = (y(2)-y(1))/h0;
    
    % Backward difference approx.
    dy(2) = (y(2)-y(1))/h0;
    
else
% There are three or more points, so we can use second order formulas
% for making the approximations. 
% h0 is the stepsize.

% For the first point, we use the forward difference.
    dy(1) = (-3*y(1)+4*y(2)-y(3))/(2*h0);
    
% For the last point, we use the backwards difference.
    dy(yn) = (3*y(yn)-4*y(yn-1)+y(yn-2))/(2*h0);
    
% Now for the points in between, we use the centered difference.
    for ii = 2:yn-1
        dy(ii) = (y(ii+1)-y(ii-1))/(2*h0);
    end
end

end % function