function spline_plot( X, S )
%SPLINE_PLOT
%
%   spline_plot( X, S ) 
% 
%   Plots each of the spline segments over the range X(0) to X(n).
%
%PARAMETERS:
%   X       Vector of x-coordinates for the dataset.
%
%   S       Cell array of the spline polynomials or 
%           a matrix of the spline's coefficients.
%
%EXAMPLE:
%-----------------------------------------------------------------
%   N = [0.0521, 0.1028, 0.2036, 0.4946, 0.9863, 2.443, 5.06 ];
%   D = [1.65, 2.10, 2.27, 2.76, 3.12, 2.92, 2.07 ];
%
%   S = cubic_spline(N,D,'notaknot');
%   figure
%   hold on
%   spline_plot(N,S)
%   scatter(N,D,40)
%   grid on
%   hold off
%-----------------------------------------------------------------
%
%AUTHOR:    Kenneth Poppie
%DATE:      Nov. 21, 2016

% Size of dataset.
n=length(X);

% What is S?
if isa(S,'double') % S is a matrix of coefficients
    W = spline2str(X, S);
elseif isa(S, 'cell') 
    W = S;
else
    error('S is invalid data type.');
end

% Plot spline segments
plot_points = 100; 
hold on
for ii = 1:n-1
    x = X(ii):(X(ii+1)-X(ii))/plot_points:X(ii+1); 
    y = eval(W{ii});
    plot(x,y,'LineWidth',2);
end
    
end % function

