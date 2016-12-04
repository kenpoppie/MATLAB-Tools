function [ f, p, rsq ] = regression( x, y, model )
%REGRESSION
%
%   [ f, rsq ] = regression( x, y, model )
%
%   Performs regression analysis on the given (x,y) data fitting 
%   it to the model specified (see parameters for list of models).
%
%PARAMETERS:
%   x       The x values of the data set.
%
%   y       The y values of the data set.
%   
%   model   The model the data should be fitted to.
%           'linear'        y = mx + b
%           'power'         y = ax^b
%           'expo'          y = ab^x
%           'quad'          y = ax^2 + bx + c
%           'cubic'         y = ax^3 + bx^2 + cx + d
%
%RETURNS:
%   f       The regression function as a vectorized string.
%
%   p       A vector of the computed paramters.
%           For 'linear', [m,b] is returned.
%           For 'power' and 'expo', [a,b] is returned.
%           For 'quad', [a,b,c] is returned. 
%           For 'cubic', [a,b,c,d] is returned.
%
%   rsq     R-squared value.
%           R-squared has the useful property that its scale is intuitive:
%           it ranges from zero to one, with zero indicating that the 
%           proposed model does not improve prediction over the mean model 
%           and one indicating perfect prediction. Improvement in the 
%           regression model results in proportional increases in
%           R-squared.
%
%AUTHOR:    Kenneth Poppie
%DATE:      Nov. 13, 2016

% Initialize return value
f = '';
p = [];

% Make sure X and Y are column vectors.
if size(x,1) < size(x,2)
    X = x';
else
    X = x;
end
if size(y,1) < size(y,2)
    Y = y';
else
    Y = y';
end

% Get size of dataset
nX = max(size(X,1),size(X,2));
nY = max(size(Y,1),size(Y,2));
if nX ~= nY
    error('X and Y not the same dimensions.');
end

if strcmpi(model,'linear')
    % Linear regression 
    % y = mx + b
    
    A = [X ones(nX,1)];
    
    % Compute m and b
    S = (A'*A)\(A'*Y);
    f = sprintf('%0.4f * x + %0.4f', S(1),S(2)); 
    p = [S(1),S(2)];
    
elseif strcmpi(model,'power')
    % Power law 
    % y = ax^b
    
    C = [log(X) ones(nX,1)];
    d = log(Y);
    S = (C'*C)\(C'*d);
    a = exp(S(2)); b = S(1);
    f = sprintf('%0.4f*x^%0.4f', a, b);
    p = [a,b];
    
elseif strcmpi(model,'expo')
    % Exponential law
    % y = ab^x
    
    C = [X ones(nX,1)];
    d = log(Y);
    S = (C'*C)\(C'*d);
    b = exp(S(1));  a = exp(S(2));
    f = sprintf('%0.4f*%0.4f^x', a, b); 
    p = [a,b];
    
elseif strcmpi(model,'quad')
    % Quadratic
    % y = ax^2 + bx + c
    
    C = [X.^2  X  ones(nX,1)];
    S = (C'*C)\(C'*Y);
    f = sprintf('%0.4f*x^2 + %0.4f*x + %0.4f', S(1),S(2),S(3));
    p = [S(1),S(2),S(3)];

elseif strcmpi(model,'cubic')
    % Cubic 
    % y = ax^3 + bx^2 + cx + d
    
    C = [X.^3  X.^2  X  ones(nX,1)];
    S = (C'*C)\(C'*Y);
    f = sprintf('%0.4f*x^3 + %0.4f*x^2 + %0.4f*x + %0.4f', ...
        S(1),S(2),S(3),S(4));
    p = [S(1),S(2),S(3),S(4)];
    
else
    error('%s is not a valid regression model.', model);
end

f = vectorize(f);

% Compute R-squared
ybar = sum(y)/nY;
sse = sum((y - eval(f)).^2);
tss = sum((y-ybar).^2);
rsq = 1 - sse/tss;

end


