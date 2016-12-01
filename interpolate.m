function [ p, T ] = interpolate( X, Y, Yprime )
%INTERPOLATE
%
%   [p, T] = interpolate(X, Y, Yprime)
%
%   Uses the method of divided differences to construct a polynomial
%   that interpolates the given data points. If a Yprime vector is 
%   given, then a hermite interpolating polynomial is constructed.
%
%PARAMETERS:
%   X       Vector of the x-cooridinates.
%
%   Y       Vector of the y-cooridinates.
%
%   Yprime  An optional vector of the slope at each (x,y) datapoint.
%
%RETURNS:
%   p       A vector of the coefficients of the interpolating 
%           polynomial (standard basis).
%
%   T       The divided difference table used to construct the polynomial.
%
%NOTES:
%   The vectors X, Y, and Yprime must all be row vectors or 
%   column vectors. They cannot be a mix of row and column vectors.
%
%AUTHOR:    Kenneth Poppie
%DATE:      Oct 26, 2016

% See if we are dealing with hermite data.
if nargin < 3
    hermite = 0;
else
    hermite = 1;
end

% Make sure size of vectors are the same.
if (size(X,1) ~= size(Y,1)) || (size(X,2) ~= size(Y,2))
    error('Vectors X and Y must have same dimensions.');
end
if hermite
    if (size(Yprime,1) ~= size(Y,1)) || (size(Yprime,2) ~= size(Y,2))
        error('Vectors X, Y, and Yprime must have same dimensions.');
    end
end

% How many data points do we have?
n = max(size(X,1),size(X,2));

if ~hermite 
    % Initialize the DD table.
    % First column of table is the Y vector.
    T = zeros(max(size(X,1),size(X,2)));
    if size(Y,2) > size(Y,1)
        T(:,1) = Y';
    else
        T(:,1) = Y;
    end
    
    % Start filling in the entries of the DD table.
    d = 1;
    while d < n
        row = 1;
        while row <= (n-d)
            T(row,d+1) = (T(row+1,d) - T(row,d))/(X(row+d)-X(row));
            row = row+1;
        end
        d = d+1;
    end
    
    % Now that we have the DD table, we can contruct the polynomial.
    p = polyfromtable(T(1,1:end),X);
    
else % hermite data
    % Initialize the hermite DD table
    T = zeros(2*n, 2*n+1);
    
    % Fill in first and second column
    toff = 1;
    ii = 1;
    while ii <= n
        T(toff:toff+1,1:2) = [X(ii), Y(ii); X(ii), Y(ii)];
        ii = ii + 1;
        toff = toff + 2;
    end
    % Fill in the third column
    row = 1;
    ii = 1;
    while row < 2*n
        if T(row,1) == T(row+1,1)
            T(row,3) = Yprime(ii);
            ii = ii + 1;
        else 
            T(row,3) = (T(row+1,2)-T(row,2))/(T(row+1,1)-T(row,1));
        end
        row = row + 1;
    end
    % Now fill in the rest of the table
    col = 4;
    d = 2;
    while col <= 2*n+1
        row = 1;
        while row <= (2*n-d) 
            T(row,col) = (T(row+1,col-1)-T(row,col-1))/...
                (T(row+d,1)-T(row,1));
            row = row + 1;
        end
        col = col + 1;
        d = d + 1;
    end
    
    % Now that we have the hermite DD table, we can build the 
    % polynomial 
    p = polyfromhermitetable(T(1,2:end), X);
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%
function p = polyfromtable( T, x )
    z = size(T,2)-1;
    p = [zeros(1,z), T(1,1)];
    z = z - 1;

    tt = 2;
    while tt <= size(T,2)
        % Compute term
        ii = 1; v = 1;
        while ii < tt
            v = conv([1, -x(ii)], v);
            ii = ii + 1;
        end
        p = p + [zeros(1,z), T(1,tt)*v];
        z = z - 1;
        tt = tt + 1;
    end
end % function

function p = polyfromhermitetable( T, x ) 
    z = size(T,2)-1; 
    p = [zeros(1,z), T(1,1)];
    z = z - 1;
    tt = 2;
    m = fix(tt/2);
    w = 1; 
    while tt <= size(T,2)
        w = conv([1,-x(m)],w);
        p = p + [zeros(1,z), w*T(1,tt)];
        z = z - 1;
        tt = tt + 1;
        m = fix(tt/2);
    end
end % function