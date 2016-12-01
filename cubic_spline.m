function coeff = cubic_spline( x, y, spline_type, slope )
%CUBIC_SPLINE
%   
%   coeff = cubic_spline( x, y, spline_type, slope )
%
%   Creates a cubic spline for the given (x,y) dataset. 
%   Uses code taken from Dr. David Hill's nspline, cspline, 
%   and nakspline functions. Email: dhill001@temple.edu
%
%PARAMETERS:
%   x               A vector of x-coordinates for the dataset.
%                   Must be in increasing order!!!
%
%   y               A vector of y-coordinates for the dataset.
%
%   spline_type     The type of spline to be generated.  The options
%                   are: 
%                   'natural', or 'nat': generate a natural spline.
%                       S''(x_0) = 0 and S''(x_n) = 0
% 
%                   'clamped', or 'clamp': generate a clamped spline.
%                       S'(x_0) = a and S'(x_n) = b where a,b are 
%                       specified by the slope paramter.
%
%                   'notaknot', or 'nak': generate a not-a-knot spline.
%                       S''' is continuous at x = x_1 and x = x_n-1
%
%   slope           The vector [a,b] where a,b are slopes for a clamped
%                   spline.  Ignored for the other spline types.
%
%RETURNS:
%   coeff           A matrix of the spline coefficients.
% 
%AUTHOR:    Kenneth Poppie
%DATE:      Nov. 21, 2016

% Check size of dataset.
n=length(x);
m=length(y);
if n~=m
   error('Input vectors of different lengths.');
end

% Make sure x-coordinates are in increasing order.
if any(sort(x)~=x)
    error('Independent variable must be in increasing order.');
end

% Make sure x and y are column vectors.  Make them column vectors
% if they're not.
[xm,xn]=size(x);
if xm~=n || xn~=1
    x=x';
end
[ym,yn]=size(y);
if ym~=m || yn~=1
    y=y';
end

if strcmpi(spline_type,'natural') || strcmpi(spline_type,'nat')
    
    % Calculate the Natural Spline %%%%%%%%%%%%%%%%%%
    
    nm1 = n-1;
    for j=1:nm1
        h(j)=x(j+1)-x(j);
    end
    
    d=[1]; d1=[0]; d2=[];
    
    for j=1:n-2
        d=[d 2*(h(j)+h(j+1))];
        d1=[d1 h(j+1)];
        d2=[d2 h(j)];
    end
    
    d=[d 1]; d2=[d2 0];
    A=diag(d)+diag(d1,1)+diag(d2,-1);b=zeros(n,1);
    
    for j=2:n-1
      b(j)=(3/h(j))*(y(j+1)-y(j)) -(3/h(j-1))*(y(j)-y(j-1));
    end
    C=A\b;
    for j=nm1:-1:1
      B(j)=(y(j+1)-y(j))/h(j) - h(j)*(C(j+1)+2*C(j))/3;
      D(j)=(C(j+1)-C(j))/(3*h(j)); 
    end
    
    coeff=[D(1:nm1)' C(1:nm1) B(1:nm1)' y(1:nm1)];

elseif strcmpi(spline_type,'clamped') || strcmpi(spline_type,'clamp')
    
    % Calculate the Clamped Cubic Spline %%%%%%%%%%%%%%%
    
    % Get the slope.
    if nargin < 4
        error('slope must be given for clamped cubic spline.');
    end
    v1 = slope(1);  v2 = slope(2);
    
    nm1=n-1;
    for j=1:nm1
       h(j)=x(j+1)-x(j);
    end

    d=[2*h(1)]; d1=[h(1)]; d2=[];
    for j=1:n-2
       d=[d 2*(h(j)+h(j+1))];
       d1=[d1 h(j+1)];
       d2=[d2 h(j)];
    end

    d=[d 2*h(nm1)]; d2=[d2 h(nm1)];

    A=diag(d)+diag(d1,1)+diag(d2,-1);
    b=zeros(n,1);
    b(1)=(3/h(1))*(y(2)-y(1))-3*v1;

    for j=2:n-1
      b(j)=(3/h(j))*(y(j+1)-y(j)) -(3/h(j-1))*(y(j)-y(j-1));
    end

    b(n) = 3*v2-(3/h(n-1))*(y(n)-y(n-1));
    C=A\b;

    for j=nm1:-1:1
      B(j)=(y(j+1)-y(j))/h(j) - h(j)*(C(j+1)+2*C(j))/3;
      D(j)=(C(j+1)-C(j))/(3*h(j)); 
    end
    coeff=[D(1:nm1)' C(1:nm1) B(1:nm1)' y(1:nm1)];

elseif strcmpi(spline_type,'notaknot') || strcmpi(spline_type,'nak')
    
    % Calculate the Not-A-Knot Cubic Spline. %%%%%%%%%%%%%%%%
    nm1=n-1;
    
    for i = 1:n-1
        hi(i) = x(i+1) - x(i);
    end
    
    for i = 1:n-2
        dd(i) = 2.0*(hi(i) + hi(i+1));
        ri(i) = (3.0/hi(i+1))*(y(i+2)-y(i+1))-...
                                    (3.0/hi(i))*(y(i+1)-y(i));
    end
    
    dd(1)   = dd(1)   + hi(1) + hi(1)^2 / hi(2);
    dd(n-2) = dd(n-2) + hi(n-1) + hi(n-1)^2 / hi(n-2);

    du = hi(2:n-2);
    dl = du;
    du(1) = du(1) - hi(1)^2 / hi(2);
    dl(n-3) = dl(n-3) - hi(n-1)^2 / hi(n-2);

    temp = tridiagonal ( dl, dd, du, ri ); %Solving the system

    c = zeros ( n,1 );
    d = c;   b = c;

    c(2:n-1) = temp;
    c(1) = ( 1 + hi(1) / hi(2) ) * c(2) - hi(1) / hi(2) * c(3);
    c(n) = ( 1 + hi(n-1) / hi(n-2) ) * c(n-1) - hi(n-1) / hi(n-2) * c(n-2);
    for i = 1 : n-1
        d(i) = (c(i+1)-c(i))/(3.0*hi(i));
        b(i) = (y(i+1)-y(i))/hi(i) - hi(i)*(c(i+1)+2.0*c(i))/3.0;
    end
    
    coeff=[d(1:nm1) c(1:nm1) b(1:nm1) y(1:nm1)];

else 
    error('%s is not a valid spline type.\n%s', ...
        spline_type, ...
        'Valid types are ''natural'', ''clamped'', and ''notaknot''');
end

end % function


function y = tridiagonal ( c, a, b, r ) %Bradie's routine

%TRIDIAGONAL  solve a linear system with a tridiagonal coefficient matrix
%
%     calling sequence:
%             x = tridiagonal ( c, a, b, r )
%             tridiagonal ( c, a, b, r )
%
%     inputs:
%             c       vector containing the entries along lower diagonal 
%                     of the coefficient matrix
%             a       vector containing the entries along main diagonal 
%                     of the coefficient matrix
%             b       vector containing the entries along upper diagonal 
%                     of the coefficient matrix
%             r       right-hand side vector
%
%     output:
%             x       solution vector
%
%     NOTE:
%             the entries in the vectors c, a and b are assumed to be
%             numbered as follows:
%
%                 | a(1)  b(1)                                 |
%                 | c(1)  a(2)  b(2)                           |
%                 |       c(2)  a(3)  b(3)                     |
%                 | 	         .     .     .                 |
%                 | 			       .     .     .           |
%                 | 				         .     .    b(n-1) |
%                 | 					        c(n-1)  a(n)   |
%

n = length ( a );

%
%   factorization step
%

for i = 1 : n-1
    b(i) = b(i) / a(i);
	a(i+1) = a(i+1) - c(i) * b(i);
end

%
%   forward substitution
%

r(1) = r(1) / a(1);
for i = 2 : n
    r(i) = ( r(i) - c(i-1) * r(i-1) ) / a(i);
end

%
%   back substitution
%

for i = n-1 : -1 : 1
    r(i) = r(i) - r(i+1) * b(i);
end

if ( nargout == 0 )
   disp ( r )
else
   y = r;
end
end % function

