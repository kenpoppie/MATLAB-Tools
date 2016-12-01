function S = spline2str( x, coeff, varletter )
%SPLINE2STR
% 
%   S = spline2str( x, coeff, varletter )
%
%   Builds a cell array of polynomials for each of the spline 
%   segments given the x-coordinates and spline coefficients.
%   The polynomials are vectorized.
%
%PARAMTERS:
%   x           A vector of the x-coordinates used to construct the 
%               spline.
%
%   coeff       A matrix of the spline's coefficients.
%
%   varletter   (Optional) The variable letter to be used in the 
%               polynomials. The default is 'x'.
%
%RETURNS:
%   S           A cell array of vectorized polynomials.
%
%AUTHOR:        Kenneth Poppie
%DATE:          Nov. 26, 2016

if nargin > 2
    vl = varletter(1);
else 
    vl = 'x';
end

S = cell(size(coeff,1),1);
for ii = 1:size(coeff,1)
    p = '';
    d = size(coeff,2)-1;
    for jj = 1:size(coeff,2)
        if coeff(ii,jj) ~= 0
            if d > 0
                p = sprintf('%s%+g*(%c%+g)^%d', ...
                        p, coeff(ii,jj), vl, -x(ii), d);
            else
                p = sprintf('%s%+g', p, coeff(ii,jj));
            end
        end
        d = d - 1;
    end
    S{ii} = vectorize(p);
end

end % function

