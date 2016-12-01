function pstr = poly2str( ply , varletter)
%POLY2STR
%
%   pstr = poly2str( ply, varletter )
%
%PARAMTERS:
%   ply         Vector of the polynomial's coefficient.
% 
%   varletter   The variable letter to use in string. The 
%               default is 'x'.
%
%RETURNS:
%   pstr        The polynomial as a vectorized string that can be 
%               evaulated using MATLAB's eval() function.
%
%EXAMPLE:
%   poly2str([2,3,5],'t') returns '+2.*t.^2+3.*t.^1+5'.
% 
%AUTHOR:    Kenneth Poppie
%DATE:      Nov. 23, 2016

if nargin > 1
    vl = varletter(1);
else 
    vl = 'x';
end
n = length(ply);
d = n-1;
pstr = '';
for ii = 1:n-1
    if ply(ii) ~= 0
        pstr = sprintf('%s%+g*%c^%d',pstr,ply(ii),vl,d);
    end
    d = d - 1;
end
pstr = sprintf('%s%+g',pstr,ply(ii+1));
pstr = vectorize(pstr);

end % function

