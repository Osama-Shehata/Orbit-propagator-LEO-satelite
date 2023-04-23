
function [ra dec] = ra_and_dec_from_r(r)
l = r(1)/norm(r);
m = r(2)/norm(r);
n = r(3)/norm(r);
dec = asind(n);
if m > 0
ra = acosd(l/cosd(dec));
else
ra = 360 - acosd(l/cosd(dec));
end
