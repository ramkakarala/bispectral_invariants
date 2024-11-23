function Y = spharm(ell,alt,azi)
% usage
%        Y = spharm(ell,alt,azi)
% computes the spherical harmonic 
%     Yelln = sqrt((ell+1/2)(ell-n)!/(ell+n)!)*Pnm(cos(alt))*exp(j*n*azi)
% where Pnm is the associated Legendre function and
% 0 <= alt <= pi and azi is in 0,2*pi
% r kakarala
x = cos(alt);
p = legendre(ell,x,'norm');  % see matlab for this normalization
P = [p(end:-1:2)' p(:)'];    % make one row vector of the valuse
Y = P.*exp(-j*(-ell:ell)*azi);


