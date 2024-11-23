n=input('n=');
% compute scale factor to apply later
alpha = sqrt(2*n+1)/factorial(n);
alpha = alpha/2^n;
nmono = (n+1)*(n+2)/2;
%%
syms x y z t;
P = expand((x-i*y + 2*z*t - (x + i*y)*t^2)^n);
Pt = collect(P,t);
%% what monomials are there?
lmn = 0;
clear Mono;
exponents = zeros(nmono,3);
for l = n:-1:0
    for m = 0:n-l
        lmn = lmn + 1;
        p = n-l-m;
        Mono(lmn)=x^l * y^m * z^p;
        exponents(lmn,:)=[l,m,p];
    end;
end;
%% use MATLAB 2009a to extract spherical harmonic terms (unnormalized)
[cterms,tterms] = coeffs(Pt,'t');  
clear Y;
for tt=1:2*n+1
    Y(tt,1)=cterms(tt);
end;
%% get the coeffs
An=zeros(2*n+1,lmn); % without the scale factors
Af=zeros(2*n+1,1);
for tt=1:2*n+1
    U=Y(tt);
    for mm = 1:lmn
        dz = exponents(mm,3); dy = exponents(mm,2); dx = exponents(mm,1);
        normfactor = factorial(dz)*factorial(dy)*factorial(dx);
        dU = diff(diff(diff(U,'z',dz),'y',dy),'x',dx);
        An(tt,mm) = subs(subs(subs(dU,'z',0),'y',0),'x',0)/normfactor;
    end;
    m = tt-1-n;
    Af(tt)=alpha*sqrt(factorial(n-m)*factorial(n+m));
end;
%%
A = diag(Af)*An;
disp(A*conj(Mono'));  % show the homogeneous polynomials

      
