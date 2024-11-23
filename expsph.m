%% make spherical harmonics Y^k and Y^ell
ell = input('ell='); k = input('k=');
if (ell < k)
    fprintf(1,'reversing ell and k')
    temp = k;
    k = ell;
    ell = temp;
end;

D = (2*ell + 1) * (2*k+1);
%% initialize
Ytensor1 = zeros(D);  Ytensor2 = Ytensor1;
Ydsum1 = zeros(D);    Ydsum2 = Ydsum1;
Ckl = zeros(D);
%% 
for dd = 1:D
    alt = pi*rand; azi = 2*pi*rand;
    Yell = spharm(ell,alt,azi);  
    Yk = spharm(k,alt,azi);
    Ytensor1(dd,:)=kron(Yell,Yk);
    index = 0;
    for p = ell-k : ell + k
        Ydsum1(dd,index+1:index+(2*p+1)) = spharm(p,alt,azi);
        index = index + 2*p + 1;
    end;
end;    
%% repeat
for dd = 1:D
    alt = pi*rand; azi = 2*pi*rand;
    Yell = spharm(ell,alt,azi);  
    Yk = spharm(k,alt,azi);
    Ytensor2(dd,:)=kron(Yell,Yk);
    index = 0;
    for p = ell-k : ell + k
        Ydsum2(dd,index+1:index+(2*p+1)) = spharm(p,alt,azi);
        index = index + 2*p + 1;
    end;
end;    
%% find the Clebsh Gordan matrix
%CG=findCG(Ytensor1,Ytensor2,Ydsum1,Ydsum2);
K = kron(eye(D),Ydsum1) - kron(Ytensor1',eye(D));
v = null(K);

% The subspace is at least one, but if bigger
% dimensional say so.

[hgt,wth]=size(v);
if (wth > 1)
    fprintf(1,'Null space is %d-dimensional\n',wth);
end;

R = reshape(v(:,1),d,d);  % this won't work for decomps with 
D = diag(R*R');           % multiplicities > 1
D2 = diag(1./sqrt(D));
C = D2*R;