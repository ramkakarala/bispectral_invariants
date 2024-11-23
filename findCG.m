function C = findCG(X12,X12n,Y12,Y12n)
% usage
%
%        C = findCG(X12,X12n,Y12,Y12n)
% where X12 is the value of a rep at (12)
% and X12n is the value of the same rep at (12..n)
% Y12 is the direct sum decomp and Y12n at (12) and (12..n)
% finds the value of the CG matrix using Jonathan Huang's algorithm from
% the paper titled
% "Fourier Theoretic Probabilistic Inference over Permutations" 
%
% r kakarala
% ntu

d = max(size(X12));

K1 = kron(eye(d),Y12) - kron(X12',eye(d));
K2 = kron(eye(d),Y12n) - kron(X12n',eye(d));
K = [K1;K2];
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
