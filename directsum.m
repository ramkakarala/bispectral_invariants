function Z = directsum(X,Y)
% usage
%        Z = directsum(X,Y)
% utility that output out a matrix direct sum of X and Y
% r kakarala
Z = zeros(size(X)+size(Y));
[M,N] = size(X);
Z(1:M,1:N) = X;
Z(M+1:end,N+1:end) = Y;