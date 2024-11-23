function T= genn11tableaux (n)
% usage
%             T = genn11tableaux(n)
% generates the Young tableaux of the (n-2,1,1) partition for n >= 5.
% any such tableau is determined by its entries in the second and third
% rows
% The tableaux has the shape
%  1 a... x
%  y 
%  z
%  Now, z > 2 
%  the number of possible tableaux is by the hook length formula
%  prod
%      n n-3 n-4 ... 1
%      2
%      1
%    Nt = n! / hook lengths = (n-1)*(n-2)/2;

if (n < 5)
  disp('only works for n >= 5');
  return;
end;
Nt = (n-1)*(n-2)/2;
T  = zeros(Nt,n);
T(:,1) = 1;  % upper left entry always 1
k = 1;
for y = n-1:-1:2
    for z = y+1:n
        row = 2:n;
        row = row(row ~= y);
        row = row(row ~= z);
        T(k,2:n-2) = row;
        T(k,n-1)= y;
        T(k,n) = z;
        k = k + 1;
    end;
end;

