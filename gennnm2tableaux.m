function T= gennnm2tableaux (n)
% usage
%             T = gennmm2tableaux(n)
% generates the Young tableaux of the (n-2,2) partition for n >= 5.
% any such tableau is determined by its entries in the second row.
% The tableaux has the shape
%  1 a... x
%  y z
%  Now, z > 3 because if z = 3, then y = 2, and a = 4.  
% the number of possible tableax is 
% l1 = n - 2 + 1 = n-1, l2 = 2;
% Frobenius 
%    Nt = n! (l1 - l2 ) / ((l1!) (l2!)) = n * (n-3)/2

if (n < 5)
  disp('only works for n >= 5');
  return;
end;
Nt = n*(n-3)/2;
T  = zeros(Nt,n);
T(:,1) = 1;
k = 1;
z = n;
while z >= 4
  for y = z-1:-1:2
    row = 2:n;
    row = row(row ~= y);
    row = row(row ~= z);
   % disp(size(row));
    %disp(row);
    T(k,2:n-2) = row;
    T(k,n-1)=y;
    T(k,n) = z;
   % disp(T(k,1:n-2));
   % disp(T(k,n-1:n));
    k = k + 1;
  end;
  z = z - 1;
end;

