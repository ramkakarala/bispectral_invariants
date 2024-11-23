function R = irrepnn1(p)
% R = irrepnn1(p)
% Returns the irreducible representation matrix in partition (n-1,1)
% for permutation p.
% r kakarala
% ntu
% 3 nov 08

n = length(p);
R = eye(n-1);
% we can decompose p into a product of adjacent swaps
% by using the bubblesort algorithm. For example, the permutation
% 1->4, 2->1, 3->2, 4->3 is in one-line notation [4 1 2 3]
% the cycle (4321). The sorting of one-line notation gives
% [4 1 2 3] -> [1 4 2 3] (swap (12)) -> [1 2 4 3] (swap (23)) -> [1 2 3 4]
% swap (34).  This means 
%          (4321) = (34)(23)(12)
% where operations are applied from right to left.  
% this idea from J. Huang's CMU tech report CMU-ML-08-108
for i = 1:n-1
    done = 1;
    for k = 1:n-1

        if (p(k) > p(k+1)) 
            t = p(k);          % swap
            p(k) = p(k+1);
            p(k+1) = t;
            %disp([k k+1 p]);
            done = 0;
            R = irrepnn1adjacent(n,k+1) * R;  % order goes right to left
        end;

    end;

    if (done)
        break;
    end;
end;

   
%-----------------------------------------------
function R = irrepnn1adjacent(n,i)
% usage
%        R = irrepnn1adjacent(n,i)
% computes the matrix of irreducible rep for partition (n-1,1)
% and adjacent transposition (i-1,i), i=2...n

R = zeros(n-1);
for k = 2:n   % step through tableau
    d = axialdistance(k,i-1,i);
    R(k-1,k-1) = 1/d;
    % now to find the tableau whose permutation by (i-1,i)
    % yields k.  Note that the k-th tableau is
    %   1 2 ... k-1 k+1 .. n
    %   k
    % if permuted by (i-1,i) it will not be a tableau
    % unless either i = k or i-1 = k and i > 2
    if (i > 2)
           if (i ~= k)
               R(k-1,i-1) = sqrt(1-1/d^2);
               R(i-1,k-1) = R(k-1,i-1);
           else
               R(k-1,i-2) = sqrt(1-1/d^2);
               R(i-2,k-1) = R(k-1,i-2);
           end;
    end;
end;

R = flipud(fliplr(R));  % ok this is a hack, but matches the printed tables
%------------------------------------------
function z = axialdistance(k,x,y)
% usage
%        z = axialdistance(k,x,y)
% computes the axial distance from x to y in the 
% standard tableau k for partition (n-1,1)

z = col(k,y) - col(k,x) - (row(k,y) - row(k,x));

%--------------------------------------------------------
function y = row(k,x)
% usage
%        y = row(k,x)
% finds the row number of x of standard tableau k in the partition (n-1,1)
% Any such tableau has the form
%   1 2 ... k-1 k+1 .. n
%   k

if (x == k)
    y = 2;
else
    y = 1;
end;
%-------------------------------------------------------
function y = col(k,x)
% usage
%        y = col(k,x)
% finds the column number of x of standard tableau k in the partition (n-1,1)
% Any such tableau has the form
%   1 2 ... k-1 k+1 .. n
%   k

if (x==k || x==1)
    y = 1; 
else
    if (x < k)
        y = x;
    else
        y = x-1;
    end;
end; 