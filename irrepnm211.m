function R = irrepnm211(p,T)
% R = irrepnm211(p)
% Returns the irreducible representation matrix in partition (n-2,1,1)
% for permutation p.
% r kakarala
% ntu
% 3 nov 08

n = length(p);
[Nt,nc] = size(T);
R = eye(Nt);
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
            R = irrepnm11adjacent(k,n,Nt,T) * R;  % order goes right to left
        end;

    end;

    if (done)
        break;
    end;
end;

function D= irrepnm11adjacent(k,n,Nt,T)
% usage
%             D = irrepnm11adjacent(k,n,T)
% r kakarala

if (k >= n)
	disp('error! k must be < n');
	return;
end;
D = zeros(Nt);
for r = 1:Nt
      d = axialdist(T(r,:),k,k+1);
      D(r,r) = 1/d;
      if (d ~= 1)
          Ts = T(r,:);
          for c = 1:n
              if T(r,c)==k
                  Ts(c)=k+1;
              elseif T(r,c) == k+1
                  Ts(c) = k;
              end;
          end;
	for c = r+1:Nt
	   if (T(c,:)==Ts)
	       dp = sqrt(1 - 1/d^2 );
	       D(r,c) = dp;
	       D(c,r) = dp;
	       break;
	   end;
	end;
      end;
end;

function d = axialdist(tab,a,b)
% finds the axial distance between a, b in tableaux (n-2,1,1)
% r kakarala
n = length(tab);
ca = find(tab == a);
if (ca <= n-2)
	ra = 1;
else
    if (ca == n-1)
        ra = 2;
    else 
        ra = 3;
    end;
end;
if (ra > 1)
    ca = 1;   % in the first column
end;
cb = find(tab == b);
if (cb <= n-2)
	rb = 1;
else
    if (cb == n-1)
        rb = 2;
    else 
        rb = 3;
    end;
end;
if (rb > 1)
    cb = 1;   % in the first column
end;
d = (cb-ca) - (rb - ra);

