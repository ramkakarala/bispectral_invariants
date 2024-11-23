function R = irrepnm22(p,T)
% R = irrepnn1(p)
% Returns the irreducible representation matrix in partition (n-1,1)
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
            R = irrepnm22adjacent(k,n,Nt,T) * R;  % order goes right to left
        end;

    end;

    if (done)
        break;
    end;
end;

function D= irrepnm22adjacent(k,n,Nt,T)
% usage
%             D = irrepnm22adjacent(k,n,T)
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
	    if Ts(c)==k
	       Ts(c)=k+1;
	    elseif Ts(c) == k+1
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
% finds the axial distance between a, b in tableaux (n-2,2)
% r kakarala
n = length(tab);
ca = find(tab == a);
if (ca <= n-2)
	ra = 1;
else
	ra = 2;
	ca = ca - (n-2);
end;
cb = find(tab == b);
if (cb <= n-2)
	rb = 1;
else
	rb = 2;
	cb = cb - (n-2);
end;
d = (cb-ca) - (rb - ra);

