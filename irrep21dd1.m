function R = irrep21dd1(p)
% R = irrep21dd1(p,match_printed_option)
% Returns the irreducible representation matrix in partition (2,1,...,1)
% for permutation p, which is in "target" notation, i.e., the permutation
% is i -> p(i);  
% If you do use this code please cite the accompanying paper
% R Kakarala, "On a simplified Fourier transform for permutations"
% http://arxiv.org/abs/0903.5129

n = length(p);
R = eye(n-1);

% we can decompose p, in "target" notation into a product of adjacent swaps
% by using the bubble sort algorithm.  For example
% 4321 -> 3421 -> 3241 -> 3214 -> 2314 -> 2134 -> 1234
% places swapped 
%         (12)    (23)     (34)    (12)    (23)    (12)
% in order, with action from right to left
%         (12)(23)(12)(34)(23)(12) = (4321)
% check the equality above
% the bubblesort idea from J. Huang, C. Guestrin, L. Guibas, 
% CMU tech report CMU-ML-08-108

for i = 1:n-1
    done = 1;
    for k = 1:n-1

        if (p(k) > p(k+1))
            t = p(k);
            p(k) = p(k+1);
            p(k+1) = t;
            done = 0;
            R = irrep21dd1adjacent(n,k) * R;  % order goes right to left
        end;

    end;

    if (done)
        break;
    end;
end;

%-----------------------------------------------
function R = irrep21dd1adjacent(n,k)
% usage
%        R = irrep21dd1adjacent(n,k)
% computes the matrix of an irreducible rep for partition (2,1,..,1)
% and adjacent transposition (k,k+1), k=1...n-1
% This is called the Yamanouchi matrix in Chen's book
% (GROUP REPRESENTATION THEORY FOR PHYSICISTS, pg 122) 
% and the Young Orthogonal Representation in other places such as Evans,
% below
% The matrix is sparse for transpositions, having
% one diagonal entry and at most one off-diagonal entry per row.


R = -eye(n-1);
% flip around the matrix so that the elements occur in
% reverse order from the ones in the n-1,1
if (k==1)
    R(1,1) = 1;
else % for k > 1
    R((k-1),(k-1)) = -1/k;
    R((k-1),k) = sqrt(k^2 - 1)/k;  
    R(k,(k-1)) = sqrt(k^2 - 1)/k;
    R(k,k) = 1/k;
end;
