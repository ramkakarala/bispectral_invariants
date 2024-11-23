function R = irrep1(p)
% R = irrep1(p,match_printed_option)
% Returns the irreducible representation matrix in partition (11...1)
% for permutation p, which is in "target" notation, i.e., the permutation
% is i -> p(i);  
% If you do use this code please cite the accompanying paper
% R Kakarala, "On a simplified Fourier transform for permutations"
% http://arxiv.org/abs/0903.5129

n = length(p);
R = 1;

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
            R = -1 * R;  % order goes right to left
        end;

    end;

    if (done)
        break;
    end;
end;
