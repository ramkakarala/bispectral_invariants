function rankings = croon1989;
% data in table from Croon (1989), as reported in
% diaconis and sturmfels ALGEBRAIC ALGORITHMS FOR SAMPLING
% FROM CONDITIONAL DISTRIBUTIONS, 1998
% r kakarala
% ntu

ranksfrompaper = [
    1234 137    2134 48     3124 330    4123 21
    1243 29     2143 23     3142 294    4132 30
    1324 309    2314 61     3214 117    4213 29
    1342 255    2341 55     3241 69     4231 52
    1423 52     2413 33     3412 70     4312 35
    1432 93     2431 59     3421 34     4321 27   
];    % correct 2431 to 59 as per                                             
      % Diaconis and Erikkson, pg 19
      
% make sure we match the data to the permutation
S4 = perms(1:4);
S4ind = S4*[1000,100,10,1]';
vind = reshape(ranksfrompaper(:,1:2:end),24,1);
vval = reshape(ranksfrompaper(:,2:2:end),24,1);
rankings = zeros(24,1);
for k = 1:24
    for t = 1:24
       if (vind(t) == S4ind(k))
           rankings(k) = vval(t);
       end;
    end;
end;


