%% collects experiments on S5 from SPL paper.
clear;
close all;
%% experiment for bispectrum on S5
S5 = perms(1:5);
Np = max(size(S5));
%% generate reps d41 and d2111, d32 and d221
d41 = zeros(4,4,Np);
d2111 = d41;
d32 = zeros(5,5,Np);
d311 = zeros(6,6,Np);
d221 = d32;
Tab32 = gennnm2tableaux(5);
Tab311 = gennn11tableaux(5);
for k = 1 : Np
    
    d41(:,:,k) = irrepn11(S5(k,:));
    d2111(:,:,k) = irrep21dd1(S5(k,:));
    d1(k) = irrep1(S5(k,:));
    d32(:,:,k) = irrepnm22(S5(k,:),Tab32);
    d221(:,:,k) = irrep2nm2(S5(k,:),Tab32);  % 2,3 really means (2,2,1) here
    d311(:,:,k) = irrepnm211(S5(k,:),Tab311);
 
end;
%% compute characters
for k = 1 : Np;
    chi41(k) = trace(d41(:,:,k));
    chi2111(k) = trace(d2111(:,:,k));
    chi32(k) = trace(d32(:,:,k));
    chi221(k) = trace(d221(:,:,k));
    chi311(k) = trace(d311(:,:,k));
    chi5(k) = 1;
    chi11111(k) = d1(k);
end;
%% decompose character of kronecker products
KD=zeros(7,7,7);  % seven reps: (5),(4,1),(3,2),(3,1,1),(2,2,1),(2,1,1,1),(1,1,1,1,1);
Chi=zeros(Np,7);
Chi(:,1)=chi5(:); % save characters in a convenient form
Chi(:,2)=chi41(:);
Chi(:,3)=chi32(:);
Chi(:,4)=chi311(:);
Chi(:,5)=chi221(:);
Chi(:,6)=chi2111(:);
Chi(:,7)=chi11111(:);
% test orthogonality
fprintf(1,'orthogonality\n');
for i=1:7
    for j=1:7
        fprintf(1,'%1d ',round(sum(Chi(:,i).*Chi(:,j))/120));
    end;
    fprintf(1,'\n');
end;
% note: orthonormality (when /120) suffices for irreducibility
fprintf(1,'\n');
% print out k-product decompositions
repstring = char(zeros(7,13));
repstring(1,:) = '      (5)    ';
repstring(2,:) = '     (4,1)   ';
repstring(3,:) = '     (3,2)   ';  
repstring(4,:) = '     (3,1^2) ';
repstring(5,:) = '     (2^2,1) ';
repstring(6,:) = '     (2,1^3) ';  
repstring(7,:) = '     (1^4)   ';
fprintf('Kronecker product decomposition\n');
fprintf(1,'            ');
for i = 1:7
    fprintf(1,' %s ',repstring(i,:));
end;
fprintf(1,'\n');
for i = 1:7
    fprintf(1,'$%s$ ',repstring(i,:));
    for j = 1:7
        kr = Chi(:,i).*Chi(:,j);  % character of k-prod 
		fprintf(1,'& $');
        for p = 1:7
            KD(i,j,p) = sum(kr.*Chi(:,p))/120; % inner product
            fprintf(1,'%1d',round(KD(i,j,p)));
        end;
        
        fprintf(1,'$    ');
    end;
    fprintf(1,'\\ \n');
end;
%% Clebsh-Gordan matrix for (4,1), (1,1,1,1,1)
X12 = kron(irrepn11([2,1,3,4,5]),irrep1([2,1,3,4,5]));
X12n = kron(irrepn11([2 3 4 5 1]),irrep1([2 3 4 5 1]));
Y12 = irrep21dd1([2,1,3,4,5]);
Y12n = irrep21dd1([2,3,4,5,1]);
C411 = findCG(X12,X12n,Y12,Y12n);
%% test the CG matrix for (4,1), (1,1,1,1,1)
for k=1:Np
    LHS = kron(d41(:,:,k),d1(k));
    RHS = C411' * d2111(:,:,k) * C411;
    e = norm(LHS-RHS);
    if (e > 1e-6)  % if we see error, then CG matrix is not right
       disp([k,norm(LHS-RHS)]);
    end;
end;
%% Clebsh-Gordan matrix for (3,2), (1,1,1,1,1)
X12 = kron(irrepnm22([2,1,3,4,5],Tab32),irrep1([2,1,3,4,5]));
X12n = kron(irrepnm22([2 3 4 5 1],Tab32),irrep1([2 3 4 5 1]));
Y12 = irrep2nm2([2,1,3,4,5],Tab32);
Y12n = irrep2nm2([2,3,4,5,1],Tab32);
C321 = findCG(X12,X12n,Y12,Y12n);
%% test the CG matrix for (3,2), (1,1,1,1,1)
for k=1:Np
    LHS = kron(d32(:,:,k),d1(k));
    RHS = C321' * d221(:,:,k) * C321;
    e = norm(LHS-RHS);
    if (e > 1e-6) % if we see error, then CG matrix is not right % error
       disp([k,norm(LHS-RHS)]);
    end;  
end;
%% Clebsh-Gordan matrix for (4,1) (3,2)
for k=1:120 
    if (S5(k,:)==[2 1 3 4 5]) 
        k12=k; 
    end; 
end;
X12 = kron(d41(:,:,k12),d32(:,:,k12)); 
for k=1:120 
    if (S5(k,:)==[2 3 4 5 1]) 
        kcycle=k; 
    end; 
end;
X12n = kron(d41(:,:,kcycle),d32(:,:,kcycle)); 
Y12 = directsum(directsum(directsum(d41(:,:,k12),d32(:,:,k12)),d311(:,:,k12)),d221(:,:,k12));
Y12n = directsum(directsum(directsum(d41(:,:,kcycle),d32(:,:,kcycle)),d311(:,:,kcycle)),d221(:,:,kcycle));
C4132 = findCG(X12,X12n,Y12,Y12n);
%% test that CG matrix for (4,1) (3,2)
for k=1:Np
    LHS = kron(d41(:,:,k),d32(:,:,k));
    MRHS =  directsum(directsum(directsum(d41(:,:,k),d32(:,:,k)),d311(:,:,k)),d221(:,:,k));
    RHS = C4132' * MRHS * C4132;
    e = norm(LHS-RHS);
    if (e > 1e-6)  % error
       disp([k,norm(LHS-RHS)]);
    end;  
end;
%% Clebsh-Gordan matrix for (4,1) (3,1,1)
for k=1:120 
    if (S5(k,:)==[2 1 3 4 5]) 
        k12=k; 
    end; 
end;
X12 = kron(d41(:,:,k12),d311(:,:,k12)); 
for k=1:120 
    if (S5(k,:)==[2 3 4 5 1]) 
        kcycle=k; 
    end; 
end;
X12n = kron(d41(:,:,kcycle),d311(:,:,kcycle)); 
Y12 = directsum(directsum(directsum(directsum(d41(:,:,k12),d32(:,:,k12)),d311(:,:,k12)),d221(:,:,k12)),d2111(:,:,k12));
Y12n = directsum(directsum(directsum(directsum(d41(:,:,kcycle),d32(:,:,kcycle)),d311(:,:,kcycle)),d221(:,:,kcycle)),d2111(:,:,kcycle));
C41311 = findCG(X12,X12n,Y12,Y12n);
%% test that CG matrix for (4,1) (3,1,1)
for k=1:Np
    LHS = kron(d41(:,:,k),d311(:,:,k));
    MRHS =  directsum(directsum(directsum(directsum(d41(:,:,k),d32(:,:,k)),d311(:,:,k)),d221(:,:,k)),d2111(:,:,k));
    RHS = C41311' * MRHS * C41311;
    e = norm(LHS-RHS);
    if (e > 1e-6)  % error
       disp([k,norm(LHS-RHS)]);
    end;  
end;
%% Clebsh-Gordan matrix for (3,2) (3,1,1) -- don't use, multiplicities
cause error
for k=1:120 
    if (S5(k,:)==[2 1 3 4 5]) 
        k12=k; 
    end; 
end;
X12 = kron(d32(:,:,k12),d311(:,:,k12)); 
for k=1:120 
    if (S5(k,:)==[2 3 4 5 1]) 
        kcycle=k; 
    end; 
end;
X12n = kron(d32(:,:,kcycle),d311(:,:,kcycle)); 
ds3112 = directsum(d311(:,:,k12),d311(:,:,k12));
Y12 = directsum(directsum(directsum(directsum(d41(:,:,k12),d32(:,:,k12)),ds3112),d221(:,:,k12)),d2111(:,:,k12));
ds3112 = directsum(d311(:,:,kcycle),d311(:,:,kcycle));
Y12n = directsum(directsum(directsum(directsum(d41(:,:,kcycle),d32(:,:,kcycle)),ds3112),d221(:,:,kcycle)),d2111(:,:,kcycle));
C32311 = findCG(X12,X12n,Y12,Y12n);
%% test that CG matrix for (3,2) (3,1,1) -- don't use
% for k=1:Np
%     LHS = kron(d32(:,:,k),d311(:,:,k));
%     ds3112 = directsum(d311(:,:,k),d311(:,:,k));
%     MRHS =  directsum(directsum(directsum(directsum(d41(:,:,k),d32(:,:,k)),ds3112),d221(:,:,k)),d2111(:,:,k));
%     RHS = C32311' * MRHS * C32311;
%     e = norm(LHS-RHS);
%     if (e > 1e-6)  % error
%        disp([k,norm(LHS-RHS)]);
%     end;  
% end;
%% load data
votes = apa1980;
%% compute spectrum
Fn = sum(votes);
F41 = zeros(4);
F2111 = zeros(4);
F32 = zeros(5);
F221 = zeros(5);
F311 = zeros(6);
F1 = 0;
for k = 1 : Np
    F41 = F41 + votes(k)*d41(:,:,k)';
    F2111 = F2111 + votes(k)*d2111(:,:,k)';
    F32 = F32 + votes(k)*d32(:,:,k)';
    F221 = F221 + votes(k)*d221(:,:,k)';
    F1 = F1 + votes(k)*d1(k)';
    F311 = F311 + votes(k)*d311(:,:,k)';
end;
%% normalize by group size to make the numbers smaller and easier to grasp
% this won't affect the rations below
Fnn = Fn/Np;
F41n = F41/Np; 
F2111n = F2111/Np; 
F1n=F1/Np;
F32n = F32/Np; 
F311n = F311/Np; 
F221n = F221/Np; 
%% print out the spectral magnitudes
% display the results
fprintf(1,' (5) (4,1)  (3,2)  (3,1^2) (2^2,1)  (2,1^3)  (1^5)\n');
fprintf(1,'%5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n',...
norm(Fnn),norm(F41n),norm(F32n),norm(F311n),norm(F221n),norm(F2111n),norm(F1n));
%% compute bispectrum
% (4,1), (3,2) value
kprod = kron(F41n,F32n);
dsum = directsum(directsum(directsum(F41n,F32n),F311n),F221n);
Bispec4132 = kprod * C4132' * dsum' * C4132;
Gtest4132 = norm(Bispec4132)^2/norm(kprod*dsum')^2;  % denominator removes effect of noise scaling
% (4,1), (3,1,1) value
kprod = kron(F41n,F311n);
dsum = directsum(directsum(directsum(directsum(F41n,F32n),F311n),F221n),F2111n);
Bispec41311 = kprod * C41311' * dsum' * C41311;
Gtest41311 = norm(Bispec41311)^2/norm(kprod*dsum')^2;  % denominator removes effect of noise scaling

% (4,1), (1,1,1,1,1) value
kprod = kron(F41n,F1n);
dsum = F2111n;
Bispec411 = kprod * C411' * dsum' * C411;
Gtest411 = norm(Bispec411)^2/norm(kprod*dsum')^2;
% (3,2), (1,1,1,1,1) value
kprod = kron(F32n,F1n);
dsum = F221n;
Bispec3223 = kprod * C321' * dsum' * C321;
Gtest3223 = norm(Bispec3223)^2/norm(kprod*dsum')^2;
% print out the result
fprintf(1,'Gamma (nonGaussianity) value for votes data = %5.3f, %5.3f, sum= %5.3f\n',Gtest4132,Gtest41311,Gtest4132+Gtest41311);
%% Compute statistics
Ntrials = 1000;
ftop = zeros(1,Np);
Gtop = 0;
Gtestrand4132 = zeros(1,Ntrials);
Gtestrand41311 = zeros(1,Ntrials);

fprintf(1,'Computing empirical pdf of Gamma statistic..\n');
for t = 1:Ntrials;
    f = randn(1,Np);
    F41 = zeros(4); 
    F2111 = zeros(4);
    F32 = zeros(5);
    F221 = zeros(5);
    F1 = 0;
    F311 = zeros(6);
    for k = 1:Np
     
      F41 = F41 + f(k)*d41(:,:,k)';
      F2111 = F2111 + f(k)*d2111(:,:,k)';
      F1 = F1 + f(k)*d1(k)';
      F32 = F32 + f(k)*d32(:,:,k)';
      F221 = F221 + f(k)*d221(:,:,k)';
      F311 = F311 + votes(k)*d311(:,:,k)';   
      
    end;
    % normalize by group length 
    F41n = F41/Np; 
    F2111n = F2111/Np; 
    F1n=F1/Np; 
    F32n = F32/Np; 
    F311n = F311/Np; 
    F221n = F221/Np; 
	
    kprod = kron(F41n,F32n);
    dsum = directsum(directsum(directsum(F41n,F32n),F311n),F221n);
    Bispec4132 = kprod*C4132'*dsum'*C4132;
    Grand4132 = norm(Bispec4132)^2/norm(kprod*dsum')^2;  
    % (4,1), (3,1,1) value
    kprod = kron(F41n,F311n);
    dsum = directsum(directsum(directsum(directsum(F41n,F32n),F311n),F221n),F2111n);
    Bispec41311 = kprod * C41311' * dsum' * C41311;
    Grand41311 = norm(Bispec41311)^2/norm(kprod*dsum')^2;  % denominator removes effect of noise scaling

    Gtestrand4132(t) = Grand4132; %
    Gtestrand41311(t) = Grand41311;

  
    if (rem(t,100)==0)
        disp(Ntrials-t);
    end;
end;
%% print out the p-value of the votes date
pval = sum(Gtestrand4132 >= Gtest4132)/length(Gtestrand4132);
fprintf(1,'Empirical pvalue of votes data is %5.3f\n',pval); 
%% fit lognormal
lognparms4132 = lognfit(Gtestrand4132);
[hh4132,hi4132] = hist(Gtestrand4132,0.1:0.2:max(Gtestrand4132));
Npdf4132 = lognpdf(hi4132,lognparms4132(1),lognparms4132(2));
epdf4132 = hh4132/sum(hh4132);
plot(hi4132,epdf4132,'--',hi4132,Npdf4132/sum(Npdf4132),'-');
legend('Experimental','Lognormal fit');
%% significance w.r.t. the lognormal fit
pvalue = logncdf(Gtest4132,lognparms4132(1),lognparms4132(2));
fprintf(1,'(4,1), (3,2) value is significant with p = %f\n',1-pvalue); 
%% load data again to test Fourier Transform
votes = apa1980;
%% recompute spectrum
Fno = sum(votes);
F41o = zeros(4);
F2111o = zeros(4);
F32o = zeros(5);
F221o = zeros(5);
F311o = zeros(6);
F1o = 0;
for k = 1 : Np
    F41o = F41o + votes(k)*d41(:,:,k)';
    F2111o = F2111o + votes(k)*d2111(:,:,k)';
    F32o = F32o + votes(k)*d32(:,:,k)';
    F221o = F221o + votes(k)*d221(:,:,k)';
    F1o = F1o + votes(k)*d1(k)';
    F311o = F311o + votes(k)*d311(:,:,k)';
end;
%% can we reconstruct the function?
f = zeros(Np,1);
for k = 1:Np
    
    f(k) = Fno/Np;
    f(k) = f(k) + max(size(F41o))/Np * trace(d41(:,:,k)*F41o);
    f(k) = f(k) + max(size(F32o))/Np * trace(d32(:,:,k)*F32o);
    f(k) = f(k) + max(size(F311o))/Np * trace(d311(:,:,k)*F311o);
    f(k) = f(k) + max(size(F221o))/Np * trace(d221(:,:,k)*F221o);
    f(k) = f(k) + max(size(F2111o))/Np * trace(d2111(:,:,k)*F2111o);
    f(k) = f(k) + 1/Np * trace(d1(k)*F1o);
 
end;
if (norm(votes - f) > 1e-6)
    disp('Cannot reconstruct original vote data\n');
else
    disp('Fourier transform and inverse are working fine\n');
end;
%% Generate filtered Gaussian by using the inverse Fourier transform
% this was not in the paper, a test of whether correlated Gaussian
% data works
Ntrials = 1000;
ftop = zeros(1,Np);
Gtop = 0;
fprintf(1,'Statistics for correlated Gaussian process..\n');
mux = zeros(Np,1);
Rxx = zeros(Np,Np);

% find the square roots of the election data coefficient matrices
[U,S,V] = svd(F41o*F41o');
H41 = U*sqrt(S)*V';
[U,S,V] = svd(F32o*F32o');
H32 = U*sqrt(S)*V';
[U,S,V] = svd(F311o*F311o');
H311 = U*sqrt(S)*V';
[U,S,V] = svd(F221o*F221o');
H221 = U*sqrt(S)*V';
[U,S,V] = svd(F2111o*F2111o');
H2111 = U*sqrt(S)*V';
Hn = Fno;
H1 = F1o;

for t = 1:Ntrials;
   
    f = zeros(1,Np);

    Fnx = Hn*randn;
    for k = 1:Np
    
      f(k) = Fnx/Np; 
      F41x = H41*randn(size(F41));
      f(k) = f(k) + max(size(F41))/Np * trace(d41(:,:,k)*F41x);
      F32x = H32*randn(size(F32));
      f(k) = f(k) + max(size(F32))/Np * trace(d32(:,:,k)*F32x);
      F311x = H311*randn(size(F311));
      f(k) = f(k) + max(size(F311))/Np * trace(d311(:,:,k)*F311x);
      F221x = H221*randn(size(F221));
      f(k) = f(k) + max(size(F221))/Np * trace(d221(:,:,k)*F221x);
      F2111x = H2111*randn(size(F2111));
      f(k) = f(k) + max(size(F2111))/Np * trace(d2111(:,:,k)*F2111x);
      F1x = H1*randn;
      f(k) = f(k) + 1/Np * trace(d1(k)*F1x);
  
      
    end;
    
        % normalize by group length 
    F41n = F41x/Np; 
    F2111n = F2111x/Np; 
    F1n=F1x/Np; 
    F32n = F32x/Np; 
    F311n = F311x/Np; 
    F221n = F221x/Np; 

    kprod = kron(F41n,F32n);
    dsum = directsum(directsum(directsum(F41n,F32n),F311n),F221n);
    Bispec4132 = kprod*C4132'*dsum'*C4132;
    Grand4132 = norm(Bispec4132)^2/norm(kprod*dsum')^2;  
   
    Gtestrandgen4132(t) = Grand4132;
    gofresult(t) = chi2gof(f);
    llresult(t) = lillietest(f);
    
    if (rem(t,100)==0)
        disp(Ntrials-t);
    end;
end;
%% print out the p-value of the votes date
pval = sum(Gtestrandgen4132 >= Gtest4132)/length(Gtestrandgen4132);
fprintf(1,'Empirical pvalue of votes data is now %5.3f\n',pval); 
%% fit lognormal
lognparmsgen4132 = lognfit(Gtestrandgen4132);
[hh4132,hi4132] = hist(Gtestrandgen4132,0.1:0.2:max(Gtestrandgen4132));
Npdfgen4132 = lognpdf(hi4132,lognparmsgen4132(1),lognparmsgen4132(2));
epdfgen4132 = hh4132/sum(hh4132);
plot(hi4132,epdfgen4132,'--',hi4132,Npdfgen4132/sum(Npdf4132),'-');
legend('Experimental','Lognormal fit');
%% significance w.r.t. the lognormal fit
pvalue = logncdf(Gtest4132,lognparmsgen4132(1),lognparmsgen4132(2));
fprintf(1,'(4,1), (3,2) value is nows significant with p = %f\n',1-pvalue);
