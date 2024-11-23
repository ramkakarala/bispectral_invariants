%% make wigner d matrices
L = input('maximum order L='); 
beta = pi*rand; alpha=2*pi*rand; gamma=2*pi*rand;
%% check if the reps have the right character
DC = dmatrixbeta(L,beta);
charok = 1;
for k=0:L-1
    ind = -k:k;
    D = exp(j*ind'*alpha)*exp(j*ind*gamma) .* DC{k+1};
    t = real(trace(D));
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)  
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else 
            ch = 1;
        end;
    end;
    if (abs (t-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
end;
if (charok)
    fprintf(1,'Character test passed\n');
end;
%% check if the tensor product decompositions have the right
% character
tenscharok = 1;
L2 = floor(L/2);
for k = 0:L2-1
    for ell = k : L2-1
        Dt = kron(DC{k+1},DC{ell+1});
        Ds = [];
        for p=abs(ell-k):ell+k
            Ds = directsum(DC{p+1},Ds);
        end;
        trDt = trace(Dt);
        trDs = trace(Ds);
        if (abs(trDt - trDs) > 1e-6)
            fprintf(1,'Error: tensor product decom is not equiv: (%d,%d)', k,ell);
        end;
    end;
end;