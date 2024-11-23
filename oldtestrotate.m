
%% experiment on rotation invariants using bispectral forms
clear;
N = 4;  % max frequency
N2 = ceil(N/2);
N3 = N-N2;
load('M2Y08.mat');% data is in cell array A
load('CG08.mat'); % data is in CG
load('stereodata/150.mat');
FACE1=d_vector;
FACE1=d_vector(find(d_vector(:,3)),:);
for mapind = 0:300
    fname=sprintf('stereodata/%1d.mat',mapind);
    disp(fname);
    load(fname);
    FACE2=d_vector(find(d_vector(:,3)),:);
    %% randomly rotate about mean
    % cent = mean(FACE1);
    % centeredFACE1 = FACE1 - ones(max(size(FACE1)),1)*cent;
    % [U,S,V]=svd(randn(3));
    % R = U*V';
    % R = det(R)*R;  % make sure it is rotation (det +1), not reflection
    % centeredFACE2 = centeredFACE1*R;  % rotate
    % FACE2 = centeredFACE2 + ones(max(size(FACE1)),1)*randn(1,3);
    %load('stereodata/250.mat');
    %FACE2=d_vector;
    %% scale
    % maxval = max(max(FACE1(:)),max(FACE2(:)));
    % FACE1 = FACE1/maxval;
    % FACE2 = FACE2/maxval;
    %% compute centroid
    xc = 0; yc = 0; zc = 0;
    xcr = 0; ycr = 0; zcr = 0;
    Npts = max(size(FACE1));
    np = 0;
    for p = 1:Npts
        
        np = np+1;
        xc = xc + FACE1(p,1);
        yc = yc + FACE1(p,2);
        zc = zc + FACE1(p,3);
        
        
    end;
    cent = [xc yc zc]/np;
    Npts2 = max(size(FACE2));
    np  = 0;
    for p = 1:Npts2
        
        np = np+1;
        xcr = xcr + FACE2(p,1);
        ycr = ycr + FACE2(p,2);
        zcr = zcr + FACE2(p,3);
        
        
    end;
    centr = [xcr ycr zcr]/np;
    
    %% compute central moments of origignal data
    for n = 1:N
        
        fprintf(1,'Evaluating moments of order %d\n',n);
        nmono = (n+1)*(n+2)/2;
        mvecx = zeros(nmono,1);
        
        lmnc = 0;
        for l = n:-1:0
            for m = 0:n-l
                
                lmnc = lmnc + 1;
                p = n-l-m;
                momx = 0; momy = 0;
                %   disp([n,l,m,p]);
                for k = 1:Npts
                    xvec = FACE1(k,1:3)-cent;
                    momx = momx + ((xvec(1)^l) * (xvec(2)^m) * (xvec(3)^p));
                end;
                
                mvecx(lmnc) = momx;
                
            end;
        end;
        
        Mcvecx{n} = mvecx;
       % disp(mvecx);
    end;
    
    
    %% compute central moments of Br rotated data
    for n = 1:N
        
        fprintf(1,'Evaluating moments for transformed data of order %d\n',n);
        nmono = (n+1)*(n+2)/2;
        mvecy = zeros(nmono,1);
        
        lmnc = 0;
        for l = n:-1:0
            for m = 0:n-l
                
                lmnc = lmnc + 1;
                p = n-l-m;
                momy = 0;
                %   disp([n,l,m,p]);
                for k = 1:Npts2
                    xrvec = FACE2(k,1:3)-centr;
                    momy = momy + ((xrvec(1)^l) * (xrvec(2)^m) * (xrvec(3)^p));
                    
                end;
                
                mvecy(lmnc) = momy;
                
            end;
        end;
        
        Mcvecy{n} = mvecy;
        
    end;
    
    %% Obtain Y vectors from moments
    for k = 1:N
        Yx{k} = A{k}*Mcvecx{k};
        Yy{k} = A{k}*Mcvecy{k};
    end;
    
    %% test rotation by power spectra
    
    for k = 2:N
        nx = norm(Yx{k});
        ny = norm(Yy{k});
        ep(k) = (nx-ny)/(nx+ny);
    end;
    if max(abs(ep)>5e-2)
        fprintf('Power spectrum test fails\n');
        disp(ep);
    else
        fprintf(1,'Power spectrum test passed \n');
    end;
    %clear e;
    
    %% test bispectral invariants
    
    xt = kron(Yx{N2}',Yx{N3}');
    yt = kron(Yy{N2}',Yy{N3}');
    
    if (N2==N3)
        xs = [1];
        ys = [1];
        startind = 1;
        indices = [1];
    else
        xs = []; ys = [];
        startind = 0;
        indices = [];
    end;
    
    for p = max((N2-N3),1):N
        dim = 2*p+1;
        Dx = zeros(dim);
        Dx(:,p+1) = Yx{p};
        Dy = zeros(dim);
        Dy(:,p+1) = Yy{p};
        xs = directsum(xs,Dx);
        ys = directsum(ys,Dy);
        indices = [indices startind+p+1];
        startind = startind + dim;
    end;
    
    for k = 1:max(size(xs))
        colnorms(k) = norm(xs(:,k))+norm(ys(:,k));
    end;
    if (indices ~= find(colnorms))
        fprintf(1,'Error! could not find nonzero indices\n');
        return;
    end;
    C = CG{N2+1,N3+1};
    LHS = xt*C*xs; %*C';
    LHSreduced = LHS(indices);
    RHS = yt*C*ys; %*C';
    RHSreduced = RHS(indices);
    eb = norm(LHSreduced-RHSreduced)/(norm(LHSreduced)+norm(RHSreduced));
    if (eb > 5e-2)
        fprintf(1,'Error! bispectral invariants do not work: error=%f\n',eb);
    else
        fprintf(1,'Bispectral invariants check out OK\n');
    end;
    ebrecord(mapind+1)=eb;
    eprecord(mapind+1,:)=ep;
end;
%% plot results
plot(0:300,ebrecord,'b',0:300,max(eprecord,[],2),'r');
save stereodata/momentbasedresults8feb;
return;
%
%% Create Y vectors directly from Data
for k=1:N
    Fx{k} = zeros(2*k+1);
    
end;
t0 = clock;
h = waitbar(0,'Computing Y vectors of data...');
for k=1:100:Npts
    waitbar(k/Npts,h);
    v = FACE1(k,1:3)-cent;
    r = norm(v);
    alpha = atan2(v(2),v(1));  % azimuth
    beta = acos(v(3)/r);    % elevation
    
    for l = 1:N
        ind = -l:l;
        Dm = zeros(2*l+1);
        Dm(:,l+1) = dmatrixmiddle(l,alpha,beta);
        %  test = exp(-j*ind'*alpha)*[zeros(1,l),1,zeros(1,l)].*DC{l+1};
        % only the middle column of D will be nonzero. Yes, I know it
        % is strange to compute all of D, but I want to try different
        % options later with the full matrix
        Fx{l} = Fx{l} + r*Dm';
        % use D', not D because of forward Fourier Trans
        
    end;
end;
close(h);
fprintf(1,'Took %6.2f min \n',etime(clock,t0)/60);
%% do the same for rotated data
for k=1:N
    Fy{k} = zeros(2*k+1);
    
end;
t0 = clock;
hr = waitbar(0,'Computing Y vectors of rotated data...');
for k=1:100:Npts2
    waitbar(k/Npts,hr);
    v = FACE2(k,1:3)-centr;
    ry = norm(v);
    alphay = atan2(v(2),v(1));  % azimuth
    betay = acos(v(3)/ry);    % elevation
    %  DC = dmatrixbetark(N,betay(k));
    for l = 1:N
        ind = -l:l;
        Dm = zeros(2*l+1);
        Dm(:,l+1) = dmatrixmiddle(l,alphay,betay);
        % D = exp(-j*ind'*alphay(k))*[zeros(1,l),1,zeros(1,l)].*DC{l+1};
        % only the middle column of D will be nonzero. Yes, I know it
        % is strange to compute all of D, but I want to try different
        % options later with the full matrix
        Fy{l} = Fy{l} + ry*Dm';
        % use D', not D because of forward Fourier Trans
        
    end;
end;
close(hr);
fprintf(1,'Took %6.2f \n',etime(clock,t0));
%% test rotation by power spectra

for k = 1:N
    nx = norm(Fx{k});
    ny = norm(Fy{k});
    epd(k) = (nx-ny)/max(1,(nx+ny));  % max with one avoids really small nums
    if (abs(epd(k)) > 1e-2)
        fprintf(1,'Error! x and y cannot be rotations; check %d!\n',k);
        return;
    end;
end;
fprintf(1,'Power spectrum test passed \n');

%% now compute bispectral invariants
xt = kron(Fx{N2},Fx{N3});
yt = kron(Fy{N2},Fy{N3});
if (N2==N3)
    xs = [0];
    ys = [0];
else
    xs = []; ys = [];
end;

for p = 1:N
    dim = 2*p+1;
    Dx = Fx{p}';
    Dy = Fy{p}';
    xs = directsum(xs,Dx);
    ys = directsum(ys,Dy);
end;
C = CG{N2+1,N3+1};
LHS = xt*C*xs*C';
RHS = yt*C*ys*C';
scale = norm(LHS,'fro')+norm(RHS,'fro');
ebd = norm(LHS-RHS,'fro')/scale;
if (ebd > 1e-2)
    fprintf(1,'Error! bispectral invariants do not work\n');
else
    fprintf(1,'Bispectral invariants check out OK\n');
end;