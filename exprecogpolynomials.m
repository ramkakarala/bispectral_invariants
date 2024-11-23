% experiment with recognition polynomials
% r kakarala 17 mar 12
clear; close all;
%%
X=randn(2,5);
Pm=X-X(:,1)*ones(1,5);
xm=Pm(1,:);
ym=Pm(2,:);
%%
Ao = randn(2);
[U,S,V] = svd(Ao);
R = U*V';
s = randn(1)^2;
A = s*R;
t = rand(2,1);
%%
fprintf(1,'Test 2-D affine invariance \n');
nlevel = input('noise level (0 for noiseless case)=');
Y = Ao*X + t*ones(1,5)+nlevel*randn(2,5);
Pn = Y - Y(:,1)*ones(1,5);
xn = Pn(1,:);
yn = Pn(2,:);
%% 
figure,plot(xm,ym,'b-+');
hold on
plot(xn,yn,'r-x');
legend('template','affine data');
%%
% compute polynmoials--note formulas are different than in BHP paper
% These don't work
% a0 = xm(2)*(xn(3)*xm(5)*ym(2) - xm(3)*xn(5)*ym(2) - xn(2)*xm(5)*ym(3) ...
%     + xm(2)*xn(5)*ym(3) + xn(2)*xm(3)*ym(5) - xm(2)*xn(3)*ym(5));
% a1 = -xm(2)*(xn(3)*xm(4)*ym(2) - xm(3)*xn(4)*ym(2) + xn(2)*xm(4)*ym(3) ...
%     - xm(2)*xn(4)*ym(3) - xn(2)*xm(3)*ym(4) + xm(2)*xn(3)*ym(4));
% 
% b0 = -xm(2)*(xm(5)*yn(2)*ym(3) - xm(5)*ym(2)*yn(3) + xm(3)*yn(2)*ym(5) ...
%     - xm(2)*yn(3)*ym(5) - xm(3)*ym(2)*ym(5) + xm(2)*ym(3)*yn(5));
% b1 = xm(2)*(xm(4)*yn(2)*ym(3) - xm(4)*ym(2)*yn(3) - xm(3)*yn(2)*ym(4) ...
%     + xm(2)*yn(3)*ym(4) + xm(3)*ym(2)*yn(4) - xm(2)*ym(3)*yn(4));
% 
% T=[a1, a0; b1, b0];
% p = det(T);
%% another form
xm = xm(:); ym = ym(:);
xn = xn(:); yn = yn(:);

M = [xm(2:4), ym(2:4), xn(2:4), yn(2:4)];
b = [xm(5), ym(5), xn(5), yn(5)];

p2 = det(M*M');
fprintf(1,'Recognition polynomial gives %f\n',p2);
fprintf(1,'Small value means affine relation\n');
%%
fprintf(1,'Test 3-D rigid motion + scaling invariance \n');
nlevel = input('noise level (0 for noiseless case)=');
Y = A*X + t*ones(1,5)+nlevel*randn(2,5);
Pn = Y - Y(:,1)*ones(1,5);
xn = Pn(1,:);
yn = Pn(2,:);
figure,plot(xm,ym,'b-+');
hold on
plot(xn,yn,'r-x');
legend('template','rigid template');
%% rigid relation
D = zeros(3);
for r=1:3
    for c=1:3
        D(r,c) = xn(r+1)*xn(c+1) + yn(r+1)*yn(c+1) ...
            - xm(r+1)*xm(c+1) - ym(r+1)*ym(c+1);
    end;
end;
p3 = det(D);
fprintf(1,'Rigid recognition polynomial gives %f\n',p3);
fprintf(1,'Small value means rigid 3-D relation\n');