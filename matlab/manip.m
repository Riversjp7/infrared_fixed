%%

% taken from ../tc/running_coupling/nf12_c0.2_v2 
%for Nf = 12, adapt for Nf = 10

% s = 2: 12 -> 24, 16 -> 32, 18 -> 36, 20 -> 40, 24 -> 48

% same as ...tuning.m now for c = 0.30

% v3: more beta values, for same set of volumes ; change interpolation to
% O(beta^4)

% v4: different data cuts, remove 30% of each run; same beta ranges and
% O(beta^4)

% v5: add 8 -> 16,10 -> 20, with 30% cuts

clear all

format long

N = 61;
M0=5;% M is updated immediately in the for loop for Q vs M
M=M0;% use M0= const value for single fit, use M0= lowest value for Q vs M.
QvsM=0 % set to 0 or 1. 0 is original single fit. 1 is range of fits 
        % from M0+1 to Mmax

if QvsM==0
    Mmax=M0;
elseif QvsM==1
    Mmax=8;
end
chi24perDof=zeros(1,Mmax);
Q24=zeros(1,Mmax);
chi12perDof=zeros(1,Mmax);
Q12=zeros(1,Mmax);
for m=1:Mmax

if QvsM==1
    M=M+1;
end 
% 12 -> 24

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7

format long

[xfit,gfit,gfit_err]=get_info(12);

figure(121)
box on
grid minor 
hold on

x = xfit;
y = gfit;
err = gfit_err;

errorbar(x,y,err,'ro','LineWidth',1);

weight=1./gfit_err;

B=gfit./gfit_err;

Ndat=length(xfit);

TEMP=tempPop(xfit,M,weight);

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p12=C*(A'*B');

p12_error=sqrt(diag(C));

p12'
p12_error'

chi12=(B-(A*p12)')*(B'-A*p12);
DoF = Ndat - M; 
chi12_dof=chi12/DoF
Pvalue=gammainc(chi12/2,DoF/2)
Qvalue=1-Pvalue
chi12perDof(m)=chi12;
Q12(m)=Qvalue;

xp=2.5:0.01:4.2;
yp= Vpn(p12,xp);
plot(xp,yp,'-b','LineWidth',1);



xroot12(1:N)=0;
error12Tuned(1:N)=0;
g12Tuned(1:N)=0;
beta12Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p12,x) - g2Target,xInit);

xroot12(i)=xr;
beta12Tuned(i)=xroot12(i);

xpp=xroot12(i);
ypp=Vpn(p12,xpp);
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error12Tuned(i)=delta;
g12Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot12,g12Tuned,error12Tuned,'bo','LineWidth',1);



[xfit,gfit,gfit_err]=get_info(24)

% figure(22)
% box on
% grid minor 
% hold on

x = xfit;
y = gfit;
err = gfit_err;

errorbar(x,y,err,'go','LineWidth',1);

weight=1./gfit_err;

B=gfit./gfit_err;

clear Ndat X1 X2 X3 X4 X5

Ndat=length(xfit);

TEMP=tempPop(xfit,M,weight);

clear A B C H U S V w winv WW

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p24=C*(A'*B');

p24_error=sqrt(diag(C));

p24'
p24_error'

chi24=(B-(A*p24)')*(B'-A*p24);
DoF = Ndat - 5; 
chi24_dof=chi24/DoF
Pvalue=gammainc(chi24/2,DoF/2)
Qvalue=1-Pvalue
chi24perDof(m)=chi24;
Q24(m)=Qvalue;

xp=2.5:0.01:4.2;
yp2=Vpn(p24,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot12(i);
ypp=Vpn(p24,xpp);
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
%xvec=[1; xpp; xpp^2; xpp^3; xpp^4];
% delta=ypp^2*sqrt(xvec'*(C*xvec));
delta=sqrt(xvec*(C*xvec'));
ypp;

error24(i)=delta;
g24stepped(i)=ypp;

end

% g24stepped'
% error24'

%errorbar(xroot12,g24stepped,error24,'ko','LineWidth',1);

text(3,5,['12^4 \rightarrow 24^4'],'fontsize',16)
text(3,4.5+(M-2).*0.5,['\chi^2/dof= ',num2str(chi12_dof,2),' and ',num2str(chi24_dof,2)],'fontsize',14)

xlim([2.5 4.2])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,7.4,'12^4','fontsize',16)
text(3.6,6.5,'24^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)
hold off

end
if QvsM==1
mrange=M0+1:M
figure(1)
box on
grid minor 
plot(mrange,Q24,'^r')
title('Q24/M')
figure(2)
box on
grid minor 
plot(mrange,Q12,'^b')
title('Q12/M')
figure(3)
box on
grid minor 
plot(mrange,chi24perDof,'*r')
title('chi^2/M (24)')
figure(4)
box on
grid minor 
plot(mrange,chi12perDof,'*b')
title('chi^2/M (12)')
end
