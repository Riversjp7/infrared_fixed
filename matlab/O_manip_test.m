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


% 12 -> 24

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7 X8 X9

format long

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L12_c0.30_ssc.dat.new3');
dat12 = textscan(fid, '%f %f %f %f %f %f');
beta12 = dat12{4};
g12 = dat12{5};
g12_err = dat12{6};
fclose(fid);


xfit=beta12;
gfit = g12;
gfit_err = g12_err;

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

X1(1:Ndat)=weight;
X2(1:Ndat)=xfit.*weight;
X3(1:Ndat)=(xfit.^2).*weight;
X4(1:Ndat)=(xfit.^3).*weight;
X5(1:Ndat)=(xfit.^4).*weight;
X6(1:Ndat)=(xfit.^5).*weight;
X7(1:Ndat)=(xfit.^6).*weight;
X8(1:Ndat)=(xfit.^7).*weight;
X9(1:Ndat)=(xfit.^8).*weight;


TEMP=[X1; X2; X3; X4; X5; X6; X7; X8; X9];

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p12=C*(A'*B');

p12_error=sqrt(diag(C));

p12'
p12_error'

chi12=(B-(A*p12)')*(B'-A*p12);
DoF = Ndat - 9; 
chi12_dof=chi12/DoF
Pvalue=gammainc(chi12/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp= p12(1) + p12(2)*xp + p12(3)*xp.^2 + p12(4)*xp.^3 + p12(5)*xp.^4 + p12(6)*xp.^5+ p12(7)*xp.^6+ p12(8)*xp.^7+ p12(9)*xp.^8;
plot(xp,yp,'-b','LineWidth',1);



xroot12(1:N)=0;
error12Tuned(1:N)=0;
g12Tuned(1:N)=0;
beta12Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) p12(1) + p12(2)*x + p12(3)*x^2 + p12(4)*x^3 + p12(5)*x^4 + p12(6)*x^5+ p12(7)*x^6+ p12(8)*x^7+ p12(9)*x^8 - g2Target,xInit);

xroot12(i)=xr;
beta12Tuned(i)=xroot12(i);

xpp=xroot12(i);
ypp=p12(1) + p12(2)*xpp + p12(3)*xpp^2 + p12(4)*xpp^3 + p12(5)*xpp^4  
            + p12(6)*xpp^5 + p12(7)*xpp^6 + p12(8)*xpp^7 + p12(9)*xpp^8 ;
xvec=[1; xpp; xpp^2; xpp^3; xpp^4;xpp^5;xpp^6;xpp^7;xpp^8];
% delta=ypp^2*sqrt(xvec'*(C*xvec));
delta=sqrt(xvec'*(C*xvec));
ypp;

error12Tuned(i)=delta;
g12Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot12,g12Tuned,error12Tuned,'bo','LineWidth',1);


fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L24_c0.30_ssc.dat.new3');
dat24 = textscan(fid, '%f %f %f %f %f %f');
beta24 = dat24{4};
g24 = dat24{5};
g24_err = dat24{6};
fclose(fid);

xfit=beta24;
gfit = g24;
gfit_err = g24_err;

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

clear Ndat X1 X2 X3 X4 X5 X6 X7 X8 X9

Ndat=length(xfit);

X1(1:Ndat)=weight;
X2(1:Ndat)=xfit.*weight;
X3(1:Ndat)=(xfit.^2).*weight;
X4(1:Ndat)=(xfit.^3).*weight;
X5(1:Ndat)=(xfit.^4).*weight;
X6(1:Ndat)=(xfit.^5).*weight;
X7(1:Ndat)=(xfit.^6).*weight;
X8(1:Ndat)=(xfit.^7).*weight;
X9(1:Ndat)=(xfit.^8).*weight;

TEMP=[X1; X2; X3; X4; X5; X6; X7; X8; X9];

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
DoF = Ndat - 9; 
chi24_dof=chi24/DoF
Pvalue=gammainc(chi24/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp2= p24(1) + p24(2)*xp + p24(3)*xp.^2 + p24(4)*xp.^3 + p24(5)*xp.^4 + p24(6)*xp.^5+ p24(7)*xp.^6+ p24(8)*xp.^7+ p24(9)*xp.^8;
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot12(i);
ypp=p24(1) + p24(2)*xpp + p24(3)*xpp^2 + p24(4)*xpp^3 + p24(5)*xpp^4 + p24(6)*xpp^5+ p24(7)*xpp^6+ p24(8)*xpp^7+ p24(9)*xpp^8;
xvec=[1; xpp; xpp^2; xpp^3; xpp^4; xpp^5; xpp^6; xpp^7; xpp^8];
% delta=ypp^2*sqrt(xvec'*(C*xvec));
delta=sqrt(xvec'*(C*xvec));
ypp;

error24(i)=delta;
g24stepped(i)=ypp;

end

% g24stepped'
% error24'

%errorbar(xroot12,g24stepped,error24,'ko','LineWidth',1);

text(3,5,['12^4 \rightarrow 24^4'],'fontsize',16)
text(3,4.5,['\chi^2/dof= ',num2str(chi12_dof,2),' and ',num2str(chi24_dof,2)],'fontsize',14)

xlim([2.5 4.2])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,7.4,'12^4','fontsize',16)
text(3.6,6.5,'24^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)