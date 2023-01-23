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
M=5;

% 12 -> 24

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7

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
DoF = Ndat - M; 
chi24_dof=chi24/DoF
Pvalue=gammainc(chi24/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp2= Vpn(p24,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot12(i);
ypp=Vpn(p24,xpp) ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
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

% 16 -> 32

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L16_c0.30_ssc.dat.new3');
dat16 = textscan(fid, '%f %f %f %f %f %f');
beta16 = dat16{4};
g16 = dat16{5};
g16_err = dat16{6};
fclose(fid);


xfit=beta16;
gfit = g16;
gfit_err = g16_err;

figure(122)
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

p16=C*(A'*B');

p16_error=sqrt(diag(C));

p16'
p16_error'

chi16=(B-(A*p16)')*(B'-A*p16);
DoF = Ndat - M; 
chi16_dof=chi16/DoF
Pvalue=gammainc(chi16/2,DoF/2)
Qvalue=1-Pvalue


xp=2.5:0.01:4.2;
yp= Vpn(p16,xp);
plot(xp,yp,'-b','LineWidth',1);



xroot16(1:N)=0;
error16Tuned(1:N)=0;
g16Tuned(1:N)=0;
beta16Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p16,x) - g2Target,xInit);

xroot16(i)=xr;
beta16Tuned(i)=xroot16(i);

xpp=xroot16(i);
ypp=Vpn(p16,xpp) ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error16Tuned(i)=delta;
g16Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot16,g16Tuned,error16Tuned,'bo','LineWidth',1);

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L32_c0.30_ssc.dat.new3');
dat32 = textscan(fid, '%f %f %f %f %f %f');
beta32 = dat32{4};
g32 = dat32{5};
g32_err = dat32{6};
fclose(fid);

xfit=beta32;
gfit = g32;
gfit_err = g32_err;

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

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p32=C*(A'*B');

p32_error=sqrt(diag(C));

p32'
p32_error'

chi32=(B-(A*p32)')*(B'-A*p32);
DoF = Ndat - M; 
chi32_dof=chi32/DoF
Pvalue=gammainc(chi32/2,DoF/2)
Qvalue=1-Pvalue



xp=2.5:0.01:4.2;
yp2= Vpn(p32,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot16(i);
ypp=Vpn(p32,xpp) ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error32(i)=delta;
g32stepped(i)=ypp;

end

% g32stepped'
% error32'

%errorbar(xroot16,g32stepped,error32,'ko','LineWidth',1);

text(3,5.5,['16^4 \rightarrow 32^4'],'fontsize',16)
text(3,5,['\chi^2/dof= ',num2str(chi16_dof,2),' and ',num2str(chi32_dof,2)],'fontsize',14)

xlim([2.5 4.2])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,7.7,'16^4','fontsize',16)
text(3.7,6.5,'32^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)

% 18 -> 36

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L18_c0.30_ssc.dat.new3');
dat18 = textscan(fid, '%f %f %f %f %f %f');
beta18 = dat18{4};
g18 = dat18{5};
g18_err = dat18{6};
fclose(fid);


xfit=beta18;
gfit = g18;
gfit_err = g18_err;

figure(123)
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

p18=C*(A'*B');

p18_error=sqrt(diag(C));

p18'
p18_error'

chi18=(B-(A*p18)')*(B'-A*p18);
DoF = Ndat - M; 
chi18_dof=chi18/DoF
Pvalue=gammainc(chi18/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp= Vpn(p18,xp);
plot(xp,yp,'-b','LineWidth',1);



xroot18(1:N)=0;
error18Tuned(1:N)=0;
g18Tuned(1:N)=0;
beta18Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p18,x) - g2Target,xInit);

xroot18(i)=xr;
beta18Tuned(i)=xroot18(i);

xpp=xroot18(i);
ypp=Vpn(p18,xpp);
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error18Tuned(i)=delta;
g18Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot18,g18Tuned,error18Tuned,'bo','LineWidth',1);

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L36_c0.30_ssc.dat.new3');
dat36 = textscan(fid, '%f %f %f %f %f %f');
beta36 = dat36{4};
g36 = dat36{5};
g36_err = dat36{6};
fclose(fid);

xfit=beta36;
gfit = g36;
gfit_err = g36_err;

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

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p36=C*(A'*B');

p36_error=sqrt(diag(C));

p36'
p36_error'

chi36=(B-(A*p36)')*(B'-A*p36);
DoF = Ndat - M; 
chi36_dof=chi36/DoF
Pvalue=gammainc(chi36/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp2= Vpn(p36,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot18(i);
ypp=Vpn(p36,xpp)  ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error36(i)=delta;
g36stepped(i)=ypp;

end

% g36stepped'
% error36'

%errorbar(xroot18,g36stepped,error36,'ko','LineWidth',1);


text(3,5,['18^4 \rightarrow 36^4'],'fontsize',16)
text(3,4.5,['\chi^2/dof= ',num2str(chi18_dof,2),' and ',num2str(chi36_dof,2)],'fontsize',14)


xlim([2.5 4.2])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,7.8,'18^4','fontsize',16)
text(3.75,6.4,'36^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)


% 20 -> 40

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L20_c0.30_ssc.dat.new3');
dat20 = textscan(fid, '%f %f %f %f %f %f');
beta20 = dat20{4};
g20 = dat20{5};
g20_err = dat20{6};
fclose(fid);


xfit=beta20;
gfit = g20;
gfit_err = g20_err;

figure(124)
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
M20_=5;
TEMP=tempPop(xfit,M20_,weight);

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p20=C*(A'*B');

p20_error=sqrt(diag(C));

p20'
p20_error'

chi20=(B-(A*p20)')*(B'-A*p20);
DoF = Ndat - M20_; 
chi20_dof=chi20/DoF
Pvalue=gammainc(chi20/2,DoF/2)
Qvalue=1-Pvalue


xp=2.5:0.01:4.2;
yp= Vpn(p20,xp);
plot(xp,yp,'-b','LineWidth',1);



xroot20(1:N)=0;
error20Tuned(1:N)=0;
g20Tuned(1:N)=0;
beta20Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p20,x) - g2Target,xInit);

xroot20(i)=xr;
beta20Tuned(i)=xroot20(i);

xpp=xroot20(i);
ypp=Vpn(p20,xpp)  ;
xvec=zeros(1,M20_);
t=1;
for j=1:M20_
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error20Tuned(i)=delta;
g20Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot20,g20Tuned,error20Tuned,'bo','LineWidth',1);

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L40_c0.30_ssc.dat.new3');
dat40 = textscan(fid, '%f %f %f %f %f %f');
beta40 = dat40{4};
g40 = dat40{5};
g40_err = dat40{6};
fclose(fid);

xfit=beta40;
gfit = g40;
gfit_err = g40_err;

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

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p40=C*(A'*B');

p40_error=sqrt(diag(C));

p40'
p40_error'

chi40=(B-(A*p40)')*(B'-A*p40);
DoF = Ndat - M; 
chi40_dof=chi40/DoF
Pvalue=gammainc(chi40/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp2=Vpn(p40,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot20(i);
ypp=Vpn(p40,xpp)  ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error40(i)=delta;
g40stepped(i)=ypp;

end

% g40stepped'
% error40'

%errorbar(xroot20,g40stepped,error40,'ko','LineWidth',1);

text(3,6,['20^4 \rightarrow 40^4'],'fontsize',16)
text(3,5.5,['\chi^2/dof= ',num2str(chi20_dof,2),' and ',num2str(chi40_dof,2)],'fontsize',14)


xlim([2.5 4.2])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,8,'20^4','fontsize',16)
text(3.8,6.2,'40^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)

% 24 -> 48

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L24_c0.30_ssc.dat.new3');
dat24 = textscan(fid, '%f %f %f %f %f %f');
beta24 = dat24{4};
g24 = dat24{5};
g24_err = dat24{6};
fclose(fid);


xfit=beta24;
gfit = g24;
gfit_err = g24_err;

figure(125)
box on
grid minor 
hold on

x = xfit;
y = gfit;
err = gfit_err;

errorbar(x,y,err,'ro','LineWidth',1);

weight=1./gfit_err;

B=gfit./gfit_err;

clear Ndat X1 X2 X3 X4 X5

Ndat=length(xfit);

TEMP=tempPop(xfit,M,weight);

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p24=C*(A'*B');

p24_error=sqrt(diag(C));

p24'
p24_error'

chi24=(B-(A*p24)')*(B'-A*p24);
DoF = Ndat - M; 
chi24_dof=chi24/DoF
Pvalue=gammainc(chi24/2,DoF/2)
Qvalue=1-Pvalue

xp=2.5:0.01:4.2;
yp= Vpn(p24,xp);
plot(xp,yp,'-b','LineWidth',1);



xroot24(1:N)=0;
error24Tuned(1:N)=0;
g24Tuned(1:N)=0;
beta24Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p24,x) - g2Target,xInit);

xroot24(i)=xr;
beta24Tuned(i)=xroot24(i);

xpp=xroot24(i);
ypp=Vpn(p24,xpp)  ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error24Tuned(i)=delta;
g24Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot24,g24Tuned,error24Tuned,'bo','LineWidth',1);

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L48_c0.30_ssc.dat.new3');
dat48 = textscan(fid, '%f %f %f %f %f %f');
beta48 = dat48{4};
g48 = dat48{5};
g48_err = dat48{6};
fclose(fid);

xfit=beta48;
gfit = g48;
gfit_err = g48_err;

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

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p48=C*(A'*B');

p48_error=sqrt(diag(C));

p48'
p48_error'

chi48=(B-(A*p48)')*(B'-A*p48);
DoF = Ndat - M; 
chi48_dof=chi48/DoF
Pvalue=gammainc(chi48/2,DoF/2)
Qvalue=1-Pvalue


xp=2.5:0.01:4.2;
% xp=2.:0.01:5;
yp2= Vpn(p48,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot24(i);
ypp=Vpn(p48,xpp)  ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error48(i)=delta;
g48stepped(i)=ypp;

end

% g48stepped'
% error48'

%errorbar(xroot24,g48stepped,error48,'ko','LineWidth',1);

text(3,6,['24^4 \rightarrow 48^4'],'fontsize',16)
text(3,5.5,['\chi^2/dof= ',num2str(chi24_dof,2),' and ',num2str(chi48_dof,2)],'fontsize',14)


xlim([2.5 4.2])
% xlim([2 5])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,8,'24^4','fontsize',16)
text(3.8,6.5,'48^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)



% 8 -> 16

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L8_c0.30_ssc.dat.new');
dat8 = textscan(fid, '%f %f %f %f %f %f');
beta8 = dat8{4};
g8 = dat8{5};
g8_err = dat8{6};
fclose(fid);


xfit=beta8;
gfit = g8;
gfit_err = g8_err;

figure(126)
box on
grid minor 
hold on

x = xfit;
y = gfit;
err = gfit_err;

errorbar(x,y,err,'ro','LineWidth',1);

weight=1./gfit_err;

B=gfit./gfit_err;

clear Ndat X1 X2 X3 X4 X5 X6

Ndat=length(xfit);
M8=5
TEMP=tempPop(xfit,M8,weight);

clear A B C H U S V w winv WW

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p8=C*(A'*B');

p8_error=sqrt(diag(C));

p8'
p8_error'

chi8=(B-(A*p8)')*(B'-A*p8);
DoF = Ndat - M8; 
chi8_dof=chi8/DoF
Pvalue=gammainc(chi8/2,DoF/2)
Qvalue=1-Pvalue

xp=2.4:0.01:4.2;
yp= Vpn(p8,xp);
plot(xp,yp,'-b','LineWidth',1);

xroot8(1:N)=0;
error8Tuned(1:N)=0;
g8Tuned(1:N)=0;
beta8Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p8,x) - g2Target,xInit);

xroot8(i)=xr;
beta8Tuned(i)=xroot8(i);

xpp=xroot8(i);
ypp=Vpn(p8,xpp) ;
xvec=zeros(1,M8);
t=1;
for j=1:M8
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error8Tuned(i)=delta;
g8Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot8,g8Tuned,error8Tuned,'bo','LineWidth',1);

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L16_c0.30_ssc.dat.new3');
dat16 = textscan(fid, '%f %f %f %f %f %f');
beta16 = dat16{4};
g16 = dat16{5};
g16_err = dat16{6};
fclose(fid);

xfit=beta16;
gfit = g16;
gfit_err = g16_err;

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

p16=C*(A'*B');

p16_error=sqrt(diag(C));

p16'
p16_error'

chi16=(B-(A*p16)')*(B'-A*p16);
DoF = Ndat - M; 
chi16_dof=chi16/DoF
Pvalue=gammainc(chi16/2,DoF/2)
Qvalue=1-Pvalue


xp=2.4:0.01:4.2;
% xp=2.:0.01:5;
yp2= Vpn(p16,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot8(i);
ypp=Vpn(p16,xpp);
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error16(i)=delta;
g16stepped(i)=ypp;

end

% g16stepped'
% error16'

%errorbar(xroot8,g16stepped,error16,'ko','LineWidth',1);

text(3,6,['8^4 \rightarrow 16^4'],'fontsize',16)
text(3,5.5,['\chi^2/dof= ',num2str(chi8_dof,2),' and ',num2str(chi16_dof,2)],'fontsize',14)


xlim([2.4 4.2])
% xlim([2 5])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,8,'8^4','fontsize',16)
text(3.8,6.5,'16^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)


% 10 -> 20

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L10_c0.30_ssc.dat.new');
dat10 = textscan(fid, '%f %f %f %f %f %f');
beta10 = dat10{4};
g10 = dat10{5};
g10_err = dat10{6};
fclose(fid);

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7


xfit=beta10;
gfit = g10;
gfit_err = g10_err;

figure(127)
box on
grid minor 
hold on

x = xfit;
y = gfit;
err = gfit_err;

errorbar(x,y,err,'ro','LineWidth',1);

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

p10=C*(A'*B');

p10_error=sqrt(diag(C));

p10'
p10_error'

chi10=(B-(A*p10)')*(B'-A*p10);
DoF = Ndat - M; 
chi10_dof=chi10/DoF
Pvalue=gammainc(chi10/2,DoF/2)
Qvalue=1-Pvalue

xp=2.4:0.01:4.2;
yp= Vpn(p10,xp);
plot(xp,yp,'-b','LineWidth',1);

resid10 = y - p10(1) - p10(2)*x - p10(3)*x.^2 - p10(4)*x.^3 - p10(5)*x.^4;


xroot10(1:N)=0;
error10Tuned(1:N)=0;
g10Tuned(1:N)=0;
beta10Tuned(1:N)=0;

for i=1:N
    
g2Target=4.+(i-1)*0.1;

xInit=3.3;


xr=fzero(@(x) Vpn(p10,x) - g2Target,xInit);

xroot10(i)=xr;
beta10Tuned(i)=xroot10(i);

xpp=xroot10(i);
ypp=Vpn(p10,xpp) ;
xvec=zeros(1,M);
t=1;
for j=1:M
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error10Tuned(i)=delta;
g10Tuned(i)=ypp;

end

% beta12Tuned'
% error12Tuned'

%errorbar(xroot10,g10Tuned,error10Tuned,'bo','LineWidth',1);

fid = fopen('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L20_c0.30_ssc.dat.new3');
dat20 = textscan(fid, '%f %f %f %f %f %f');
beta20 = dat20{4};
g20 = dat20{5};
g20_err = dat20{6};
fclose(fid);

clear xbeta gfit gfit_err xfit weight B Ndat X1 X2 X3 X4 X5 X6 X7

xfit=beta20;
gfit = g20;
gfit_err = g20_err;

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

clear Ndat X1 X2 X3 X4

Ndat=length(xfit);

TEMP=tempPop(xfit,M20_,weight);

clear A B C H U S V w winv WW

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD; 

B=(gfit.*weight)';

p20=C*(A'*B');

p20_error=sqrt(diag(C));

p20'
p20_error'

chi20=(B-(A*p20)')*(B'-A*p20);
DoF = Ndat - M20_; 
chi20_dof=chi20/DoF
Pvalue=gammainc(chi20/2,DoF/2)
Qvalue=1-Pvalue


xp=2.4:0.01:4.2;
% xp=2.:0.01:5;
yp2= Vpn(p20,xp);
plot(xp,yp2,'-c','LineWidth',1);

for i=1:N

xpp=xroot10(i);
ypp=Vpn(p20,xpp) ;
xvec=zeros(1,M20_);
t=1;
for j=1:M20_
    xvec(j)=t;
    t=t*xpp;
end
delta=sqrt(xvec*(C*xvec'));
ypp;

error20(i)=delta;
g20stepped(i)=ypp;

end

% g20stepped'
% error20'

%errorbar(xroot10,g20stepped,error20,'ko','LineWidth',1);

text(3,6,['10^4 \rightarrow 20^4'],'fontsize',16)
text(3,5.5,['\chi^2/dof= ',num2str(chi10_dof,2),' and ',num2str(chi20_dof,2)],'fontsize',14)


xlim([2.4 4.2])
% xlim([2 5])
% ylim([0 9])
xlabel('{\fontsize{16}\beta}','fontsize',14)
ylabel('{ g^2(L)}','fontsize',14)
text(3,8,'10^4','fontsize',16)
text(3.8,6.5,'20^4','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('N_f = 10   SSC  c = 0.30  interpolation','fontsize',14)


figure(227)
box on
grid minor 
hold on

xlim([2.4 4.2])
ylim([-0.1 0.1])

resid20 = y - p20(1) - p20(2)*x - p20(3)*x.^2 - p20(4)*x.^3 - p20(5)*x.^4;
errorbar(beta10,resid10,g10_err,'bo','LineWidth',1);
errorbar(beta20,resid20,g20_err,'ro','LineWidth',1);
stepped=[g48stepped; g32stepped; g36stepped; g40stepped; g24stepped; g20stepped; g16stepped];
spouseErrH=[error48; error32; error36; error40; error24; error20; error16];
spouseErrL=[error24Tuned; error16Tuned; error18Tuned; error20Tuned; error12Tuned; error10Tuned; error8Tuned];

%%
%%%%%%%%%%%%%%%%%%
scale=2.0;
scale2=scale^2;

% FitIndex=[2 3 4];
% FitIndex=[2 3 4 5];
% FitIndex=[1 2 3 4 5];
% FitIndex=[1 2 3 5];
% FitIndex=[1 2 3 4 5 6 7];
FitIndex=[1 2 3 5 6 7];
% FitIndex=[1 2 3 5 7];

offset = 1.2;

pflag = 1;

for indx=1:N
% for indx=6:46
% for indx=31:31

clear xfit XX1 XX2 XX3  

g2=4.+ 0.1*(indx-1);
        
xL(1)=1/12^2;
xL(2)=1/16^2;
xL(3)=1/18^2;
xL(4)=1/20^2;
xL(5)=1/24^2;
xL(6)=1/8^2;
xL(7)=1/10^2;


sigma(1)=g24stepped(indx);
sigma(2)=g32stepped(indx);
sigma(3)=g36stepped(indx);
sigma(4)=g40stepped(indx);
sigma(5)=g48stepped(indx);
sigma(6)=g16stepped(indx);
sigma(7)=g20stepped(indx);


sigma_err(1)=sqrt( error24(indx)^2 + error12Tuned(indx)^2 );
sigma_err(2)=sqrt( error32(indx)^2 + error16Tuned(indx)^2 );
sigma_err(3)=sqrt( error36(indx)^2 + error18Tuned(indx)^2 );
sigma_err(4)=sqrt( error40(indx)^2 + error20Tuned(indx)^2 );    
sigma_err(5)=sqrt( error48(indx)^2 + error24Tuned(indx)^2 );    
sigma_err(6)=sqrt( error16(indx)^2 + error8Tuned(indx)^2 );    
sigma_err(7)=sqrt( error20(indx)^2 + error10Tuned(indx)^2 );    

xfit=xL(FitIndex);
gfit = sigma(FitIndex);
gfit_err = sigma_err(FitIndex);
weight=1./gfit_err;

B=gfit./gfit_err;

Ndat=length(xfit);

XX1(1:Ndat)=weight;
XX2(1:Ndat)=xfit.*weight;
XX3(1:Ndat)=(xfit.^2).*weight;

% if(Ndat == 4 && pflag == 1)
if(pflag == 1)
 TEMP=[XX1; XX2; XX3];
else
 TEMP=[XX1; XX2];   
end

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD;

B=(gfit.*weight)';

p=C*(A'*B);

p_error=sqrt(diag(C));

Nf=10;

b0=(11-2*Nf/3)/(4*pi)^2;
b1=(102-38*Nf/3)/(4*pi)^4;

s0=2*b0*log(scale);
s1=(2*b0*log(scale))^2 + 2*b1*log(scale);

yp  = g2 + s0*g2^2 + s1*g2^3;
ypp = g2 + s0*g2^2;

g2target(indx)   = g2;
g2step(indx)     = p(1);
g2step_err(indx) = p_error(1);
g2_2loop(indx)   = yp;
g2_1loop(indx)   = ypp;

end


figure(100)
box on
grid minor
hold on

xp=0.0:0.00001:15.;
yp  = xp + s0*xp.^2 + s1*xp.^3;
ypp = xp + s0*xp.^2;

errorbar(g2target,(g2step - g2target)./log(scale2) ,g2step_err./log(scale2),'ro','LineWidth',1);
hold on

% plot(xp,(ypp - xp)./log(scale2) ,'m-','LineWidth',1);
% plot(xp,(yp - xp)./log(scale2) ,'b-','LineWidth',1);
% hold on
xb=linspace(0,max(g2target),10000)
yb=p_beta(10,xb);
plot(xb,yb./log(scale2),'k-');

xlim([3.5  10.5]);
ylim([0.2  1.3]);


% hleg1=legend('beta function','1 loop','2 loop');
% %set(hleg1,'Location','SouthEast')
% set(hleg1,'Location','NorthWest')
% set(hleg1,'Interpreter','none')


xlabel('\fontsize{16}g^2(L)','fontsize',14)
ylabel('{ ( g^2(sL) - g^2(L) )/log(s^2) }','fontsize',14)
% text(1,0.25,'C = 0.2','fontsize',16)
% title('SSC gradient flow with unimproved s=1.5 step scaling function','fontsize',14)
title('SSC  c = 0.30  s = 0.75   beta function','fontsize',14)

for i = 1:N
fprintf('%6.3f %12.6f %12.6f\n',g2target(i),(g2step(i) - g2target(i))/log(scale2) ,g2step_err(i)/log(scale2));
end



% print -depsc2 beta_fn_ssc_c0.30.eps

%%


% DoF=3;

scale=2.0;

clear xfit XX1 XX2 XX3

% FitIndex=[2 3 4];
% FitIndex=[2 3 4 5];
% FitIndex=[1 2 3 4 5];
% FitIndex=[1 2 3 5];
% FitIndex=[1 2 3 4 5 6 7];
% FitIndex=[1 2 3 5 6 7];
% FitIndex=[1 2 3 5 7];
% FitIndex=[1 2 3 5];
% FitIndex=[2 3 5];

datflag = 4;

a4 = 1;

offset = 0.4;

if (a4 == 1)
  FitIndex=[1 2 3 4 5 6 7];
  pflag = 1;
else
  FitIndex=[2 3 5];
  pflag = 0;
end

gsqrd(1:N) = 0.0;
slope(1:N) = 0.0;
slope_err(1:N) = 0.0;

% for indx=41:41;
% for indx=6:46;
for indx=6:2:41;

g2=4.+ 0.1*(indx-1); 

gsqrd(indx) = g2;

figure(indx)
box on
grid minor
hold on
    
xL(1)=1/12^2;
xL(2)=1/16^2;
xL(3)=1/18^2;
xL(4)=1/20^2;
xL(5)=1/24^2;
xL(6)=1/8^2;
xL(7)=1/10^2;

sigma(1)=g24stepped(indx);
sigma(2)=g32stepped(indx);
sigma(3)=g36stepped(indx);
sigma(4)=g40stepped(indx);
sigma(5)=g48stepped(indx);
sigma(6)=g16stepped(indx);
sigma(7)=g20stepped(indx);

sigma_err(1)=sqrt( error24(indx)^2 + error12Tuned(indx)^2 );
sigma_err(2)=sqrt( error32(indx)^2 + error16Tuned(indx)^2 );
sigma_err(3)=sqrt( error36(indx)^2 + error18Tuned(indx)^2 );
sigma_err(4)=sqrt( error40(indx)^2 + error20Tuned(indx)^2 );    
sigma_err(5)=sqrt( error48(indx)^2 + error24Tuned(indx)^2 );    
sigma_err(6)=sqrt( error16(indx)^2 + error8Tuned(indx)^2 );    
sigma_err(7)=sqrt( error20(indx)^2 + error10Tuned(indx)^2 );    

xfit=xL(FitIndex);
gfit = sigma(FitIndex);
gfit_err = sigma_err(FitIndex);
weight=1./gfit_err;

% DoF = length(xfit)-2

xfit
gfit
gfit_err


B=gfit./gfit_err;

Ndat=length(xfit);

XX1(1:Ndat)=weight;
XX2(1:Ndat)=xfit.*weight;
XX3(1:Ndat)=(xfit.^2).*weight;

% if(Ndat == 4 && pflag == 1)
if(pflag == 1)
 TEMP=[XX1; XX2; XX3];
 DoF = length(xfit)-3
else
 TEMP=[XX1; XX2]; 
 DoF = length(xfit)-2
end

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD;

B=(gfit.*weight)';

p=C*(A'*B);

p_error=sqrt(diag(C));

chi2=(B-A*p)'*(B-A*p);

chi2_dof=chi2/DoF

Pvalue=gammainc(chi2/2,DoF/2)
Qvalue=1-Pvalue


xp=-0.0001:0.00001:0.007;
xp=-0.0001:0.00001:0.017;

% if(Ndat == 4 && pflag == 1)
if(pflag == 1)
 yp= p(1)+p(2)*xp + p(3)*xp.^2;
else
 yp= p(1)+p(2)*xp;
end


slope(indx) = p(2);
slope_err(indx) = p_error(2);


xlim([-0.001  0.008 ]);
xlim([-0.001  0.018 ]);

yfit=(p(1)-g2)/log(2.^2);
% ymin = yfit - 1.5*abs((sigma(1)-g2)/log(1.5^2));
% ymax = yfit + 1.5*abs((sigma(1)-g2)/log(1.5^2));
% ymin = yfit - 1.0*abs((sigma(1)-g2)/log(1.5^2));
% ymax = yfit + 1.0*abs((sigma(1)-g2)/log(1.5^2));
% ymin = yfit - 4.5*abs((sigma(1)-g2)/log(1.5^2));
% ymax = yfit + 1.0*abs((sigma(1)-g2)/log(1.5^2));
ymin = 0.2;
ymax = 1.1;
ylim([ymin  ymax]);

Nf=10;
scale=2.;

b0=(11-2*Nf/3)/(4*pi)^2;
b1=(102-38*Nf/3)/(4*pi)^4;

s0=2*b0*log(scale);
s1=(2*b0*log(scale))^2 + 2*b1*log(scale);

p_2loop  = g2 + s0*g2^2 + s1*g2^3;
p_1loop  = g2 + s0*g2^2;

% plot(-0.0002,(p_2loop-g2)/log(2.^2),'bo','LineWidth',1);
% plot(-0.0002,(p_1loop-g2)/log(2.^2),'mo','LineWidth',1);

% hleg1=legend('2 loop','1 loop');
% set(hleg1,'Location','SouthWest')
% set(hleg1,'Interpreter','none')

if (a4 == 1)
    errorbar(0.0000,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'mo','LineWidth',1);
else
    errorbar(-5e-4,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'co','LineWidth',1);
end

% if(Ndat > 4)
%    errorbar(xL,(sigma-g2)/log(2.^2),sigma_err/log(2.^2),'ro','LineWidth',1);
% elseif(Ndat == 4)
%    if(pflag == 4)
%        errorbar(xL(2:5),(sigma(2:5)-g2)/log(2.^2),sigma_err(2:5)/log(2.^2),'ro','LineWidth',1);
%        errorbar(xL(1),(sigma(1)-g2)/log(2.^2),sigma_err(1)/log(2.^2),'co','LineWidth',1);
%    elseif(pflag == 5)
%         errorbar(xL(1:3),(sigma(1:3)-g2)/log(2.^2),sigma_err(1:3)/log(2.^2),'ro','LineWidth',1);
%         errorbar(xL(5),(sigma(5)-g2)/log(2.^2),sigma_err(5)/log(2.^2),'ro','LineWidth',1);
%         errorbar(xL(4),(sigma(4)-g2)/log(2.^2),sigma_err(4)/log(2.^2),'co','LineWidth',1);
%    end
% elseif(Ndat == 3)
%    errorbar(xL(2:4),(sigma(2:4)-g2)/log(2.^2),sigma_err(2:4)/log(2.^2),'ro','LineWidth',1);    
%    errorbar(xL(1),(sigma(1)-g2)/log(2.^2),sigma_err(1)/log(2.^2),'co','LineWidth',1);    
% end
if (datflag == 4)
    errorbar(xL(1:3),(sigma(1:3)-g2)/log(2.^2),sigma_err(1:3)/log(2.^2),'ro','LineWidth',1);
    errorbar(xL(5:7),(sigma(5:7)-g2)/log(2.^2),sigma_err(5:7)/log(2.^2),'ro','LineWidth',1);
    errorbar(xL(4),(sigma(4)-g2)/log(2.^2),sigma_err(4)/log(2.^2),'co','LineWidth',1);
else
    errorbar(xL,(sigma-g2)/log(2.^2),sigma_err/log(2.^2),'ro','LineWidth',1);
end   
    
if (a4 == 1)
    plot(xp,(yp-g2)/log(2.^2),'-b','LineWidth',1);
else
    plot(xp,(yp-g2)/log(2.^2),'--c','LineWidth',1);
end

xlabel('\fontsize{16} a^2/L^2')
ylabel('{ (g^2(sL) - g^2(L))/log(s^2) }','fontsize',16)

% if(p(2) > 0 )
%   text(0.003,ymin + 0.3*(yfit-ymin),['g^2 = ',num2str(g2)],'fontsize',14)
%   text(0.003,ymin + 0.9*(yfit-ymin),[ 'c_{fit}       = ',num2str((p(1)-g2)/log(2.^2),3),' \pm ',num2str(p_error(1)/log(2.^2),2), ],'fontsize',14)
%   text(0.003,ymin + 0.6*(yfit-ymin),[ 'slope  = ',num2str(p(2)/log(2.^2),3),' \pm ',num2str(p_error(2)/log(2.^2),2), ],'fontsize',14)
%   text(0.003,ymin + 0.75*(yfit-ymin),[ 'c_{2loop} = ',num2str((p_2loop-g2)/log(2.^2),3), ],'fontsize',14)
%   text(0.003,ymin + 0.45*(yfit-ymin),[ '\chi^2/dof= ',num2str(chi2_dof,2) ],'fontsize',14)
% else
%   text(0.003,ymin + 0.3*(yfit-ymin),['g^2 = ',num2str(g2)],'fontsize',14)
%   text(0.003,ymax - 0.2*(ymax-yfit),['c_{fit}       = ',num2str((p(1)-g2)/log(2.^2),3),' \pm ',num2str(p_error(1)/log(2.^2),2), ],'fontsize',14)
%   text(0.003,ymax - 0.5*(ymax-yfit),['slope  = ',num2str(p(2)/log(2.^2),3),' \pm ',num2str(p_error(2)/log(2.^2),2), ],'fontsize',14)
%   text(0.003,ymax - 0.35*(ymax-yfit),['c_{2loop} = ',num2str((p_2loop-g2)/log(2.^2),3), ],'fontsize',14)
%   text(0.003,ymax - 0.65*(ymax-yfit),['\chi^2/dof= ',num2str(chi2_dof,2) ],'fontsize',14)
% end

if (a4 == 1)
text(0.006,ymin + 0.9*(ymax-ymin),['g^2 = ',num2str(g2)],'fontsize',14)
text(0.006,ymin + 0.8*(ymax-ymin),['c_0       = ',num2str((p(1)-g2)/log(2.^2),3),' \pm ',num2str(p_error(1)/log(2.^2),2), ],'fontsize',14)
text(0.006,ymin + 0.7*(ymax-ymin),['slope  = ',num2str(p(2)/log(2.^2),3),' \pm ',num2str(p_error(2)/log(2.^2),2), ],'fontsize',14)
% text(0.003,ymax - 0.35*(ymax-yfit),['c_{2loop} = ',num2str((p_2loop-g2)/log(2.^2),3), ],'fontsize',14)
text(0.006,ymin + 0.6*(ymax-ymin),['\chi^2/dof= ',num2str(chi2_dof,2) ],'fontsize',14)

title('N_f = 10  c=0.30  SSC  s=2  beta function','fontsize',16)
end

% eval(['print -depsc2 Fit_ssc/g2_',num2str(indx),'.eps']);
% eval(['print -depsc2 Fit_ssc_c0.3_L8.12.16.20.24/g2_',num2str(1.+0.1*(indx-1)),'.eps']);
% eval(['print -depsc2 Fit_ssc_c0.2_L16.20.24_strong_bracket2/g2_v2_',num2str(1.+0.1*(indx-1)),'.eps']);
% eval(['print -depsc2 Fit_ssc_strong/g2_',num2str(1.+0.1*(indx-1)),'.eps']);

end

% figure(314)
% box on
% grid minor
% hold on
% 
% title('N_f = 10  c=0.30  SSC  s=2  beta function','fontsize',16)
% ylabel('\fontsize{16} extrapolation slope')
% xlabel('{ g^2 }','fontsize',16)
% xlim([4  9]);
% ylim([-150  150]);
% errorbar(gsqrd,slope,slope_err,'ro','LineWidth',1);


%%

save_stuff = [xL; (sigma-g2)/log(4.); sigma_err/log(4.)]';
save_stuff2 = [0; (p(1)-g2)/log(4.); p_error(1)/log(4.)]';
save_stuff3 = [xp;(yp-g2)/log(2.^2)]';


% save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8_quad.dat  save_stuff -ascii
% save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8_continuum_quad.dat  save_stuff2 -ascii
% save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8_fit_quad.dat  save_stuff3 -ascii

save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8.dat  save_stuff -ascii
save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8_continuum.dat  save_stuff2 -ascii
save  /Users/khollan1/talks/lattice2018/dani/nf10_ssc_c030_gsq8_fit.dat  save_stuff3 -ascii

%%
