


% DoF=3;

scale=2.0;

clear xfit XX1 XX2 XX3 M

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



pflag = 1;


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
OFI=[1 2 3 4 5 6 7]
for i=1:3
switch i
    case 1
        FitIndex=OFI;
        M=3;
    case 2 
        FitIndex=OFI(1:5);
        M=2;
    otherwise
        M=2;
        FitIndex=OFI(1:end-i);       
end

xfit=xL(FitIndex);
gfit = sigma(FitIndex);
gfit_err = sigma_err(FitIndex);
weight=1./gfit_err;

% DoF = length(xfit)-2

%xfit
%gfit
%gfit_err
DoF = length(xfit)-M;
[p,p_error,chi2]=analyze(M,xfit,gfit,gfit_err);
chi2_dof=chi2/DoF;
Pvalue=gammainc(chi2/2,DoF/2);
Qvalue=1-Pvalue


xp=-0.0001:0.00001:0.007;
xp=-0.0001:0.00001:0.017;

% if(Ndat == 4 && pflag == 1)
yp=Vpn(p,xp);

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
%{
Nf=10;
scale=2.;

b0=(11-2*Nf/3)/(4*pi)^2;
b1=(102-38*Nf/3)/(4*pi)^4;

s0=2*b0*log(scale);
s1=(2*b0*log(scale))^2 + 2*b1*log(scale);

p_2loop  = g2 + s0*g2^2 + s1*g2^3;
p_1loop  = g2 + s0*g2^2;
%}
% plot(-0.0002,(p_2loop-g2)/log(2.^2),'bo','LineWidth',1);
% plot(-0.0002,(p_1loop-g2)/log(2.^2),'mo','LineWidth',1);

% hleg1=legend('2 loop','1 loop');
% set(hleg1,'Location','SouthWest')
% set(hleg1,'Interpreter','none')

% This plots the predicted point
if (i == 1)
    errorbar(0.0000,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'mo','LineWidth',1);
elseif(i==2)
    errorbar(0.0000,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'co','LineWidth',1);
    %errorbar(-5e-4,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'co','LineWidth',1);
else
    errorbar(0.0000,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'LineWidth',1);
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

if (i == 1)
    plot(xp,(yp-g2)/log(2.^2),'-b','LineWidth',1);
else
    plot(xp,(yp-g2)/log(2.^2),'--c','LineWidth',1);
end

end

if (datflag == 4)
    errorbar(xL(1:3),(sigma(1:3)-g2)/log(2.^2),sigma_err(1:3)/log(2.^2),'ro','LineWidth',1);
    errorbar(xL(5:7),(sigma(5:7)-g2)/log(2.^2),sigma_err(5:7)/log(2.^2),'ro','LineWidth',1);
    errorbar(xL(4),(sigma(4)-g2)/log(2.^2),sigma_err(4)/log(2.^2),'co','LineWidth',1);
else
    errorbar(xL,(sigma-g2)/log(2.^2),sigma_err/log(2.^2),'ro','LineWidth',1);
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

title(sprintf('N_f = 10  c=0.30  SSC  s=%g  beta function',S),'fontsize',16)
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