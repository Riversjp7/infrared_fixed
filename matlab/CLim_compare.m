function [p,p_error,chi2] = CLim_compare(FI,m,indx)

[xL,sigma,sigma_err]=get_xl_sig;
scale=2.0;

clear xfit XX1 XX2 XX3 M

%indx
g2=4.+ 0.1*(indx-1) %indx should be returned by slider

  
OFI=[1 2 3 4 5 6 7];

FitIndex=OFI(1:FI);
M=m; %this M is for quadratic vs linear fits. NOT M for previous section

xfit=xL(FitIndex);
gfit = sigma(FitIndex);
gfit_err = sigma_err(FitIndex);
weight=1./gfit_err;

DoF = length(xfit)-M;
[p,p_error,chi2]=analyze(M,xfit,gfit,gfit_err);
chi2_dof=chi2/DoF;
Pvalue=gammainc(chi2/2,DoF/2);
Qvalue=1-Pvalue;



end