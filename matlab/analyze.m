function [p,p_error,chi2]= analyze(M,xfit,gfit,gfit_err)

weight=1./gfit_err;
%B=gfit./gfit_err;

Ndat=length(xfit);
B=(gfit.*weight)';
% if(Ndat == 4 && pflag == 1)
TEMP=tempPop(xfit,M,weight);

A=TEMP';

H=A'*A;

[U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); C = V*WW*U'; %C=inv(H) replaced by SVD;

p=C*(A'*B);

p_error=sqrt(diag(C));

chi2=(B-A*p)'*(B-A*p);

end