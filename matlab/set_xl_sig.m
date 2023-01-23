function set_xl_sig(indx,stepp,Herr,Lerr)
global xL;
xL(1)=1/24^2;
xL(2)=1/16^2;
xL(3)=1/18^2;
xL(4)=1/20^2;
xL(5)=1/12^2;
xL(6)=1/10^2;
xL(7)=1/8^2;

global sigma;
sigma(1)=stepp(1,indx);
sigma(2)=stepp(2,indx);
sigma(3)=stepp(3,indx);
sigma(4)=stepp(4,indx);
sigma(5)=stepp(5,indx);
sigma(6)=stepp(6,indx);
sigma(7)=stepp(7,indx);


global sigma_err;
sigma_err(1)=sqrt( Herr(1,indx)^2 + Lerr(1,indx)^2 ); 
sigma_err(2)=sqrt( Herr(2,indx)^2 + Lerr(2,indx)^2 );
sigma_err(3)=sqrt( Herr(3,indx)^2 + Lerr(3,indx)^2 );
sigma_err(4)=sqrt( Herr(4,indx)^2 + Lerr(4,indx)^2 ); 
sigma_err(5)=sqrt( Herr(5,indx)^2 + Lerr(5,indx)^2 );   
sigma_err(6)=sqrt( Herr(6,indx)^2 + Lerr(6,indx)^2 );
sigma_err(7)=sqrt( Herr(7,indx)^2 + Lerr(7,indx)^2 ); 
end