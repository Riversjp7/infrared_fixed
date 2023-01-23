function bet = p_beta(nf,grange)
% This function returns y as the beta function. it takes in
% 2 parameters the number of flavors and the  g^2 range.
   format long
    %This switch allows for default arguments
   switch nargin
       case 0
           nf=10;
           brange=linspace(0,1,1000);
       case 1
           brange=linspace(0,1,1000);
   end
   %beta coefficients
   b0=(0.25)*(11-2*nf/3);
    b1=(0.25^2)*(102-38*nf/3);
	b2=(0.25^3)*(2857/2-5033*nf/18+325*(nf^2)/54);
	b3=(0.25^4)*(149753/6+3564*zeta(3)-(1078361/162+6508*zeta(3)/27)*nf...
        +(50056/162+6472*zeta(3)/81)*nf^2 +(1093/729)*nf^3);
	b4=(0.25^5)*(8157455/16+621885*zeta(3)/2-88209*zeta(4)/2-288090*zeta(5)...
					 +(-336460813/1944-4811164*zeta(3)/81+33935*zeta(4)/6 + 1358995*zeta(5)/27)*nf...
					 +(25960913/1944+698531*zeta(3)/81-10526*zeta(4)/9-381760*zeta(5)/81)*nf^2 ...
					 +(-630559/5832-48722*zeta(3)/243+1618*zeta(4)/27+460*zeta(5)/9)*nf^3 ...
					 +(1205/2916-152*zeta(3)/81)*nf^4);
    bCo=[0,0,b0,b1,b2,b3,b4];
    %must be modified for beta(g^2) rather than beta(a)
    gCo=[4*pi^2,1,1/(4*pi^2),1/(4*pi^2)^2,1/(4*pi^2)^3,1/(4*pi^2)^4,1/4*(pi^2)^5];
    bet=Vpn((bCo.*gCo),grange);
    %plot(brange,bet,'k-');
end