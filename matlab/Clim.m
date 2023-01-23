classdef Clim
    properties  
     indx=[];
     g2=[];
     M=[];
     Q=[];
     C=[];
     C_inv=[];
     HAT=[];
     chi2=[];
     OFI=[];
     FI=[];
     stepp=[];
     HPerr=[];
     LPerr=[];
     DoF=[];
     c0=[];
     p=[];
     p_error=[];
    xL=[1/24^2 1/16^2 1/18^2 1/20^2 1/12^2 1/10^2 1/8^2];
    sigma=[];
    sigma_err=[];
    end
    methods
        function obj = Clim(step,herr,lerr)
            obj.indx=1;
            obj.stepp=step;
            obj.HPerr=herr;
            obj.LPerr=lerr;
            obj.OFI=1:7;
            obj.sigma(1)=obj.stepp(1,obj.indx);
            obj.sigma(2)=obj.stepp(2,obj.indx);
            obj.sigma(3)=obj.stepp(3,obj.indx);
            obj.sigma(4)=obj.stepp(4,obj.indx);
            obj.sigma(5)=obj.stepp(5,obj.indx);
            obj.sigma(6)=obj.stepp(6,obj.indx);
            obj.sigma(7)=obj.stepp(7,obj.indx);
            obj.sigma_err(1)=sqrt(obj.HPerr(1,obj.indx)^2 + obj.LPerr(1,obj.indx)^2 ); 
            obj.sigma_err(2)=sqrt(obj.HPerr(2,obj.indx)^2 + obj.LPerr(2,obj.indx)^2 );
            obj.sigma_err(3)=sqrt(obj.HPerr(3,obj.indx)^2 + obj.LPerr(3,obj.indx)^2 );
            obj.sigma_err(4)=sqrt(obj.HPerr(4,obj.indx)^2 + obj.LPerr(4,obj.indx)^2 ); 
            obj.sigma_err(5)=sqrt(obj.HPerr(5,obj.indx)^2 + obj.LPerr(5,obj.indx)^2 );   
            obj.sigma_err(6)=sqrt(obj.HPerr(6,obj.indx)^2 + obj.LPerr(6,obj.indx)^2 );
            obj.sigma_err(7)=sqrt(obj.HPerr(7,obj.indx)^2 + obj.LPerr(7,obj.indx)^2 ); 
        obj.g2=4.+ 0.1*((obj.indx)-1);
        obj.FI= obj.OFI;
        obj=lin_reg(obj);
        end
        function obj=indx_change(obj,indx)
            obj.indx=indx;
            obj.sigma(1)=obj.stepp(1,obj.indx);
            obj.sigma(2)=obj.stepp(2,obj.indx);
            obj.sigma(3)=obj.stepp(3,obj.indx);
            obj.sigma(4)=obj.stepp(4,obj.indx);
            obj.sigma(5)=obj.stepp(5,obj.indx);
            obj.sigma(6)=obj.stepp(6,obj.indx);
            obj.sigma(7)=obj.stepp(7,obj.indx);
            obj.sigma_err(1)=sqrt(obj.HPerr(1,obj.indx)^2 + obj.LPerr(1,obj.indx)^2 ); 
            obj.sigma_err(2)=sqrt(obj.HPerr(2,obj.indx)^2 + obj.LPerr(2,obj.indx)^2 );
            obj.sigma_err(3)=sqrt(obj.HPerr(3,obj.indx)^2 + obj.LPerr(3,obj.indx)^2 );
            obj.sigma_err(4)=sqrt(obj.HPerr(4,obj.indx)^2 + obj.LPerr(4,obj.indx)^2 ); 
            obj.sigma_err(5)=sqrt(obj.HPerr(5,obj.indx)^2 + obj.LPerr(5,obj.indx)^2 );   
            obj.sigma_err(6)=sqrt(obj.HPerr(6,obj.indx)^2 + obj.LPerr(6,obj.indx)^2 );
            obj.sigma_err(7)=sqrt(obj.HPerr(7,obj.indx)^2 + obj.LPerr(7,obj.indx)^2 ); 
            obj.g2=4.+ 0.1*((obj.indx)-1);
        end
        function obj=fi_block(obj,num1)
                            temp=obj.FI;
                             obj.FI= temp;
                             obj.FI(num1)=[];
        end
        function obj=change_sig(obj, sig, sig_err)
            obj.sigma=sig;
            obj.sigma_err=sig_err;
        end
       function obj=lin_reg(obj)
            xfit=obj.xL(obj.FI);
            gfit = obj.sigma(obj.FI);
            gfit_err = obj.sigma_err(obj.FI);
            weight=1./gfit_err;
            Ndat=length(xfit);
            B=(gfit.*weight)';

            TEMP=obj.tempPop(xfit,obj.M,weight);

            A=TEMP';

            H=A'*A;
            obj.C_inv=H;
            obj.HAT=(TEMP./weight)'*inv(TEMP*(TEMP./weight)')*TEMP;
            %obj.HAT=(TEMP./weight)'*inv((TEMP./weight)*(TEMP./weight)')*TEMP./weight

            [U,S,V] = svd(H,0); w=diag(S); winv = 1./w; WW = diag(winv); obj.C = V*WW*U'; %C=inv(H) replaced by SVD;

            obj.p=obj.C*(A'*B);

            obj.p_error=sqrt(diag(obj.C));
             
            obj.chi2=(B-A*obj.p)'*(B-A*obj.p);
            obj.DoF = length(xfit)-obj.M;
            %chi2_dof=obj.chi2/obj.DoF;
        Pvalue=gammainc(obj.chi2/2,obj.DoF/2);
        obj.Q=1-Pvalue;
            %xp=-0.0001:0.00001:0.017;
            %yp=obj.Vpn(obj.p,xp);
            obj.c0=(obj.p(1)-(obj.g2))/log(2.^2);
        end     
        
     function [mj, sj, avg_chi2] = jackknife(obj)
         obj = obj.lin_reg;
         origCl=obj.c0;
         MJiCL=zeros(1,length(obj.OFI));
         MJiC_D=zeros(1,length(obj.OFI));
         for run = 1:length(obj.OFI)
             obj=fi_block(obj,run);
             %obj.FI
             obj = obj.lin_reg;
             MJiCL(run)=obj.c0;
             MJiC_D(run)=obj.chi2/obj.DoF;
             obj.FI=obj.OFI;
         end
         avg_chi2=mean(MJiC_D);
         sj= sqrt((length(obj.OFI)-1).*sum((MJiCL-origCl).^2)/length(obj.OFI));
        mj=origCl-(length(obj.OFI)-1).*(mean(MJiCL,'all')-origCl);
     end
     function ch=bootstrap(obj,smpl_n)
         %[c_val, hist_width]
         c_val_hold=zeros(1,smpl_n);
         for n=1:smpl_n
            pkr=randi([1,7],1,7);
            obj.FI=pkr;
            obj=obj.lin_reg;
            %obj.c0
            c_val_hold(n)=obj.c0;
            
         end
         %obcat(c_val_hold>1)
         %obcat(c_val_hold<0.1)
         %c_val_hold
         %c_val_hold(c_val_hold>1)
         %c_val_hold(c_val_hold<0.1)
         %bins=round(sqrt(length(c_val_hold)));
         bins=22;
         %length(c_val_hold)
         %round(sqrt(length(c_val_hold)))
         histfit(c_val_hold,500);
         ch=fitdist(c_val_hold(:),'normal');
         ch.mean
         ch.sigma
         
     end
     function mat=d_mat(obj)
         %design matrix for polynomials
         mat=[(obj.xL').^0,(obj.xL').^1,(obj.xL').^2];
     end
     
     function RRBS=RRbootstrap(obj,smpl_n)
         obj.lin_reg;
         res=((obj.sigma-obj.g2)/log(2.^2))-((Vpn(obj.p,obj.xL)-obj.g2)/log(2.^2));
         hii=diag(obj.HAT);
         sig_hat=sqrt((1/(length(obj.FI)-obj.M))*(res*res'));
         stand_res=res./((sig_hat)*sqrt(1-hii'));
         op=obj.p;
         osig=obj.sigma;
         osig_err=obj.sigma_err;
         %o_sig_err=obj.sig_err;
         %new_res=zeros(1,length(stand_res));
         ccat=zeros(1,smpl_n);
         %ccaterr=zeros(1,smpl_n);
         
         for n=1:smpl_n
            pkr=randi([1,length(stand_res)],1,length(stand_res)); %change from stand_res
            new_res=res(pkr);
            new_err=osig_err(pkr);
            %obj.sigma_err=o_sig_err(pkr);
            y_0=Vpn(op,obj.xL)+new_res;
            %b=XX_inv*X'*y_0';
            wt=1./osig_err(pkr);
            B=(y_0(pkr).*wt)';
            %X=obj.d_mat.*[wt; wt; wt]';
            X=obj.d_mat;
            X_ip=X'*X;
            [U,S,V] = svd(X_ip,0); w=diag(S); winv = 1./w; WW = diag(winv); XX_inv = V*WW*U'; %C=inv(H) replaced by SVD;
            %b=XX_inv*(X'*B);
            b=XX_inv*X'*y_0';
            
            
            %obj=obj.change_sig(y_0,new_err);
            
            %obj.sigma=y_0;
            %obj.sigma_err=new_err;
            %obj.sigma
            %obj.sigma_err
            %obj.lin_reg;
            %ccat(n)=obj.c0;
            ccat(n)=(b(1)-(obj.g2))/log(2.^2);
            
         end
         histfit(ccat)
         RRBS=fitdist(ccat(:),'normal');
         RRBS.mean
         RRBS.sigma
         
     end
     
    end
    methods(Static)
     function y = Vpn(C,x)
        %simple function that takes in a coefficient array and a smooth 
        %range of x and returns y as a polynomial
        y=zeros(size(x));
        for i=1:length(C)
             y=y + C(i) * x.^(i-1);
        end
     end
     function temp=tempPop(n,m,w)
    %Xhold=zeros(length(m));
    temp=zeros(m,length(w));
    temp(1,:)=w;
    for i=2:m
        Xhold=(n.^(i-1)).*w;
        temp(i,:)=Xhold;
    end
     end
    end
end
