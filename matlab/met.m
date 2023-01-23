y_len=10000;
step=0.5;
rej_count=0;
isreject=0;

x=rand;
%y=x;
xtest=linspace(-2,2,10000);
%xProp=zeros(1,y_len);
yacc=zeros(1,y_len);
%ydec=zeros(1,y_len);

sig=0.35;
mu=0.5;

%P=NormProb(xtest,mu,sig);
%plot(xtest,P,'-b');


yacc(1)=x;

%xProp=yacc(1)+step*(2*rand-1);
for i=2:y_len
    isreject=1;
xProp=yacc(i-1)+step*(2*rand-1);
while isreject==1
rat=NormProb(xProp,mu,sig)/NormProb(yacc(i-1),mu,sig);
if rat>1
   yacc(i)=xProp;
   isreject=0;
else
    z=rand;
    if z<rat
        yacc(i)=xProp;
        isreject=0;
    elseif z>rat
        rej_count= rej_count+1;
        isreject=1;
        xProp=yacc(i-1)+step*(2*rand-1);
        
    end
end   
end
end
%plot(linspace(0,1,length(yacc)),yacc,'.k')
histfit(yacc,round(sqrt(y_len)));
p=fitdist(yacc(:),'normal');
p.mean
p.sigma
mean(yacc)
