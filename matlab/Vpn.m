function y = Vpn(C,x)
%simple function that takes in a coefficient array and a smooth 
%range of x and returns y as a polynomial
    y=zeros(size(x));
    for i=1:length(C)
       y=y + C(i) * x.^(i-1);
    end
end