function P=NormProb(x,mu,sig)
    P=(exp(-(((x-mu)./sig).^2)./(2)))./(sig*sqrt(2*pi));
end