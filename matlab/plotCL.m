function plotCL(ax,p,p_error,g2);
    
    xp=-0.0001:0.00001:0.017;
    yp=Vpn(p,xp);
    errorbar(ax,0.0000,(p(1)-g2)/log(2.^2),p_error(1)/log(2.^2),'LineWidth',1);
    plot(ax,xp,(yp-g2)/log(2.^2),'LineWidth',1,'Tag','cp');
    
end