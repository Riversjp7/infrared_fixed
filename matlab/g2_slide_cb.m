function g2_slide_cb(slid,par,BOB)
    indx=round(slid.Value);
    BOB=indx_change(BOB,indx);
    BOB=lin_reg(BOB);
    fit_num=4;
    DN=["All","-8/10","5L","4L"];
    
    xpa=-0.0001:0.00001:0.017;
    yp1=Vpn(BOB.p,xpa);
    
    c0_ray=zeros(1,fit_num);
    ypa=zeros(fit_num,length(yp1));
    perray=zeros(fit_num,length(BOB.p_error));
    
    c0_ray(1)=BOB.c0;
    ypa(1,1:end)=yp1;
    perray(1,1:end)=BOB.p_error;
    
    ob_hold=BOB;
    for i=1:3
        ob_hold=fi_block(ob_hold,length(BOB.xL)+1-i);
        switch i
            case 1
                ob_hold.M=3;
            otherwise
                ob_hold.M=2;
        end
        ob_hold=lin_reg(ob_hold);
        c0_ray(i+1)=ob_hold.c0;
        perray(i+1,1:end)=ob_hold.p_error(1);
        ypa(i+1,1:end)=Vpn(ob_hold.p,xpa);
    end
    
    y_dat1=[ypa(1,1:end);ypa(1,1:end);ypa(1,1:end);ypa(2,1:end);ypa(2,1:end);ypa(3,1:end)];
    y_dat2=[ypa(2,1:end);ypa(3,1:end);ypa(4,1:end);ypa(3,1:end);ypa(4,1:end);ypa(4,1:end)];
    cplots1=[c0_ray(1),c0_ray(1),c0_ray(1),c0_ray(2),c0_ray(2),c0_ray(3)];
    cplots2=[c0_ray(2),c0_ray(3),c0_ray(4),c0_ray(3),c0_ray(4),c0_ray(4)];
    cerrs1=[perray(1,1),perray(1,1),perray(1,1),perray(2,1),perray(2,1),perray(3,1)]/log(2.^2);
    cerrs2=[perray(2,1),perray(3,1),perray(4,1),perray(3,1),perray(4,1),perray(4,1)]/log(2.^2);
    DN1=[DN(1),DN(1),DN(1),DN(2),DN(2),DN(3)];
    DN2=[DN(2),DN(3),DN(4),DN(3),DN(4),DN(4)];
    set(par,'ColorOrderIndex',1); 
    lines1=findobj(par,'Type','line','Tag','l1');
    lines2=findobj(par,'Type','line','Tag','l2');
    eb1=findobj(par,'Type','errorbar','Tag','l1');
    eb2=findobj(par,'Type','errorbar','Tag','l2');
    
    %size(lines1(1).XData)
    %size(y_dat1(1,1:end))
    %size(ypa)
    
    for i=1:6
    set(lines1(i),'ydata',(y_dat1(i,1:end)-BOB.g2)/log(2.^2),...
             'DisplayName',DN1(i));
    set(lines2(i),'ydata',(y_dat2(i,1:end)-BOB.g2)/log(2.^2),...
            'DisplayName',DN2(i));
    set(eb1(i),'ydata',cplots1(i),...
               'YPositiveDelta',cerrs1(i),...
               'YNegativeDelta',cerrs1(i),...
               'DisplayName',DN1(i),...
               'Marker','o');
    set(eb2(i),'ydata',cplots2(i),...
                'YPositiveDelta',cerrs2(i),...
                'YNegativeDelta',cerrs2(i),...
                'DisplayName',DN2(i),...
                'Marker','o');
    legend(par(i))
    end
    set(findobj(par,'Type','errorbar','Tag','perr'),'ydata',(BOB.sigma-BOB.g2)/log(2.^2),...
                'YPositiveDelta',BOB.sigma_err/log(2.^2),...
                'YNegativeDelta',BOB.sigma_err/log(2.^2),...
                'Marker', 'o',...
                'Color','b');
    
    
    %{
    a rather complicated relationship between grid positions 1-6
    and the algorithm that compares the first beta object with the
    consecutive beta objects until the end of the array is reached
    and then moves on to the second beta object and compares it to
    the end of the ray and so forth.
  
     
    axi=[-1,1,2];
    for pair1=1:3
        
         pair2=pair1+1;
         while pair2<=4
            
            axpos=pair2+axi(pair1);
            delete(findall(par(axpos).Children))
            set(par(axpos),'ColorOrderIndex',1);
            
           
            errorbar(par(axpos),0.0000,b_obj(pair1).c0,b_obj(pair1).p_error(1)/log(2.^2),...
                'LineWidth',1,...
                'DisplayName',DN(pair1));
            plot(par(axpos),xpa,(ypa(pair1,1:end)-b_obj(pair1).g2)/log(2.^2),...
                'LineWidth',1,...
                'DisplayName',DN(pair1));
            
            errorbar(par(axpos),0.0000,b_obj(pair2).c0,b_obj(pair2).p_error(1)/log(2.^2),...
                'LineWidth',1,....
                'DisplayName',DN(pair2));
            plot(par(axpos),xpa,(ypa(pair2,1:end)-b_obj(pair2).g2)/log(2.^2),...
                'LineWidth',1,...
                'DisplayName',DN(pair2));
       
            leg=legend(par(axpos));
            errorbar(par(axpos),b_obj(1).xL,(b_obj(1).sigma-b_obj(1).g2)/log(2.^2),b_obj(1).sigma_err/log(2.^2),'bo','LineWidth',1);
            pair2=pair2 +1;
         end
    end
      %}
end