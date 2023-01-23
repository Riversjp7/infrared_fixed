function g2_slide(slid,par,BOB)
    indx=round(slid.Value);
    BOB=indx_change(BOB,indx);
    BOB=lin_reg(BOB);
    xpa=-0.0001:0.00001:0.017;

    b_obj(1)=BOB;
    l=Vpn(BOB.p,xpa);
    ypa=zeros(4,length(l));

   ypa(1, 1:end)=Vpn(BOB.p,xpa);
  
    DN=["All","-8/10","5L","4L"];
    
    
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
        b_obj(i+1)=ob_hold;
        ypa(i+1,1:end)=Vpn(ob_hold.p,xpa);
    end
    %{
    a rather complicated relationship between grid positions 1-6
    and the algorithm that compares the first beta object with the
    consecutive beta objects until the end of the array is reached
    and then moves on to the second beta object and compares it to
    the end of the ray and so forth.
    %}
    axi=[-1,1,2];
    for pair1=1:3
        
         pair2=pair1+1;
         while pair2<=4
            
            axpos=pair2+axi(pair1);
            delete(findall(par(axpos).Children))
            set(par(axpos),'ColorOrderIndex',1);
            
            switch axpos
                case 1
                    set(par(axpos),'xticklabel',[]); 
                case 4
                    
                case 5
                    set(par(axpos),'yticklabel',[]);
                case 6
                    set(par(axpos),'yticklabel',[]);
                otherwise
                    set(par(axpos),'xticklabel', [],...
                                'yticklabel', []);
            end
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
end