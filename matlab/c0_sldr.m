function c0_sldr(event, ax,BOB)

indx=round(event.Value);
BOB=indx_change(BOB,indx);
comp_num=6;



delete(findall(ax.Children));
delete(findall(ax.Parent,'Tag','ax3'));
delete(findall(ax.Parent,'Tag','ax2'));

sldpan=findobj(ax.Parent.Parent.Children(1));
sldpan(1).Title=['g' char(0178) '= ' num2str(BOB.g2)];
c0_cat=zeros(1,comp_num);
err_cat=zeros(1,comp_num);
chi_cat=zeros(1,comp_num);

BOB=lin_reg(BOB);
%bsfit=bootstrap(BOB,500);
bsfit=RRbootstrap(BOB,500);
mb=bsfit.mean;
sb=bsfit.sigma;
c0_cat(1)=mb;
err_cat(1)=sb;
chi_cat(1)=0;

[mj,sj,avgC_D]=BOB.jackknife;
c0_cat(2)=mj;
err_cat(2)=sj;
chi_cat(2)=avgC_D;
c0_cat(3)=BOB.c0;
err_cat(3)=BOB.p_error(1);
chi_cat(3)=BOB.chi2/BOB.DoF;
BOB=fi_block(BOB, 7);
BOB=lin_reg(BOB);
c0_cat(4)=BOB.c0;
err_cat(4)=BOB.p_error(1);
chi_cat(4)=BOB.chi2/BOB.DoF;

BOB.M=2;
BOB=fi_block(BOB, 6);
BOB=lin_reg(BOB);
c0_cat(5)=BOB.c0;
err_cat(5)=BOB.p_error(1);
chi_cat(5)=BOB.chi2/BOB.DoF;

BOB=fi_block(BOB, 5);
BOB=lin_reg(BOB);
c0_cat(6)=BOB.c0;
err_cat(6)=BOB.p_error(1);
chi_cat(6)=BOB.chi2/BOB.DoF;


label_ray=["Bootstrap" "Jackknife" "Original" "8/10 removed" "5 linear" "4 linear"];
set(ax,'ColorOrderIndex',1);

color_ray=["#636363" "black", "magenta", "blue","green","red"];
eplots=zeros(1,comp_num);
eplots2=zeros(1,comp_num);
eplots3=zeros(1,comp_num);
for i=1:comp_num
eplots(i)=errorbar(ax,0.1*i,c0_cat(i),err_cat(i)./log(2.^2),'linestyle','none',...
            'Marker','o',...
            'DisplayName',label_ray(i),...
            'Color',color_ray(i));
eplots2(i)=errorbar(ax,0.1*i,c0_cat(i),err_cat(i)./log(2.^2),'linestyle','none',...
            'Marker','o',...
            'DisplayName',label_ray(i),...
            'Color',color_ray(i));
eplots3(i)=errorbar(ax,0.1*i,c0_cat(i),err_cat(i)./log(2.^2),'linestyle','none',...
            'Marker','o',...
            'DisplayName',label_ray(i),...
            'Color',color_ray(i));
end
%text(ax,linspace(0.1,0.5,5),c0_cat-0.1,label_ray);
leg1=legend(ax,eplots,categorical(string(c0_cat)),'Location','southeast');
leg2_ax=uiaxes(ax.Parent,'position',get(ax,'position'),...
                        'visible','off',...
                        'tag','ax2');
leg2=legend(leg2_ax,eplots2,categorical(string(chi_cat)),'Location','East');
title(leg1,'C0');
title(leg2,texlabel('chi^2/DoF'));

leg3_ax=uiaxes(ax.Parent,'position',get(ax,'position'),...
                        'visible','off',...
                        'tag','ax3');

leg3=legend(leg3_ax,eplots3,label_ray,'Location','northeast');

%findall(ax.Parent)
%delete(leg2_ax);
%findobj(ax,'type','errorbar')
%errorbar(ax,0.5,0.1*indx,0.2);
%errorbar(ax,0.4,BOB.c0, 0.2);

end