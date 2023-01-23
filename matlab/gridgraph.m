function g=gridgraph(par,indx,isslide)
    
    %indx %needs to be replaced by slider here and in Clim_compare
    g2=4.+ 0.1*(indx-1); 
    [xL,sigma,sigma_err]=get_xl_sig;
     
for pos=1:6
   %{ 
if isslide==1   
    switch pos
        case 1
            delete(findall(par,'Tag','All vs. All-1'));
        case 2
            delete(findall(par,'Tag','All vs 5L'));
        case 3
            delete(findall(par,'Tag','All vs 4L'));
        case 4
            delete(findall(par,'Tag','All-1 vs 5L'));
        case 5
            delete(findall(par,'Tag','All-1 vs 4L'));
        case 6
            delete(findall(par,'Tag','5L vs 4L'));
    end
end
    
    par(pos)=uipar(pos)es(par,'xlim',[-0.001  0.017],'ylim',[0.2 1.1]);
    if pos<=3
    %par(pos)=uipar(pos)es(par,'xlim',[-0.001  0.017],'ylim',[0.2 1.1]);
    par(pos).Layout.Row= 1;
    par(pos).Layout.Column=pos;
    elseif pos>3
    %par(pos)=uipar(pos)es(par,'xlim',[-0.001  0.017],'ylim',[0.2 1.1]);
    par(pos).Layout.Row=2;
    par(pos).Layout.Column=pos-3;
    end
    par(pos).XGrid='on';
    par(pos).YGrid='on';
    par(pos).Box='on';
   
    par(pos).NextPlot='add';
    %cla reset%par(pos).NextPlot='replace';
    %}
    set(par(pos),'ColorOrderIndex',1);
    switch pos
        case 1
            set(par(pos),'xticklabel',[],'Tag','All vs. All-1');
            [p,p_error,chi2]=CLim_compare(7,3,indx);
            plotCL(par(pos),p,p_error,g2);
            set(findobj(par(pos),'Tag','cp'),'Tag','all','DisplayName','All');
            clear p p_error chi2;
            [p,p_error,chi2]=CLim_compare(6,3,indx);
            plotCL(par(pos),p,p_error,g2);
            findobj(gca)
            set(findobj(par(pos),'Tag','cp'),'Tag','prmv','DisplayName','no 8/10');
        case 2
             set(par(pos),'xticklabel', [],...
                'yticklabel', [], ...
                'Tag', 'All vs 5L');
            [p,p_error,chi2]=CLim_compare(7,3,indx);
            plotCL(par(pos),p,p_error,g2);
            %par(pos).NextPlot='add';
            clear p p_error chi2;
            [p,p_error,chi2]=CLim_compare(5,2,indx);
            plotCL(par(pos),p,p_error,g2);
        case 3
             set(par(pos),'xticklabel', [],...
                'yticklabel', [],...
                'Tag','All vs 4L');
            [p,p_error,chi2]=CLim_compare(7,3,indx);
            plotCL(par(pos),p,p_error,g2);
            %par(pos).NextPlot='add';
            clear p p_error chi2;
            [p,p_error,chi2]=CLim_compare(4,2,indx);
            plotCL(par(pos),p,p_error,g2);
        case 4
            set(par(pos),'Tag','All-1 vs 5L');
           [p,p_error,chi2]=CLim_compare(6,3,indx);
           plotCL(par(pos),p,p_error,g2);
           %par(pos).NextPlot='add';
           clear p p_error chi2;
           [p,p_error,chi2]=CLim_compare(5,2,indx);
           plotCL(par(pos),p,p_error,g2);
        case 5
            set(par(pos),'yticklabel',[],...
                'Tag','All-1 vs 4L');
            [p,p_error,chi2]=CLim_compare(6,3,indx);
            plotCL(par(pos),p,p_error,g2);
            %par(pos).NextPlot='add';
            clear p p_error chi2;
            [p,p_error,chi2]=CLim_compare(4,2,indx);
            plotCL(par(pos),p,p_error,g2);
        case 6
            set(par(pos),'yticklabel',[],...
                'Tag','5L vs 4L');
           [p,p_error,chi2]= CLim_compare(5,2,indx);
           plotCL(par(pos),p,p_error,g2);
          % par(pos).NextPlot='add';
           clear p p_error chi2;
           [p,p_error,chi2]= CLim_compare(4,2,indx);
           plotCL(par(pos),p,p_error,g2);
        %otherwise
               % set(par(pos),'xticklabel', [],...
                %'yticklabel', []);
    end  
    leg=legend(par(pos));
    errorbar(par(pos),xL,(sigma-g2)/log(2.^2),sigma_err/log(2.^2),'bo','LineWidth',1);
end
%par(pos).NextPlot='replaceall';
end
            