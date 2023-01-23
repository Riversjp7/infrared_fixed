

fig=uifigure('Name','Continuum limit','Position',[0,0,1375,800]);

BOB=Clim(stepped,spouseErrH,spouseErrL);
BOB.M=3;

tele_grid=uigridlayout(fig,[1,2]);
tele_grid.ColumnWidth={100,'1x'};
pan1=uipanel(tele_grid);
c0_button=uibutton(pan1,'push',...
                    'Text','Compare C0',...
                    'Position',[5,650,90,25],...
                    'ButtonPushedFcn',@(btn,event) c0_window(btn,BOB));
pan2=uipanel(tele_grid);

sld_grid=uigridlayout(pan2,[2,1]);
sld_grid.RowHeight = {'fit', 'fit'};
sld_grid.ColumnWidth={'fit'};

grd_pan=uipanel(sld_grid);
sld_pan=uipanel(sld_grid);
sld_pan_dim=get(sld_pan,'Position');
init_sld=40;

graph_grid = uigridlayout(grd_pan,[2,3],...
                            'RowSpacing',0,...
                            'ColumnSpacing',0,...
                            'Padding',[1,1,1,1]);
graph_grid.RowHeight={'1x','1x'};
graph_grid.ColumnWidth={'1x','1x','1x'};
%{
sldr=uicontrol(sld_pan,...
    'Style','Slider',...
    'Min',0,'Max',61,...
    'SliderStep',[1 10],...
    'Position', [25,75,.9*sld_pan_dim(3),25],...
    'Value',init_sld,...
    'Callback', {@g2slide,graph_grid,stepped,spouseErrH,spouseErrL});
%}


            %graph_grid.RowSpacing=0
%graph_grid.ColumnSpacing=0
%graph_grid.Padding=
set_xl_sig(indx,stepped,spouseErrH,spouseErrL);
ax_tags=["All vs All - 8/10", "All vs 5L", "All vs 4L", "-8/10 vs 5L", "-8/10 vs 4L", "5L vs 4L"]
for pos=1:6
    ax(pos)=uiaxes(graph_grid,'xlim',[-0.001  0.017],'ylim',[0.2 1.1],...
        'Tag',ax_tags(pos));
    if pos<=3
    %ax=uiaxes(par,'xlim',[-0.001  0.017],'ylim',[0.2 1.1]);
    ax(pos).Layout.Row= 1;
    ax(pos).Layout.Column=pos;
    elseif pos>3
    %ax=uiaxes(par,'xlim',[-0.001  0.017],'ylim',[0.2 1.1]);
    ax(pos).Layout.Row=2;
    ax(pos).Layout.Column=pos-3;
    end
    ax(pos).XGrid='on';
    ax(pos).YGrid='on';
    ax(pos).Box='on';
    ax(pos).NextPlot='add';
    

end

  gridgraph(ax,init_sld,0);

itog_convert=4.+ 0.1*([1:61]-1);
sldr= uislider(sld_pan,...
                'Position', [25,75,.9*sld_pan_dim(3),25],...
                'Value',init_sld,...   
                'ValueChangedFcn',@(sldr,event) g2slide(sldr,ax,stepped,spouseErrH,spouseErrL),...
                'Limits', [1,61],...
                'MajorTicks',1:61,...
                'MajorTickLabels', categorical(string(itog_convert)));

