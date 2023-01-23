function g2slide(slid,par,stepped,spouseErrH,spouseErrL)
    %par.Children.DeleteFcn;
    isslide=1;
    indx=round(slid.Value);
    delete(findall(par,'type','line'));
    delete(findall(par,'type','errorbar'));
    set_xl_sig(indx,stepped,spouseErrH,spouseErrL); 
    gridgraph(par,indx,isslide)
end