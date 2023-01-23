function temp=tempPop(n,m,w)
    %Xhold=zeros(length(m));
    temp=zeros(m,length(w));
    temp(1,:)=w;
    for i=2:m
        Xhold=(n.^(i-1)).*w;
        temp(i,:)=Xhold;
    end