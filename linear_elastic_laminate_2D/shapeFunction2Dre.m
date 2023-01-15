function [NN,dNdedx] = shapeFunction2Dre(N,dNdx,dNde,order)
if(order==1)
    NN=zeros(2,8);
    for i=1:4
        NN(1,2*i-1)=N(i);
        NN(2,2*i)=N(i);
    end
    
    dNdedx=zeros(4,8);
    for i=1:4
        dNdedx(1,2*i-1)=dNdx(i);
        dNdedx(2,2*i-1)=dNde(i);
        dNdedx(3,2*i)=dNdx(i);
        dNdedx(4,2*i)=dNde(i);
    end
end

if(order==2)
    NN=zeros(2,16);
    for i=1:8
        NN(1,2*i-1)=N(i);
        NN(2,2*i)=N(i);
    end
    
    dNdedx=zeros(4,16);
    for i=1:8
        dNdedx(1,2*i-1)=dNdx(i);
        dNdedx(2,2*i-1)=dNde(i);
        dNdedx(3,2*i)=dNdx(i);
        dNdedx(4,2*i)=dNde(i);
    end
end
end