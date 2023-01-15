function [N, dNdx, dNde]=shapeFunction2D(xi,eta,order)
%xi-¦Î(x), eta-¦Ç(y), dNdx-derivative of N to xi, dNde-derivative of N to eta
%here only the Q4 element(order 1) is temporarily included
%Q4 element only has 4 nodes so only four shape function will appears
if(order==1)%Concepts and Applications of Finite Element Analysis 4ed-(6.2-3)
    
    N(1)=(1-xi)*(1-eta);
    N(2)=(1+xi)*(1-eta);
    N(3)=(1+xi)*(1+eta);
    N(4)=(1-xi)*(1+eta);
    N=N/4;
    
    dNdx(1)=-(1-eta);
    dNdx(2)=(1-eta);
    dNdx(3)=(1+eta);
    dNdx(4)=-(1+eta);
    dNdx=dNdx/4;
    
    dNde(1)=-(1-xi);
    dNde(2)=-(1+xi);
    dNde(3)=(1+xi);
    dNde(4)=(1-xi);
    dNde=dNde/4;
end

%Q8 element has 8 nodes so only four shape function will appears
if(order==2)
    N(5)=(1-xi^2)*(1-eta)/2;
    N(6)=(1+xi)*(1-eta^2)/2;
    N(7)=(1-xi^2)*(1+eta)/2;
    N(8)=(1-xi)*(1-eta^2)/2;
    N(1)=(1-xi)*(1-eta)/4-(N(5)+N(8))/2;
    N(2)=(1+xi)*(1-eta)/4-(N(5)+N(6))/2;
    N(3)=(1+xi)*(1+eta)/4-(N(6)+N(7))/2;
    N(4)=(1-xi)*(1+eta)/4-(N(7)+N(8))/2;
    
    dNdx(5)=-xi*(1-eta);
    dNdx(6)=(1-eta^2)/2;
    dNdx(7)=-xi*(1+eta);
    dNdx(8)=-(1-eta^2)/2;
    dNdx(1)=-(1-eta)/4-(dNdx(5)+dNdx(8))/2;
    dNdx(2)=(1-eta)/4-(dNdx(5)+dNdx(6))/2;
    dNdx(3)=(1+eta)/4-(dNdx(6)+dNdx(7))/2;
    dNdx(4)=-(1+eta)/4-(dNdx(7)+dNdx(8))/2;
    
    dNde(5)=-(1-xi^2)/2;
    dNde(6)=-eta*(1+xi);
    dNde(7)=(1-xi^2)/2;
    dNde(8)=-eta*(1-xi);
    dNde(1)=-(1-xi)/4-(dNde(5)+dNde(8))/2;
    dNde(2)=-(1+xi)/4-(dNde(5)+dNde(6))/2;
    dNde(3)=(1+xi)/4-(dNde(6)+dNde(7))/2;
    dNde(4)=(1-xi)/4-(dNde(7)+dNde(8))/2;
end
end