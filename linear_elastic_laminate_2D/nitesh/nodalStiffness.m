function[K]=nodalStiffness(dNdX,C,A,B)
%nodal stiffness matrix
K=zeros(2,2);
dNdX_t=dNdX';

for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                K(i,k)= K(i,k)+ dNdX_t(A,l)*C(i,j,k,l)*dNdX(j,B);
            end
        end
    end
end
end