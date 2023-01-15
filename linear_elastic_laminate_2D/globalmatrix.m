function [kGlobal,mGlobal,fGlobal]=globalmatrix(nDof,nNodes,nNoEl,nElements,ElementsM,ElementsL,wq,order,connArray,coorNo,E1,E2,nu,rho,t)
gaussMatrix=Quadrature(wq); %gaussian  quadrature matirx

kGlobal=zeros(nDof*nNodes,nDof*nNodes);%empty global stiffness matrix

mGlobal=zeros(nDof*nNodes,nDof*nNodes);%empty global mass matrix

fGlobal=zeros(nDof*nNodes,1);%empty global load matrix

for e=1:nElements%loop of every element
    if ismember(e,ElementsM)
        [kLocal,mLocal,fLocal]=localmatrix(nDof,nNoEl,wq,gaussMatrix,order,connArray,coorNo,E1,nu,rho,e,t);
    end
    if ismember(e,ElementsL)
        [kLocal,mLocal,fLocal]=localmatrix(nDof,nNoEl,wq,gaussMatrix,order,connArray,coorNo,E2,nu,rho,e,t);
    end
    %Assembly of local matrics into global matrics
    for i=1:nNoEl
        for j=1:nNoEl
            I=2*(connArray(e,i)-1)+1:2*connArray(e,i);
            J=2*(connArray(e,j)-1)+1:2*connArray(e,j);
            kGlobal(I,J)=kGlobal(I,J)+kLocal(2*(i-1)+1:2*i,2*(j-1)+1:2*j);
            mGlobal(I,J)=mGlobal(I,J)+mLocal(2*(i-1)+1:2*i,2*(j-1)+1:2*j);
        end
%       I=2*(connArray(e,i)-1)+1:2*connArray(e,i);
        fGlobal(I,1)=fGlobal(I,1)+fLocal(2*(i-1)+1:2*i);
    end
end
end