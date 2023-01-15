function [kLocal,mLocal,fLocal]=localmatrix(nDof,nNoEl,wq,gaussMatrix,order,connArray,coorNo,E,nu,rho,e,t)
%generate following matrices
kLocal=zeros(nDof*nNoEl,nDof*nNoEl);%local stiffness matrix
mLocal=zeros(nDof*nNoEl,nDof*nNoEl);%local mass matrix
fLocal=zeros(nDof*nNoEl,1);%local load matrix
Em=(E/(1-nu^2))*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
T=[1,0,0,0;0,0,0,1;0,1,1,0];
for wqX=1:wq
    for wqY=1:wq%Concepts and Applications of Finite Element Analysis 4ed-(6.2-6)
        [N,dNdx,dNde] = shapeFunction2D(gaussMatrix(wqX,1),gaussMatrix(wqY,1),order);
        [NN,dNdedx] = shapeFunction2Dre(N,dNdx,dNde,order);
        Jac=[dNdx*coorNo(connArray(e,:),1), dNdx*coorNo(connArray(e,:),2); dNde*coorNo(connArray(e,:),1), dNde*coorNo(connArray(e,:),2)];%Jacobian matrix
        B=T*[inv(Jac),zeros(size(Jac,1),size(Jac,1));zeros(size(Jac,1),size(Jac,1)),inv(Jac)]*dNdedx;%B matrix
        kLocal=kLocal+B'*E*B*det(Jac)*t*gaussMatrix(wqX,2)*gaussMatrix(wqY,2);
        mLocal=mLocal+NN'*NN*rho*det(Jac)*t*gaussMatrix(wqX,2)*gaussMatrix(wqY,2);
        fLocal=fLocal+0;    %no force term in the question
    end
end