function [M]=nodalMass(N,rho,delta,A,B)
% this is nodal mass matrix, not the element mass matrix
% in 3D M is 3-3, in 2D M is 2-2
% rho is density
% delta is identical matrix
M=zeros(2,2);
for i=1:2
    for k=1:2
        M(i,k)=N(A)*rho*delta(i,k)*N(B);
    end
end
end