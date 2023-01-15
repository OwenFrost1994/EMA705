function [XE] = indexT(X,nDof,nNodes)
XE=zeros(nNodes,2);
for i=1:1:nNodes
    XE(i,1)=X(nDof*(i-1)+1);
    XE(i,2)=X(nDof*(i-1)+2);
end