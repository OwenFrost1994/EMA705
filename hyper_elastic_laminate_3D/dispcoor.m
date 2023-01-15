function [wv] = dispcoor(w,nDof,nNodes)
%
%Transform 1D displacement works with Stiffness matrix into same form
%node coordinate storage: coor[a,i] a node, ith coordinate
%
wv=zeros(nNodes,2);
if nDof == 2
    for i=1:1:nNodes
        wv(i,1)=w(nDof*(i-1)+1);
        wv(i,2)=w(nDof*(i-1)+2);
    end
end
if nDof == 3
    for i=1:1:nNodes
        wv(i,1)=w(nDof*(i-1)+1);
        wv(i,2)=w(nDof*(i-1)+2);
        wv(i,3)=w(nDof*(i-1)+3);
    end
end
end