function [kD,mD,fD]=boundarycondition(nDof,nNodes,nDirich,dirich,nNeuman,neuman,kGlobal,mGlobal,fGlobal,D)
%for some nodes, the displacement might be fixed in some direction where
%the k disappear
kD=kGlobal;
mD=mGlobal;
fD=fGlobal;
aDof=nDof*nNodes;

for i=1:1:nDirich
    mD(dirich(i,1),:)=0;
    mD(:,dirich(i,1))=0;
    kD(dirich(i,1),:)=0;
    kD(:,dirich(i,1))=0;
    kD(dirich(i,1),dirich(i,1))=1;
    fD(dirich(i,1))=dirich(i,2);
end
for i=1:1:nNeuman
    mD(neuman(i,1),:)=0;
    mD(:,neuman(i,1))=0;
    kD(neuman(i,1),:)=0;
    kD(:,neuman(i,1))=0;
    kD(neuman(i,1),neuman(i,1))=1;
    fD(neuman(i,1))=neuman(i,2)*D;
end
for i=1:1:aDof
    if ismember(i,dirich(:,1))
        
    else
        for j =1:1:nDirich
            fD(i)=fD(i)-kGlobal(i,dirich(j,1))*dirich(j,2);
        end
    end
    if ismember(i,neuman(:,1))
        
    else
        for j =1:1:nNeuman
            fD(i)=fD(i)-kGlobal(i,neuman(j,1))*neuman(j,2)*D;
        end
    end
end

for i=1:1:aDof
    
end
end