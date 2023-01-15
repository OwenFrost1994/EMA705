function [kBD,mBD,fBD]=boundarycondition(nDof,nNodes,nDirich,dirich,nNeuman,neuman,kGlobal,mGlobal,fGlobal,H)
%for some nodes, the displacement might be fixed in some direction where
%the k disappear
kD=kGlobal;
mD=mGlobal;
fD=fGlobal;
aDof=nDof*nNodes;
IN1=eye(aDof,aDof);
IN2=-eye(aDof,aDof);
vN=zeros(aDof,1);

for i=1:1:nDirich
    mD(dirich(i,1),:)=0;
    mD(:,dirich(i,1))=0;
    kD(dirich(i,1),:)=0;
    kD(:,dirich(i,1))=0;
    kD(dirich(i,1),dirich(i,1))=1;
    fD(dirich(i,1))=dirich(i,2);
end
for i=1:1:aDof
    if ismember(i,neuman(:,1))
        
    else
        for j =1:1:nDirich
            fD(i)=fD(i)-kGlobal(i,dirich(j,1))*dirich(j,2);
        end
    end
end
for i=1:1:nNeuman
    IN2(neuman(i,1),neuman(i,1))=0;
    vN(neuman(i,1))=neuman(i,2)*H;
end
mD(:,dirich(:,1))=[];
mD(dirich(:,1),:)=[];
kD(:,dirich(:,1))=[];
kD(dirich(:,1),:)=[];
fD(dirich(:,1))=[];
IN1(:,neuman(:,1))=[];
IN1(neuman(:,1),:)=[];
IN2(:,neuman(:,1))=[];
IN2(neuman(:,1),:)=[];
vN(neuman(:,1))=[];

kBD=[zeros(aDof-nNeuman,aDof-nDirich),IN2;kD,zeros(aDof-nDirich,aDof-nNeuman)];
mBD=[IN1,zeros(aDof-nNeuman,aDof-nDirich);zeros(aDof-nDirich,aDof-nNeuman),mD];
fBD=[vN;fD];
end