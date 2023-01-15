%================= Stress ==================================
%
%   Computes stress Kirchhoff_{ij} given B_{ij}
%
function stress=Kirchhoffstress(nDof,B,J,materialprops)

stress=zeros(nDof,nDof);
dl=[[1,0,0];[0,1,0];[0,0,1]];%¦Äij, identical matrix

mu1=materialprops(1);%¦Ì
K1=materialprops(2);%¦Ê

Bkk=trace(B);
if nDof==2
    Bkk=Bkk+1;
end 
for i=1:1:nDof
    for j =1:1:nDof
        stress(i,j)=mu1*(B(i,j)-Bkk*dl(i,j)/3.)/J^(2/3)+K1*J*(J-1)*dl(i,j);
    end
end
end
