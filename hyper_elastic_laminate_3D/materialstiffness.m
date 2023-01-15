%================= Material Stiffness ==================================
%
% Computes material elastic tensor C_{ijkl} 
% Currently coded either for plane strain or general 3D.
%
function C=materialstiffness(nDof,B,J,materialprops)

mu1=materialprops(1);%¦Ì
K1=materialprops(2);%¦Ê
dl=[[1,0,0];[0,1,0];[0,0,1]];  

if (nDof == 2)
    C=zeros(2,2,2,2);
    Bqq = B(1,1)+B(2,2)+1.;  
    for i=1:1:2 
        for j=1:1:2 
            for k=1:1:2 
                for l=1:1:2 
                    C(i,j,k,l)=mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k)-(2/3)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l))...
                        +(2/3)*Bqq*dl(i,j)*dl(k,l)/3 )/J^(2/3)+K1*(2*J-1)*J*dl(i,j)*dl(k,l);
                end
            end
        end
    end
elseif (nDof == 3)
    C=zeros(3,3,3,3);
    Bqq = trace(B);
    for i=1:1:3
        for j=1:1:3
            for k=1:1:3
                for l=1:1:3
                    C(i,j,k,l)=mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k)-(2/3)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l)) ...
                         +(2/3)*Bqq*dl(i,j)*dl(k,l)/3 )/J^(2/3)+K1*(2*J-1)*J*dl(i,j)*dl(k,l);
                end
            end
        end
    end
end
end