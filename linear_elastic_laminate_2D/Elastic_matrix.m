syms E nu
lambda=E*nu/(1+nu)/(1-nu);
mu=E/2/(1+nu);

delta=eye(2);
delta=sym(delta);

C=zeros(2,2,2,2);
C=sym(C);

%%% Cijkl Elastic moduli
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                C(i,j,k,l)=lambda*delta(i,j)*delta(k,l)+mu*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
            end
        end
    end
end

Ei=[1,1,2,2,1,2];
E=zeros(3,3);
E=sym(E);
for i=1:3
    for j=1:3
        E(i,j)=C(Ei(2*i-1),Ei(2*i),Ei(2*j-1),Ei(2*j));
    end
end
E
C(1,1,1,1)