function [C]=Celasticm(lambda,mu,delta)
C=zeros(2,2,2,2);
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                C(i,j,k,l)=lambda*delta(i,j)*delta(k,l)+mu*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
            end
        end
    end
end