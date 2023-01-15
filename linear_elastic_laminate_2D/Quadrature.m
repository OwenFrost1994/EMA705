function[A]=Quadrature(n)
%gauss quadrature, points and weights
if(n==1)
    A(1,1:2)=[0 2];
end

if(n==2)
    A(1,1:2)=[-1/sqrt(3) 1];
    A(2,1:2)=[1/sqrt(3) 1];
end

if(n==3)
    A(1,1:2)=[-sqrt(3/5) 5/9];
    A(2,1:2)=[0 8/9];
    A(3,1:2)=[sqrt(3/5) 5/9];
end

end