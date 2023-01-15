%================= Material Stress ==================================
%
%   Computes stress sigma_{ij} given strain epsilon_{ij}
%
function stress = materialstress(ndof,ncoord,strain,materialprops)

%  Bulk modulus
K = materialprops(3)*(materialprops(1)/materialprops(2))/ ...
    (3*(1-2*materialprops(4)));
stress = zeros(ndof,ncoord);
e = zeros(ndof,ncoord);
dl = [ [1,0,0];[0,1,0];[0,0,1] ];  

evol = strain(1,1) + strain(2,2);
if (ndof==3) evol = evol + strain(3,3); end

ee = 0;
for i = 1 : ndof
    for j = 1 : ncoord
        e(i,j) = strain(i,j) - dl(i,j)*evol/3;
        ee = ee + e(i,j)*e(i,j);
    end
end
ee = sqrt(2.*ee/3.);
se = effectivestress(ee,materialprops);

for i = 1 : ndof
    for j = 1 : ncoord
        stress(i,j) = K*evol*dl(i,j)/3;
        if (ee>0) stress(i,j) = stress(i,j) + (2/3)*se*e(i,j)/ee; end
    end
end

end