%======================= Uniaxial tangent ===================
%
function dsde = hardeningslope(ee,materialprops)

s0 = materialprops(1);
e0 = materialprops(2);
n = materialprops(3);

if (ee < e0)
    if (n-1<10^(-12))
        dsde = s0/e0;
    else
        dsde = s0*(n/(n-1)-ee/e0)/e0;
        dsde = dsde/sqrt( (1+n^2)/(n-1)^2 - (n/(n-1)-ee/e0)^2 );
    end
else
    dsde = s0*( (ee/e0)^(1/n)  )/(n*ee);
end

end