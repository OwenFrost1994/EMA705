%
%====================== ELEMENT DISTRIBUTED LOAD VECTOR ==============
%
function r = eldload(nDof,nfacenodes,coords,traction)

  npoints = numberofintegrationpoints(nDof-1,nfacenodes);
  xi = zeros(nDof-1,1);
  dxdxi = zeros(nDof,nDof-1);
  r = zeros(nDof*nfacenodes,1);
   
  xilist = integrationpoints(nDof-1,nfacenodes,npoints);
  w = integrationweights(nDof-1,nfacenodes,npoints);

  for intpt=1:1:npoints

    for i=1:1:nDof-1
      xi(i) = xilist(i,intpt);
    end

    N = shapefunctions(nfacenodes,nDof-1,xi);
    dNdxi = shapefunctionderivs(nfacenodes,nDof-1,xi);
%
%     Compute the jacobian matrix && its determinant
%
    for i=1:1:nDof
      for j=1:1:nDof-1
        dxdxi(i,j) = 0.;
        for a = 1:nfacenodes
          dxdxi(i,j) = dxdxi(i,j) + coords(a,i)*dNdxi(a,j);
        end
      end
    end
    if (nDof == 2) 
      dt = sqrt(dxdxi(1,1)^2+dxdxi(2,1)^2);
    elseif (nDof == 3) 
      dt = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))^2 );
    end
   for a=1:1:nfacenodes
      for i=1:1:nDof
        row = nDof*(a-1)+i;
        r(row) = r(row) + N(a)*traction(i)*w(intpt)*dt;
      end
    end
  end
end
