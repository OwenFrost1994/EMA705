%
%================= ELEMENT RESIDUAL VECTOR ================================
%
function rel = elresid(nDof,nNoEl,coord,materialprops,disp)
%
%  Assemble the element residual force
%
% Inputs:
%nDof                node degree of freedom 2-2D,3-3D
%nNoEl               number of nodes in an elements
%coord               coordinates of ith node in previous status
%materialprops       Material properties passed on:constitutive equation
%displacements(a,i)  ith displacement component at ath node

%
% Local variables
%npoints            No. integration points
%xi[i,inpt]         ith local coord of integration point no. intpt
%w[intpt]           weight for integration point no. intpt
%N[a]               Shape function associated with ath node on element
%dNdxi[a,i]         Derivative of ath shape function wrt ith local coord
%dNdx[a,i]          Derivative of ath shape function wrt ith global coord
%dxdxi[i,j]         Derivative of ith global coord wrt jth local coord
%dxidx[i,j]         Derivative of ith local coord wrt jth global coord
%det                Determinant of jacobian
%strain[i,j]        strain_ij components
%stress[i,j]        stress_ij components
%r[row]             Residual vector
%
   npoints = numberofintegrationpoints(nDof,nNoEl);
   dxdxi = zeros(nDof,nDof);
   dxidx = zeros(nDof,nDof);
   dNdxs = zeros(nNoEl,nDof);
   rel = zeros(nDof*nNoEl,1);
%
%  Set up integration points and weights    
%
   xilist = integrationpoints(nDof,nNoEl,npoints);
   w = integrationweights(nDof,nNoEl,npoints);
   %
%  Loop over the integration points
%
   for intpt=1:1:npoints

%     Compute shape functions and derivatives wrt local coords
%
      for i=1:1:nDof
        xi(i) = xilist(i,intpt);
      end      
      N = shapefunctions(nNoEl,nDof,xi);
      dNdxi = shapefunctionderivs(nNoEl,nDof,xi);
%
%     Compute the jacobian matrix and its determinant
%
      for i=1:1:nDof
        for j=1:1:nDof
          dxdxi(i,j) = 0.;
          for a = 1:1:nNoEl
            dxdxi(i,j) = dxdxi(i,j) + coord(a,i)*dNdxi(a,j);
          end
        end
      end

      dxidx = inv(dxdxi);
      dt = det(dxdxi);
%
%     Convert shape function derivatives to derivatives wrt global coords
%
      for a=1:1:nNoEl
        for i=1:1:nDof
          dNdx(a,i) = 0.;
          for j=1:1:nDof
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end
%
%     Compute the deformation gradients by differentiating displacements
%
      for i=1:1:nDof
         for j=1:1:nDof
            F(i,j) = 0.;
            if i==j
                F(i,i) = 1.;
            end
            for a=1:1:nNoEl
              F(i,j) = F(i,j) + (disp(a,i)*dNdx(a,j));
            end
         end
      end
%
%     Compute Bbar and J
%
      J = det(F);
      B = F*transpose(F);
%
%     Convert shape function derivatives to derivatives wrt spatial coords
%
      Finv = inv(F);
      for a=1:1:nNoEl
        for i=1:1:nDof
          dNdxs(a,i) = 0.;
          for j=1:1:nDof
          dNdxs(a,i) = dNdxs(a,i) + dNdx(a,j)*Finv(j,i);  
          end
        end
      end
%
%     Compute the stress
%
      stress = Kirchhoffstress(nDof,B,J,materialprops);
%
%     Compute the element residual
%             
      for a=1:1:nNoEl
        for i=1:1:nDof
          row = nDof*(a-1)+i;
          for j=1:1:nDof
            rel(row) = rel(row) + stress(i,j)*dNdxs(a,j)*w(intpt)*dt;
          end
        end
      end
   end

end