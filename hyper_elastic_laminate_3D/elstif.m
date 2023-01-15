%================= ELEMENT STIFFNESS MATRIX ================================
%
function kel=elstif(nDof,nNoEl,coord,disps,materialprops)
%
%  Assemble the element stiffness
%  Input variables
%nDof                Numbers of degrees of freedom per node (2 in 2D and 3 in 3D)
%nNoEl               Numbers nodes on the element
%coord(a,i)          ith coordinate of ath node
%materialprops       Material properties passed on:constitutive equation
%displacements(a,i)  ith displacement component at ath node
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%det                Determinant of jacobian
%strain(i,j)        strain_ij components
%dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%kel(row,col)       element stiffness matrix, row！！row index, col！！column
%index
%
%
npoints=numberofintegrationpoints(nDof,nNoEl);
dNdx=zeros(nNoEl,nDof);
dxdxi=zeros(nDof,nDof);
strain=zeros(nDof,nDof);
kel=zeros(nDof*nNoEl,nDof*nNoEl);
%
%  Set up integration points && weights    
%
xilist=integrationpoints(nDof,nNoEl,npoints);
w=integrationweights(nDof,nNoEl,npoints);
%
%  Loop over the integration points(intpt)
%
 for intpt=1:1:npoints
     %Compute shape functions && derivatives wrt local coords
     for i=1:1:nDof
         xi(i)=xilist(i,intpt);
     end
     N=shapefunctions(nNoEl,nDof,xi);
     dNdxi=shapefunctionderivs(nNoEl,nDof,xi);
     
     %Compute the jacobian matrix && its determinant
     for i=1:1:nDof
        for j=1:1:nDof
            dxdxi(i,j) = 0.;
            for a=1:1:nNoEl
                dxdxi(i,j)=dxdxi(i,j)+coord(a,i)*dNdxi(a,j);
            end
        end
     end
     
     dxidx=inv(dxdxi);%dζ/dx
     dt=det(dxdxi);
     
     %Convert shape function derivatives:derivatives wrt global coords
     for a=1:1:nNoEl
         for i=1:1:nDof
             dNdx(a,i) = 0.;
             for j=1:1:nDof
                 dNdx(a,i)=dNdx(a,i)+dNdxi(a,j)*dxidx(j,i);
             end
         end
     end
     
     %Compute the deformation gradients by differentiating displacements
     for i=1:1:nDof
         for j=1:1:nDof
            F(i,j)=0.;
            if (i==j)
                F(i,i)=1.;
            end
            for a=1:1:nNoEl
              F(i,j)=F(i,j)+(disps(a,i)*dNdx(a,j));
            end%deformation gradient is formed with displacement and the derivative of shape function wrt global coord
         end
     end
     
     %Compute Bbar and J
     J = det(F);
     B = F*F';%F*F^T
     
     %Convert shape function derivatives to derivatives wrt spatial coords
     Finv=inv(F);
     for a=1:nNoEl
         for i=1:1:nDof
             dNdxs(a,i)=0.;
             for j=1:1:nDof
                 dNdxs(a,i)=dNdxs(a,i)+dNdx(a,j)*Finv(j,i);  
             end
         end
     end
     
     %Compute the Kirchhoff stress
     stress=Kirchhoffstress(nDof,B,J,materialprops);
     
     %Compute the material tangent stiffness (d stress/d strain)
     %C is just C_ijkl for linear elasticity - dsde notation is used
     %to allow extension to nonlinear problems
     C = materialstiffness(nDof,B,J,materialprops);
     
     %Compute the element stiffness
     for a=1:1:nNoEl
         for i=1:1:nDof
             for b=1:1:nNoEl
                 for k=1:1:nDof
                     row=nDof*(a-1)+i;
                     col=nDof*(b-1)+k;
                     for j=1:1:nDof
                         for l=1:1:nDof
                             kel(row,col)=kel(row,col)+C(i,j,k,l)*dNdxs(b,l)*dNdxs(a,j)*w(intpt)*dt;
                         end
                         kel(row,col)=kel(row,col)-stress(i,j)*dNdxs(a,k)*dNdxs(b,j)*w(intpt)*dt;
                     end
                 end
             end
         end
     end
 end
end