 %====================== No. nodes on element faces ================
%
%   This procedure returns the number of nodes on each element face
%   for various element types.  This info is needed for computing
%   the surface integrals associated with the element traction vector
%
function n = nfacenodes(nDof,nNoEl)
   if nDof == 2
     if nNoEl == 4
         n = 2;
     elseif nNoEl == 8
         n = 3;
     end
   elseif nDof == 3
     if nNoEl == 8
         n = 4;
     elseif nNoEl == 20
         n = 8;
     end
   end
end
