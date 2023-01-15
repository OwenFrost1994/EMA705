%======================= Lists of nodes on element faces =============
%
%    This procedure returns the list of nodes on an element's face
%    The nodes are ordered so that the element face forms
%    a line for 2D problems or 2D surface 3D problems
%
function list = facenodes(nDof,nNoEl,face)

   i4 = [2,3,4,1]; 

   list = zeros(nfacenodes(nDof,nNoEl),1);

   if nDof == 2
     if nNoEl==4
       list(1) = face;
       list(2) = i4(face);
     elseif nNoEl==8
       list(1) = face;
       list(2) = i4(face);
       list(3) = face+4;
     end
   elseif nDof == 3
     if nNoEl == 8
       if   face == 1
           list = [1,2,3,4];
       elseif face == 2
           list = [5,8,7,6];
       elseif face == 3
           list = [1,5,6,2];
       elseif face == 4
           list = [2,3,7,6];
       elseif face == 5
           list = [3,7,8,4];
       elseif face == 6
           list = [4,8,5,1];
       end
     elseif nNoEl == 20
       if   face == 1
           list = [1,2,3,4,9,10,11,12];
       elseif face == 2
           list = [5,8,7,6,16,15,14,13];
       elseif face == 3
           list = [1,5,6,2,17,13,18,9];
       elseif face == 4
           list = [2,6,7,3,18,14,19,10];
       elseif face == 5
           list = [3,7,8,4,19,15,20,11];
       elseif face == 6
           list = [4,8,5,1,20,16,17,12];
       end
     end
   end
end