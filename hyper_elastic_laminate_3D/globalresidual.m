%
%====================== Assemble the global residual vector =================
%
function resid = globalresidual(nDof,nEle,nNodes,nNoEl,connArray,coordNo,materialProps_M,materialProps_L,nEleM,nEleL,disp)
%nDof            node degree of freedom 2-2D,3-3D
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%coorNo          coordinates of ith node in previous status
%materialprops   number of traction loads
%disp            the displacements of nodes
%
%Assemble the global residual matrix
%

resid = zeros(nDof*nNodes,1);
lmncoord = zeros(nNoEl,nDof);
lmndisp = zeros(nNoEl,nDof);
rel = zeros(nDof*nNoEl,nDof*nNoEl);
%
%   Loop over all the elements
%
   for lmn=1:1:nEle
%
%   Extract coords of nodes, DOF for the current element
%
      for a=1:1:nNoEl
        for i=1:1:nDof
          lmncoord(a,i) = coordNo(connArray(lmn,a),i);
        end
        for i=1:1:nDof
          lmndisp(a,i) = disp(nDof*(connArray(lmn,a)-1)+i);
        end
      end
      if ismember(lmn,nEleM)
          rel=elresid(nDof,nNoEl,lmncoord,materialProps_M,lmndisp);
      end
      if ismember(lmn,nEleL)
          rel=elresid(nDof,nNoEl,lmncoord,materialProps_L,lmndisp);
      end

%
%   Add the current element residual to the global residual
%
      for a=1:1:nNoEl
        for i=1:1:nDof
          rw = nDof*(connArray(lmn,a)-1)+i;
          resid(rw) = resid(rw) + rel(nDof*(a-1)+i);
        end
     end
  end
 
end