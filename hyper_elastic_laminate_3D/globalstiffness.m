%====================== Assemble the global stiffness matrix =================
%
function K=globalstiffness(nDof,nEle,nNodes,nNoEl,connArray,coorNo,materialProps_M,materialProps_L,nEleM,nEleL,disp)%elident,nelnodes,
%nDof            node degree of freedom 2-2D,3-3D
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%coorNo          coordinates of ith node in previous status
%materialprops   material property
%disp            the displacements of nodes
%
%Generate the element stiffness matrix
%Assemble the global stiffness matrix
%
K = zeros(nDof*nNodes,nDof*nNodes);%empty global stiffness
lmncoord = zeros(nNoEl,nDof);
lmndis = zeros(nNoEl,nDof);
%
%Loop over all the elements to generate every stiffness 
for lmn = 1:nEle
%
%Extract coords of nodes, DOF for the current element
%
    for a = 1:nNoEl%loop of element nodes
        for i = 1:nDof%pick out node coordinates
          lmncoord(a,i)=coorNo(connArray(lmn,a),i);
        end
        for i = 1:nDof%pick out node displacement
          lmndis(a,i)=disp(nDof*(connArray(lmn,a)-1)+i);
        end
    end
    
    if ismember(lmn,nEleM)
        kel=elstif(nDof,nNoEl,lmncoord,lmndis,materialProps_M);%k elastic
    end
     if ismember(lmn,nEleL)
        kel=elstif(nDof,nNoEl,lmncoord,lmndis,materialProps_L);%k elastic
     end
     
    % Add the current element stiffness:the global stiffness
    for a = 1:nNoEl
      for i = 1:nDof
        for b = 1:nNoEl
          for k = 1:nDof
            rw = nDof*(connArray(lmn,a)-1)+i;
            cl = nDof*(connArray(lmn,b)-1)+k;
            K(rw,cl) = K(rw,cl) + kel(nDof*(a-1)+i,nDof*(b-1)+k);
          end
        end
      end
    end
end
end