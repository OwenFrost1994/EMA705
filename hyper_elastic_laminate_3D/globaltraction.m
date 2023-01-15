%
%===================== Assemble the global traction vector =============
%
function F=globaltraction(nDof,nNodes,nNoEl,connArray,coorNo,nDload,dLoads,disp)
%nDof            node degree of freedom 2-2D,3-3D
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%coorNo          coordinates of ith node in previous status
%nDload          number of traction loads
%dLoads          list of load: nodes, directions and amplitude
%disp            the displacements of nodes
%
%Generate a traction vector through discretization
%
F=zeros(nDof*nNodes,1);
traction=zeros(nDof,1);

for load=1:1:nDload
%
%     Extract the coords of the nodes on the appropriate element face
%
      lmn = dLoads(load,1);
      face = dLoads(load,2);
      nfnodes = nfacenodes(nDof,nNoEl); 
      nodelist = facenodes(nDof,nNoEl,face);     
      lmncoord = zeros(nfnodes,nDof);
      for a = 1:nfnodes
        for i = 1:nDof
          lmncoord(a,i) = coorNo(connArray(nodelist(a),dLoads(load,1)),i);
        end
        for i = 1:nDof
          lmndof(a,i) = disp(nDof*(connArray(nodelist(a),dLoads(load,1))-1)+i);
        end
      end
%
%    Compute the element load vector
%
     for i = 1:nDof
       traction(i) = dLoads(load,i+2);
     end

     rel = eldload(nDof,nfnodes,lmncoord,traction);
%
%    Assemble the element load vector into global vector
%
     for a = 1:nfnodes
       for i = 1:nDof
         rw = (connArray(nodelist(a),dLoads(load,1))-1)*nDof+i;
         F(rw) = F(rw) + rel((a-1)*nDof+i);
       end
     end

end
end       