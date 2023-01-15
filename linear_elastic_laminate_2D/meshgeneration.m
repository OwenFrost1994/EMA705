function [nElements,nNodes,connArray,coorNoUD,ElementsM,ElementsL]=meshgeneration(nElx_M,nElx_L,nEly,W,W1,H,order)
%this function works as the generator of mesh
%this function will calculate the number of nodes, number of elements, coordinates of nodes before the deformation,
%the connection array,in-matrix element index, in-layer element index
if(order==2)%contemporarily don't have
    nNodes=(nElx_M+nElx_L+1)*(nEly+1)+(nElx_M+nElx_L+1)*(nEly)+(nElx_M+nElx_L)*(nEly+1);%number of nodes
    nElements=(nElx_M+nElx_L)*nEly;%number of elements
    elemW_M=(W-W1)/nElx_M;%width of element in matrix
    elemH_M=H/nEly;%hight of element in matrix
    elemW_L=W1/nElx_L;%width of element in layer
    elemH_L=H/nEly;%hight of element in layer
    %generate a matrix store the coordinates of nodes before the deformation
    coorNoUD=zeros(nNodes,2);
    for i=1:2*(nElx_M+nElx_L)+1
        if mod(i,2) ~= 0
            for j=1:2*nEly+1
                idxnode=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2))+j;
                if i<=nElx_M+1
                    coorNoUD(idxnode,1)=(i-1)*elemW_M/2;
                else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                         coorNoUD(idxnode,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                    else
                         coorNoUD(idxnode,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                    end
                end
                coorNoUD(idxnode,2)=(2*nEly+1-j)*elemH_L/2;
            end
        else
            for j=1:nEly+1
                idxnode=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2)-1)+j;
                if i<=nElx_M+1
                    coorNoUD(idxnode,1)=(i-1)*elemW_M/2;
                else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                         coorNoUD(idxnode,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                    else
                         coorNoUD(idxnode,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                    end
                end
                coorNoUD(idxnode,2)=(nEly+1-j)*elemH_L;
            end
        end
    end
    
    % generate the connectivity array
    connArray=zeros(nElements,8);
    nEleM=1;
    nEleL=1;
    ElementsM=zeros(nElx_M*nEly,1);
    ElementsL=zeros(nElx_L*nEly,1);
    for i=1:nElx_M+nElx_L
        for j=1:nEly
            k=j+(i-1)*nEly;
            connArray(k,1)=1+2*j+(3*nEly+2)*(i-1);
            connArray(k,2)=1+2*j+(3*nEly+2)*i;
            connArray(k,3)=1+2*(j-1)+(3*nEly+2)*i;
            connArray(k,4)=1+2*(j-1)+(3*nEly+2)*(i-1);
            connArray(k,5)=1+j+(2*nEly+1)*i+(nEly+1)*(i-1);
            connArray(k,6)=2*j+(3*nEly+2)*i;
            connArray(k,7)=j+(2*nEly+1)*i+(nEly+1)*(i-1);
            connArray(k,8)=2*j+(3*nEly+2)*(i-1);
            if nElx_M/2<i && i<=nElx_M/2+nElx_L
                ElementsL(nEleL,1)=k;
                nEleL=nEleL+1;
            else
                ElementsM(nEleM,1)=k;
                nEleM=nEleM+1;
            end
        end
    end
end
if(order==1)
     nNodes=(nElx_M+nElx_L+1)*(nEly+1);%number of nodes
     nElements=(nElx_M+nElx_L)*nEly;%number of elements
     elemW_M=(W-W1)/nElx_M;%width of element in matrix
     elemH_M=H/nEly;%hight of element in matrix
     elemW_L=W1/nElx_L;%width of element in layer
     elemH_L=H/nEly;%hight of element in layer
     %generate a matrix store the coordinates of nodes before the deformation
     coorNoUD=zeros(nNodes,2);
     for i=1:nElx_M+nElx_L+1
         for j=1:nEly+1
             if i<=nElx_M/2+1
                 coorNoUD(j+(i-1)*(nEly+1),1)=(i-1)*elemW_M;
             else if nElx_M/2+1<i && i<nElx_M/2+nElx_L+1
                     coorNoUD(j+(i-1)*(nEly+1),1)=(W-W1)/2+(i-nElx_M/2-1)*elemW_L;
                 else
                     coorNoUD(j+(i-1)*(nEly+1),1)=(W+W1)/2+(i-nElx_M/2-nElx_L-1)*elemW_M;
                 end
             end
             coorNoUD(j+(i-1)*(nEly+1),2)=(nEly+1-j)*elemH_L;
         end
     end
     % generate the connectivity array
     connArray=zeros(nElements,4);
     nEleM=1;
     nEleL=1;
     ElementsM=zeros(nElx_M*nEly,1);
     ElementsL=zeros(nElx_L*nEly,1);
     for i=1:nElx_M+nElx_L
         for j=1:nEly
             k=j+(i-1)*nEly;
             connArray(k,1)=1+j+(i-1)*(nEly+1);
             connArray(k,2)=1+j+i*(nEly+1);
             connArray(k,3)=j+i*(nEly+1);
             connArray(k,4)=j+(i-1)*(nEly+1);
             if nElx_M/2<i && i<=nElx_M/2+nElx_L
                 ElementsL(nEleL,1)=k;
                 nEleL=nEleL+1;
             else                 
                 ElementsM(nEleM,1)=k;
                 nEleM=nEleM+1;
             end
         end
     end
end
end