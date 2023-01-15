function [nElements,nNodes,connArray,coorNoUD,ElementsM,ElementsL]=meshgeneration(nElx_M,nElx_L,nEly,nElz,W,W1,H,T,order,nDof)
%this function works as the generator of mesh
%this function will calculate the number of nodes, number of elements, coordinates of nodes before the deformation,
%the connection array,in-matrix element index, in-layer element index
if nDof ==2 %2D problem
    if(order==2)
        nNodes=(nElx_M+nElx_L+1)*(nEly+1)+(nElx_M+nElx_L+1)*(nEly)+(nElx_M+nElx_L)*(nEly+1);%number of nodes
        nElements=(nElx_M+nElx_L)*nEly;%number of elements
        elemW_M=(W-W1)/nElx_M;%width of element in matrix
        elemW_L=W1/nElx_L;%width of element in layer
        elemH=H/nEly;%hight of element
        %generate a matrix store the coordinates of nodes before the deformation
        coorNoUD=zeros(nNodes,2);
        for i=1:2*(nElx_M+nElx_L)+1
            if mod(i,2) ~= 0
                for j=1:2*nEly+1
                    nid=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2))+j;
                    if i<=nElx_M+1
                        coorNoUD(nid,1)=(i-1)*elemW_M/2;
                    else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                            coorNoUD(nid,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                        else
                            coorNoUD(nid,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                        end
                    end
                    coorNoUD(nid,2)=(2*nEly+1-j)*elemH/2;
                end
            else
                for j=1:nEly+1
                    nid=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2)-1)+j;
                    if i<=nElx_M+1
                        coorNoUD(nid,1)=(i-1)*elemW_M/2;
                    else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                            coorNoUD(nid,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                        else
                            coorNoUD(nid,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                        end
                    end
                    coorNoUD(nid,2)=(nEly+1-j)*elemH;
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
                eid=j+(i-1)*nEly;
                connArray(eid,1)=1+2*j+(3*nEly+2)*(i-1);
                connArray(eid,2)=1+2*j+(3*nEly+2)*i;
                connArray(eid,3)=1+2*(j-1)+(3*nEly+2)*i;
                connArray(eid,4)=1+2*(j-1)+(3*nEly+2)*(i-1);
                connArray(eid,5)=1+j+(2*nEly+1)*i+(nEly+1)*(i-1);
                connArray(eid,6)=2*j+(3*nEly+2)*i;
                connArray(eid,7)=j+(2*nEly+1)*i+(nEly+1)*(i-1);
                connArray(eid,8)=2*j+(3*nEly+2)*(i-1);
                if nElx_M/2<i && i<=nElx_M/2+nElx_L
                    ElementsL(nEleL,1)=eid;
                    nEleL=nEleL+1;
                else
                    ElementsM(nEleM,1)=eid;
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
        elemH=H/nEly;%hight of element in layer
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
                coorNoUD(j+(i-1)*(nEly+1),2)=(nEly+1-j)*elemH;
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
                eid=j+(i-1)*nEly;
                connArray(eid,1)=1+j+(i-1)*(nEly+1);
                connArray(eid,2)=1+j+i*(nEly+1);
                connArray(eid,3)=j+i*(nEly+1);
                connArray(eid,4)=j+(i-1)*(nEly+1);
                if nElx_M/2<i && i<=nElx_M/2+nElx_L
                    ElementsL(nEleL,1)=eid;
                    nEleL=nEleL+1;
                else
                    ElementsM(nEleM,1)=eid;
                    nEleM=nEleM+1;
                end
            end
        end
    end
else%3D problem
    if (order==2)
        nNodes=(nElx_M+nElx_L+1)*(nEly+1)*(nElz+1)+(nElx_M+nElx_L)*(nEly+1)*(nElz+1)+(nElx_M+nElx_L+1)*(nEly)*(nElz+1)+(nElx_M+nElx_L+1)*(nEly+1)*(nElz);%number of nodes
        nElements=(nElx_M+nElx_L)*nEly*nElz;%number of elements
        elemW_M=(W-W1)/nElx_M;%width of element in matrix
        elemW_L=W1/nElx_L;%width of element in layer
        elemH=H/nEly;%hight of element
        elemT=T/nElz;%hight of element
        %generate a matrix store the coordinates of nodes before the deformation
        coorNoUD=zeros(nNodes,3);
        for i=1:2*(nElx_M+nElx_L)+1%%loop of x axis
            if mod(i,2) ~= 0%%ÆæÊý
                if i<=nElx_M+1
                    nX=(i-1)*elemW_M/2;
                else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                        nX=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                    else
                        nX=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                    end
                end
                for k=1:2*nElz+1%%loop of z axis
                    nZ=(2*nElz+1-k)*elemT/2;
                    if mod(k,2) ~= 0%%ÆæÊý
                        for j=1:2*nEly+1%%loop of y axis
                            nY=(2*nEly+1-j)*elemH/2;
                            nid=((2*nEly+1)*(nElz+1)+(nEly+1)*nElz)*(floor(i/2))+(nEly+1)*(nElz+1)*(floor(i/2))...
                                +(2*nEly+1)*(floor(k/2))+(nEly+1)*(floor(k/2))...
                                +j;
                            coorNoUD(nid,1:3)=[nX,nY,nZ];
                        end
                    else%%Å¼Êý
                        for j=1:nEly+1%%loop of y axis
                            nY=(nEly+1-j)*elemH;
                            nid=((2*nEly+1)*(nElz+1)+(nEly+1)*nElz)*(floor(i/2))+(nEly+1)*(nElz+1)*(floor(i/2))...
                                +(2*nEly+1)*(floor(k/2))+(nEly+1)*(floor(k/2)-1)...
                                +j;
                            coorNoUD(nid,1:3)=[nX,nY,nZ];
                        end
                    end
                end
            else
                if i<=nElx_M+1
                    nX=(i-1)*elemW_M/2;
                else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                        nX=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                    else
                        nX=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                    end
                end
                for k=1:nElz+1%%loop of z axis
                    nZ=(nElz+1-k)*elemT;
                    for j=1:nEly+1
                        nY=(nEly+1-j)*elemH;
                        nid=((2*nEly+1)*(nElz+1)+(nEly+1)*nElz)*(floor(i/2))+(nEly+1)*(nElz+1)*(floor(i/2)-1)...
                            +(nEly+1)*(k-1)+j;
                        coorNoUD(nid,1:3)=[nX,nY,nZ];
                    end
                end
            end
        end
        % generate the connectivity array
        connArray=zeros(nElements,20);
        nEleM=1;
        nEleL=1;
        ElementsM=zeros(nElx_M*nEly*nElz,1);
        ElementsL=zeros(nElx_L*nEly*nElz,1);
        for i=1:nElx_M+nElx_L
            for j=1:nEly
                for k=1:nElz
                    eid=j+nEly*(k-1)+nEly*nElz*(i-1);
                    connArray(eid,1)=1+2*j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,2)=1+2*j+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,3)=1+2*(j-1)+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,4)=1+2*(j-1)+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    
                    connArray(eid,5)=1+2*j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,6)=1+2*j+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,7)=1+2*(j-1)+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,8)=1+2*(j-1)+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    
                    connArray(eid,9)=2*nEly+1+1+j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,10)=2*j+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,11)=2*nEly+1+j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    connArray(eid,12)=2*j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1);
                    
                    connArray(eid,13)=2*nEly+1+1+j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,14)=2*j+(3*nEly+2)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,15)=2*nEly+1+j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    connArray(eid,16)=2*j+(3*nEly+2)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*i;
                    
                    connArray(eid,17)=1+j+(nEly+1)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1)+(2*nEly+1)*(nElz+1)+(nEly+1)*nElz;
                    connArray(eid,18)=1+j+(nEly+1)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1)+(2*nEly+1)*(nElz+1)+(nEly+1)*nElz;
                    connArray(eid,19)=j+(nEly+1)*k+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1)+(2*nEly+1)*(nElz+1)+(nEly+1)*nElz;
                    connArray(eid,20)=j+(nEly+1)*(k-1)+((2*nEly+1)*(nElz+1)+(nEly+1)*nElz+(nEly+1)*(nElz+1))*(i-1)+(2*nEly+1)*(nElz+1)+(nEly+1)*nElz;
                    if nElx_M/2<i && i<=nElx_M/2+nElx_L
                        ElementsL(nEleL,1)=eid;
                        nEleL=nEleL+1;
                    else
                        ElementsM(nEleM,1)=eid;
                        nEleM=nEleM+1;
                    end
                end
            end
        end
    end
    if(order==1)
        nNodes=(nElx_M+nElx_L+1)*(nEly+1)*(nElz+1);%number of nodes
        nElements=(nElx_M+nElx_L)*nEly*nElz;%number of elements
        elemW_M=(W-W1)/nElx_M;%width of element in matrix
        elemW_L=W1/nElx_L;%width of element in layer
        elemH=H/nEly;%hight of element
        elemT=T/nElz;%hight of element
        %generate a matrix store the coordinates of nodes before the deformation
        coorNoUD=zeros(nNodes,3);
        for i=1:1:nElx_M+nElx_L+1
            for j=1:1:nEly+1
                for k=1:1:nElz+1
                    nid=j+(k-1)*(nEly+1)+(i-1)*(nEly+1)*(nElz+1);
                    if i<=nElx_M/2+1
                        coorNoUD(nid,1)=(i-1)*elemW_M;
                    else if nElx_M/2+1<i && i<nElx_M/2+nElx_L+1
                            coorNoUD(nid,1)=(W-W1)/2+(i-nElx_M/2-1)*elemW_L;
                        else
                            coorNoUD(nid,1)=(W+W1)/2+(i-nElx_M/2-nElx_L-1)*elemW_M;
                        end
                    end
                    coorNoUD(nid,2)=(nEly+1-j)*elemH;
                    coorNoUD(nid,3)=(nElz+1-k)*elemT;
                end
            end
        end
        % generate the connectivity array
        connArray=zeros(nElements,8);
        nEleM=1;
        nEleL=1;
        ElementsM=zeros(nElx_M*nEly*nElz,1);
        ElementsL=zeros(nElx_L*nEly*nElz,1);
        for i=1:nElx_M+nElx_L
            for j=1:nEly
                for k=1:1:nElz
                    eid=j+(k-1)*nEly+(i-1)*nEly*nElz;
                    connArray(eid,1)=1+j+(k-1)*(nEly+1)+(i-1)*(nEly+1)*(nElz+1);
                    connArray(eid,2)=1+j+k*(nEly+1)+(i-1)*(nEly+1)*(nElz+1);
                    connArray(eid,3)=j+k*(nEly+1)+(i-1)*(nEly+1)*(nElz+1);
                    connArray(eid,4)=j+(k-1)*(nEly+1)+(i-1)*(nEly+1)*(nElz+1);
                    connArray(eid,5)=1+j+(k-1)*(nEly+1)+i*(nEly+1)*(nElz+1);
                    connArray(eid,6)=1+j+k*(nEly+1)+i*(nEly+1)*(nElz+1);
                    connArray(eid,7)=j+k*(nEly+1)+i*(nEly+1)*(nElz+1);
                    connArray(eid,8)=j+(k-1)*(nEly+1)+i*(nEly+1)*(nElz+1);
                    if nElx_M/2<i && i<=nElx_M/2+nElx_L
                        ElementsL(nEleL,1)=eid;
                        nEleL=nEleL+1;
                    else
                        ElementsM(nEleM,1)=eid;
                        nEleM=nEleM+1;
                    end
                end
            end
        end
    end
end