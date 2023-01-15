clear
clc
close all
%% material and workpiece variables
W=0.03;%width of workpiece
H=W;%height of workpiece
T=0.03;%thickness of workpiece
volf=0.06;%volume fracture of layer
W1=W*volf;%width of layer in workpiece
u1=1;%shear modules of matrix
uc=10;%shear contrast of layer to matrix
u2=u1*uc;%shear modules of layer
ku1=1000;%incompressibility of matrix
ku2=1000;%incompressibility of layer
k1=ku1*u1;%bulk modulus of matrix
k2=ku2*u2;%bulk modulus of layer
%Total no. material parameters, && list of parameters
materialProps_M=[u1;k1;];%List of matrix material parameters
materialProps_L=[u2;k2;];%List of layer material parameters
%% shape functions and boundary conditions
% order of element (1-Linear or 2-Quadratic)
order=2;

%Neumann boundary condition at y=L(top of the workpiece), this will be
%tranformed into displacement boundary in computation.
yVeloc=-1;
% Dirichlet boundary condition at y=0(bottom of the workpiece)
zDisp=0;%at z=0
yDisp=0;%at y=0
xDisp=0;%at x=0,y=0

%% Meshing parameters
nElx_M=8; %%number of elements in x in matrix
nElx_L=2; %%number of elements in x in layer
nEly=8;   %%number of elements in y(the division of elements in y are same in two phase)
nElz=4;   %%number of elements in z(the division of elements in z are same in two phase)

nDof=2;        %% nodal degree of freedom;%Number degrees of freedom per node (2 for 2D, 3 for 3D)
if nDof==2
    if(order==1)
        nNoEl=4;    %% number of Node in an 2D Element for order 1(Q4)
    end
    if(order==2)
        nNoEl=8;     %% number of Node in an 2D Element for order 2(Q8)
    end
else
    if(order==1)
        nNoEl=8;    %% number of Node in an 3D Element for order 1(H8)
    end
    if(order==2)
        nNoEl=20;     %% number of Node in an 3D Element for order 2(H20)
    end
end
%(here ndof=ncoord, but the program allows them to be different to allow extension to plate & beam elements with C^1 continuity)

%No. elements && connectivity
[nEle,nNodes,connArray,coorNoUD,nEleM,nEleL]=meshgeneration(nElx_M,nElx_L,nEly,nElz,W,W1,H,T,order,nDof);
%nNodes       Number of nodes
%nEle       Number of elements
%connArray       The node number of ith element's jth node in global node marking
%coorNoUD       The coordinates of ith node in Undeformed status, jth column is the jth coordinate
%nEleM       Number of elements in matrix
%nEleL       Number of elements in layer

%% Boundary conditions only 2D is included
%Number of nodes with prescribed displacements, with the prescribed displacements
if nDof==2
    j=1;
    for i=1:nNodes
        if(coorNoUD(i,1)==0 && coorNoUD(i,2)==0)  % at x=0,y=0 node
            dirich(j,1)=i;
            dirich(j,2)=1;
            dirich(j,3)=xDisp;
            dirich(j+1,1)=i;
            dirich(j+1,2)=2;
            dirich(j+1,3)=yDisp;
            j=j+2;
        end
        if(coorNoUD(i,1)~=0 && coorNoUD(i,2)==0)  %dirichlet at other y=0 nodes
            dirich(j,1)=i;
            dirich(j,2)=2;
            dirich(j,3)=yDisp;
            j=j+1;
        end
    end
    %Neuman boundary condition, which will be transformed into
    %displacement boudanry later
    j=1;
    for i=1:nNodes
        if(coorNoUD(i,2)==H)  %Neumann at y=H nodes
            neuman(j,1)=i;
            neuman(j,2)=2;
            neuman(j,3)=yVeloc;
            j=j+1;
        end
    end
    nDiri=length(dirich);
    nNeu=length(neuman);
else
    j=1;
    for i=1:nNodes
        if(coorNoUD(i,1)==0)  % fix x=0 nodes to fix a plane
            dirich(j,1)=i;
            dirich(j,2)=1;
            dirich(j,3)=xDisp;
            j=j+1;
        end
        if(coorNoUD(i,2)==0 )  % fix y=0 nodes to fix a plane
            dirich(j,1)=i;
            dirich(j,2)=2;
            dirich(j,3)=yDisp;
            j=j+1;
        end
        if(coorNoUD(i,3)==0 )  % fix z=0 nodes to fix a plane
            dirich(j,1)=i;
            dirich(j,2)=3;
            dirich(j,3)=zDisp;
            j=j+1;
        end
    end
    %Neuman boundary condition, which will be transformed into
    %displacement boudanry later
    j=1;
    for i=1:nNodes
        if(coorNoUD(i,2)==H)  %Neumann at y=H nodes
            neuman(j,1)=i;
            neuman(j,2)=2;
            neuman(j,3)=yVeloc;
            j=j+1;
        end
    end
    nDiri=length(dirich);
    nNeu=length(neuman);
end

%No. loaded element faces, with the loads
nDload=0;%Total number of element faces subjected to tractions
dLoads=[];%List of element tractions
%dloads(j,1) Element number
%dloads(j,2) face number
%dloads(j,3), dloads(j,4), dloads(j,5) Components of traction(assumed
%uniform)£¬5th component will not appears in 2D problem

% Plot the initial mesh for check
meshplotFunc(nDof,nEle,nNoEl,connArray,coorNoUD,1)

%% ============================ MAIN FEM ANALYSIS PROCEDURE ========================
%
%   w           Nodal displacements.  Let w_i^a be ith displa           cement component
%               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
%               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
%   dw          Correction to nodal displacements.  Let w_i^a be ith displacement component
%               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
%               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
%   K           Global stiffness matrix.  Stored as [K_1111 K_1112  K_1121  K_1122...
%                                                    K_1211 K_1212  K_1221  K_1222...
%                                                    K_2111 K_2112  K_2121  K_2122...]
%               for 2D problem and similarly for 3D problem
%   F           Force vector.  Currently only includes contribution from tractions
%               acting on element faces (i.e. body forces are neglected)
%   R           Volume contribution to residual
%   b           RHS of equation system
%
w = zeros(nNodes*nDof,1);
%
%  Newton Raphson iteration
%  Load is applied in nsteps increments
%  tol is the tolerance used in checking Newton-Raphson convergence
%  maxit is the max no. Newton-Raphson iterations
%  relax is the relaxation factor (Set to 1 unless big convergence problems)
%
%  Substeps for loading in each step
nsteps = 5;
tol = 0.0001;%tolerance
maxit = 30;%maximum iteration
relax = 1.;

%  time step parameters
tDelta=0.02;   %% ¦¤t
tSteps=10;    %% number of time steps

wtstep=zeros(nNodes*nDof,tSteps*nsteps);
kmin=zeros(1,tSteps*nsteps);

for tk=1:1:tSteps
    fprintf(1,'\n Time %f \n',tk*tDelta);
    if tk == 1
        coorNo = coorNoUD;
    else
        coorNo = coorNoUD + dispcoor(wtstep(:,nsteps*(tk-1)),nDof,nNodes);
    end
    for step=1:1:nsteps
        loadfactor=step/nsteps;
        err1=1.;
        nit=0;%number of iteration
        
        fprintf(1,'\n Step %f Load %f\n',step,loadfactor);
        while ((err1>tol) && (nit<maxit))          % Newton Raphson loop
            nit=nit + 1;
            K = globalstiffness(nDof,nEle,nNodes,nNoEl,connArray,coorNo,materialProps_M,materialProps_L,nEleM,nEleL,w);
            F = globaltraction(nDof,nNodes,nNoEl,connArray,coorNo,nDload,dLoads,w);
            R = globalresidual(nDof,nEle,nNodes,nNoEl,connArray,coorNo,materialProps_M,materialProps_L,nEleM,nEleL,w);
            b = loadfactor*F - R;
            
            %Fix constrained nodes to apply Dirichlet boundary condition
            for n=1:1:nDiri%loop over every Dirichlet boundary condition
                rw = nDof*(dirich(n,1)-1) + dirich(n,2);
                for cl=1:1:nDof*nNodes%loop to set all columns of row as 0
                    K(rw,cl) = 0;
                end
                K(rw,rw) = 1.;%set K(rw,rw) as 1
                b(rw) = loadfactor*dirich(n,3)-w(rw);
            end
            %Transfer Neuman boundary condition into node displacement to apply
            for n=1:1:nNeu
                rw = nDof*(neuman(n,1)-1) + neuman(n,2);
                for cl=1:1:nDof*nNodes
                    K(rw,cl) = 0;
                end
                K(rw,rw) = 1.;%set K(rw,rw) as 1
                b(rw) = loadfactor*neuman(n,3)*H*tk*tDelta-w(rw);
            end
            %Solve for the correction
            dw = K\b;
            %Check convergence
            w = w + relax*dw;
            wnorm = dot(w,w);
            err1 = dot(dw,dw);
            err2 = dot(b,b);
            err1 = sqrt(err1/wnorm);
            err2 = sqrt(err2)/(nDof*nNodes);
            fprintf(1,'Iteration number %d Correction %f Residual %f tolerance %f\n',nit,err1,err2,tol);
        end
        
        wtstep(:,nsteps*(tk-1)+step) = w;
        kmin(1,nsteps*(tk-1)+step) = min(real(eig(K)));
    end
    %
    %================================= IN-PROCESS DISPLAY =================================
    %
    % Create a plot of the deformed mesh
    meshdeformplotFunc(nDof,nEle,nNoEl,connArray,coorNoUD,(coorNoUD + dispcoor(w,nDof,nNodes)),tk+1,tk*tDelta)
end
%
%================================= POST-PROCESSING =================================
%
% Create a plot of the deformed mesh
figure(tk+2)
plot(1:1:tSteps*nsteps,kmin,'r','LineWidth',3);
grid on;
xlabel({'loop'},'FontSize',20);
ylabel({'min(eig(K))'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
title('Eigen value of stiffness matrix');
