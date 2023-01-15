close all
clear
clc

%% material and workpiece variables
W=0.03;%width of workpiece
H=W;%hight of workpiece
t=0.01;%thickness
volf=0.1;%volume fracture of layer
W1=W*volf;%width of layer in workpiece

E1=100;  % youngs modulus for matrix
E2=1000; % youngs modulus for layer
nu=0.3;   % poisson's ratio
rho=1;    % density is same for two phases

%% C matrix for linear elasticity
% Lame's parameters for matrix
%lambda1=E1/2/(1+nu);
%mu1=E1*nu/(1+nu)/(1-2*nu);
% Lame's parameters for layer
%lambda2=E2/2/(1+nu);
%mu2=E2*nu/(1+nu)/(1-2*nu);

%Cijkl Elastic moduli
% [C1]=Celasticm(lambda1,mu1,delta);
% [C2]=Celasticm(lambda2,mu2,delta);

%% shape functions and boundary conditions
% order of shape function (1-Linear or 2-Quadratic), only have Linear now
order=2;

% Gaussian quadrature  %%% select from 1, 2 or 3
wq=3;

%Neumann boundary condition at y=L(top of the workpiece)
yVeloc=-1;
% Dirichlet boundary condition at y=0(top of the workpiece)
yDisp=0;%at y=0
xDisp=0;%at x=0,y=0

%% time step parameters
beta=0.255;
gamma=2*beta;
tDelta=0.01;   %% deltaT
tSteps=10;    %% number of time steps

%% Meshing parameters
nElx_M=4; %%number of elements in x in matrix
nElx_L=1; %%number of elements in x in layer
nEly=4;   %%number of elements in y(the division of elements in y are same in two phase)

if(order==2)
    nNoEl=8;    %% number of Node in an Element for order 2
else
    nNoEl=4;     %% number of Node in an Element for order 1
end

nDof=2;        %% degree of freedom;

%% Mesh generation
[nElements,nNodes,connArray,coorNoUD,ElementsM,ElementsL]=meshgeneration(nElx_M,nElx_L,nEly,W,W1,H,order);

%Dirichlet boundary condition
j=1;
for i=1:nNodes
    if(coorNoUD(i,1)==0 && coorNoUD(i,2)==0)  % at x=0,y=0 node
        dirich(j,1:2)=[2*(i-1)+1,xDisp];
        dirich(j+1,1:2)=[2*(i-1)+2,yDisp];
        j=j+2;
    end
    if(coorNoUD(i,1)~=0 && coorNoUD(i,2)==0)  %dirichlet at other y=0 nodes
        dirich(j,1:2)=[2*(i-1)+2,yDisp];
        j=j+1;
    end
end
%Neumann boundary condition
j=1;
for i=1:nNodes
    if(coorNoUD(i,2)==H)  %Neumann at y=H nodes
        neuman(j,1:2)=[2*(i-1)+2,yVeloc];
        j=j+1;
    end
end
nDirich=length(dirich);
nNeuman=length(neuman);

%% Initial conditions: displacement and velocity
aDof=nDof*nNodes;
u0=zeros(aDof,1);
v0=zeros(aDof,1);
X0=[v0;u0];
u=zeros(aDof,tSteps);%store the solved node displacements and velocities
v=zeros(aDof,tSteps);

IndexD=zeros(aDof,1);
IndexN=zeros(aDof,1);
for i=1:1:aDof
    IndexD(i)=i;
    IndexN(i)=i;
end
IndexD(dirich(:,1))=[];
IndexN(neuman(:,1))=[];
%% Time-dependent solution
utemp=zeros(aDof,1);
vtemp=zeros(aDof,1);
for tk=1:1:tSteps
    if tk ==1
        Xtk=[v0(IndexN);u0(IndexD)];
        coorNoD_tk=coorNoUD+indexT(u0,nDof,nNodes);
    else
        Xtk=[v(IndexN,tk-1);u(IndexD,tk-1)];
        coorNoD_tk=coorNoUD+indexT(u(tk-1),nDof,nNodes);
    end
    %local stiffness matrix and global stiffness matrix
    [kGlobal,mGlobal,fGlobal]=globalmatrix(nDof,nNodes,nNoEl,nElements,ElementsM,ElementsL,wq,order,connArray,coorNoD_tk,E1,E2,nu,rho,t);
    %kD, mD and fD, add Dirichlet boundary condition
    [kBD,mBD,fBD]=boundarycondition(nDof,nNodes,nDirich,dirich,nNeuman,neuman,kGlobal,mGlobal,fGlobal,H);

    Xtk1b=Xtk+tDelta*inv(mBD)*(-kBD*Xtk+fBD);
    
    utemp(IndexD)=Xtk1b(aDof-nNeuman+1:2*aDof-nNeuman-nDirich);
    utemp(dirch(:,1))=0;
    utemp(neuman(:,1))=yVeloc*tk*tDelta;
    coorNoD_tk=coorNoUD+indexT(utemp,nDof,nNodes);
    
    [kGl,mGl,fGl]=globalmatrix(nDof,nNodes,nNoEl,nElements,ElementsM,ElementsL,wq,order,connArray,coorNoD_tk,E1,E2,nu,rho,t);
    [kBD1,mBD1,fBD1]=boundarycondition(nDof,nNodes,nDirich,dirich,nNeuman,neuman,kG,mGl,fGl,H);
    Xtk1=Xtk+0.5*tDelta*(inv(mBD)*(-kBD*Xtk+fBD)+inv(mBD1)*(-kBD1*Xtk1b+fBD1));
    
end