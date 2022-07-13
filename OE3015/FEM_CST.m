set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesFontName','Times New Roman')


L = 5;                          % Length in meters
h = 1;                          % Depth of beam in meters

E = 200 * 1e9;                  % Youngs Modulus
nu = 0.3;                       % Poisson's Ratio
t = 5/1000;                    % Thickness
q = 10 * 1e3;                   % Uniformly distributed loading
G = E/2/(1+nu);                 % Shear Modulus

hsd = 20;                       % Number of horizontal subdivisions
vsd = 60;                       % Number of vertical subdivisions

xsd = linspace(0,L,vsd+1);      % x-coordinates along length
ysd = linspace(0,h,hsd+1);      % y-coordinates along height

nnodes = (hsd+1)*(vsd+1);       % Number of nodes
nmeb = hsd*2*vsd;               % Number of members

nodes = zeros(2,nnodes);        % Initialize nodes matrix
conn = zeros(3,nmeb);           % Initialize connectivity matrix

for i = 1:vsd+1
    for j = 1:hsd+1
        
        % Populate the nodes matrix
        
        node_ind = (i-1)*(hsd+1) + j;
        nodes(:,node_ind) = [xsd(i) ysd(j)];
        
        % Populate the connectivity matrix
        
        if i <= vsd && j <= hsd
            conn_ind1 = (i-1)*(2*hsd) + 2*(j-1) + 1;
            conn_ind2 = (i-1)*(2*hsd) + 2*(j-1) + 2;
            
            conn(:,conn_ind1) = [node_ind node_ind+hsd+1 node_ind+1];
            conn(:,conn_ind2) = [node_ind+1 node_ind+hsd+1 node_ind+hsd+2];
            
        end
    end
end

ndof = 2;       % Number of degrees of freedom per node

% Initialize DOF connectivity matrix

dof_conn = zeros(3*ndof,nmeb);

% Initialize global stiffness matrix

K = zeros(nnodes*ndof);


for i = 1:nmeb
    
    % Degree of Freedom Connectivity Matrix - Ith column of the matrix
    % corresponds to the global degrees of freedom of the Ith element/member
    
    dof_conn(:,i) = [2*conn(1,i)-1 2*conn(1,i) 2*conn(2,i)-1 2*conn(2,i) 2*conn(3,i)-1 2*conn(3,i)];
    
    % Find the nodes of the Ith member
    
    node1 = nodes(:,conn(1,i));
    node2 = nodes(:,conn(2,i));
    node3 = nodes(:,conn(3,i));
    
    % Calculate the stiffness matrix using the function we previously
    % computed
    
    Ke = cst(node1,node2,node3,E,nu,t);
    
    % Assemble the global stiffness matrix
    
    K(dof_conn(:,i),dof_conn(:,i)) = K(dof_conn(:,i),dof_conn(:,i)) + Ke;
    
end

% Formulate the force vector by assuming a discrete trapezoidal
% distribution across the vertical degrees of freedom of the top edge of
% the beam

F = zeros(ndof*nnodes,1);
f_ind = ndof*(hsd+1):ndof*(hsd+1):nnodes*ndof;

f_temp = ones(1,vsd+1);
f_temp(1) = 0.5; f_temp(end) = 0.5;

F(f_ind) = -q*L/vsd * f_temp;

% Constraints and Structure's degrees of freedom

constraints = [ndof*(hsd/2+1)-1 2:2:ndof*(hsd+1) (nnodes-hsd)*ndof:ndof:ndof*nnodes];
struct_dof = setdiff(1:nnodes*ndof,constraints);

% Compute the displacement vector

U = zeros(nnodes*ndof,1);
U(struct_dof) = K(struct_dof,struct_dof)\F(struct_dof);

fem_disp = U(ndof*(hsd/2+1):ndof*(hsd+1):ndof*nnodes);

I = h^3*t/12;
sbt_disp = -q*xsd/(24*E*I).*(L^3 - 2*L*xsd.^2 + xsd.^3);
et_disp = -q*xsd/(24*E*I).*(L^3 - 2*L*xsd.^2 + xsd.^3 - 18/h/t/G*E*I * (xsd - L));

figure(1)
plot(xsd,sbt_disp*1e3,'k',xsd,et_disp*1e3,'b--',xsd,fem_disp*1e3,'r')
xlabel('x in meters')
ylabel('\delta in mm')
title('FEM with CST elements')
legend('Simple Beam Theory','Theory of Elasticity','FEM','location','best')
legend boxoff
hold all
