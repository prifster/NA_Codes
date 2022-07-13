set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesFontName','Times New Roman')

A = 100 * 1e-6;             % Area of Cross Section
L = 2;                      % Dimension L
E = 200 * 1e9;              % Young's Modulus of Elasticity
F13 = 20 * 1e3;             % Force in horizontal direction - F13
F14 = 20 * 1e3;             % Force in vertical direction - F14

nmeb = 15;                  % Number of members
ndof = 4;                   % Degrees of freedom per member

% Degree of Freedom Connectivity Matrix - Jth column of the matrix
% corresponds to the global degrees of freedom of the Jth element/member

b = zeros(ndof,nmeb);

b(:,1)  = [1;2;7;8];
b(:,2)  = [3;4;5;6];
b(:,3)  = [1;2;5;6];
b(:,4)  = [3;4;7;8];
b(:,5)  = [5;6;7;8];
b(:,6)  = [5;6;11;12];
b(:,7)  = [7;8;9;10];
b(:,8)  = [5;6;9;10];
b(:,9)  = [7;8;11;12];
b(:,10) = [9;10;11;12];
b(:,11) = [9;10;15;16];
b(:,12) = [11;12;13;14];
b(:,13) = [9;10;13;14];
b(:,14) = [11;12;15;16];
b(:,15) = [13;14;15;16];

% Lengths and Orientations of members - Len and chi are 1 x 15 vectors
% containing the length and oreintation (in rad) of each member of the
% structure

Len = L*[1.414 1.414 1 1 1 1.414 1.414 1 1 1 1.414 1.414 1 1 1];
chi = pi*[0.25 0.75 0.5 0.5 0 0.25 0.75 0.5 0.5 0 0.25 0.75 0.5 0.5 0];

% Constraints

constraint = [1 2 3 4];

% Initialize matrices

ke = zeros(ndof,ndof,nmeb);     % Element stiffness matrix in element coordinates
Ke = zeros(ndof,ndof,nmeb);     % Element stiffness matrix in structure coordinates
Ndof = max(max(b));             % Total degrees of freedom of structure
K = zeros(Ndof);                % Global stiffness matrix in structure coordinates

% Transformation Matrix as a function of angle
    
T = @(x) [cos(x) sin(x) 0 0; -sin(x) cos(x) 0 0; 0 0 cos(x) sin(x); 0 0 -sin(x) cos(x)]


for i = 1:nmeb
    
    % Formulate the element stiffness matrix in element coordinates
    
    ke(1,1,i) = (A*E)/Len(i);
    ke(1,3,i) = -(A*E)/Len(i);
    ke(3,1,i) = -(A*E)/Len(i);
    ke(3,3,i) = (A*E)/Len(i);
    
    % Calculate the element stiffness matrix in structure coordinates

    Ke(:,:,i) = T(chi(i))'*ke(:,:,i)*T(chi(i));   
    
end

for i = 1:nmeb
    
    % Assemble the global stiffness matrix
    
    K(b(:,i),b(:,i))  = K(b(:,i),b(:,i)) + Ke(:,:,i);
end

% Reduce the stiffness matrix

red_ind = setdiff(1:Ndof,constraint);       % Non prescribed degrees of freedom
K_red = K(red_ind,red_ind);                 % Reduced stiffness matrix

% External Force Vector (reduced)

F_red = [zeros(8,1); F13; F14; zeros(2,1);];

% Compute the displacement vector

U = zeros(Ndof,1);
U(red_ind) = K_red\F_red;

% Calculate the element displacement and internal force

delta = zeros(ndof,nmeb);
F_int = zeros(nmeb,1);
for i = 1:nmeb
    delta(:,i) = T(chi(i))*U(b(:,i));
    F_int(i) = ((E*A)/Len(i))*(delta(3,i) - delta(1,i));
end

% Which member has the maximum internal force? Indicate the member number
% that has the maximum internal force

max_stress_member_number = max(F_int);


% Visualization ( You do not need to change any code below)

nodes = [0 0; L 0; 0 L; L L; 0 2*L; L 2*L; 0 3*L; L 3*L; ]';
conn = [1 4; 2 3; 1 3; 2 4; 3 4; 3 6; 4 5; 3 5; 4 6; 5 6; 5 8; 6 7; 5 7; 6 8; 7 8;]';

scale = 5;

figure(1)
for i = 1:nmeb
    temp1 = [nodes(:,conn(1,i)) nodes(:,conn(2,i))];
    temp2 = [nodes(:,conn(1,i)) + scale*U(2*(conn(1,i)-1)+(1:2)) nodes(:,conn(2,i)) + scale*U(2*(conn(2,i)-1)+(1:2))];
    plot(temp1(1,:),temp1(2,:),'k')
    hold all
    plot(temp2(1,:),temp2(2,:),'r')    
end
axis equal
box off
