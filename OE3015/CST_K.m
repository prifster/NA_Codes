%% Code to call function:
E = 200 * 1e9;
nu = 0.3;
t = 5/1000;

node1 = [0 0];
node2 = [1 0];
node3 = [1 1];

Ke = cst(node1, node2, node3, E, nu, t)

%% Code:

function Ke = cst(node1, node2, node3, E, nu, t)

x1 = node1(1); y1 = node1(2);
x2 = node2(1); y2 = node2(2);
x3 = node3(1); y3 = node3(2);

% Area of the element
one = [1 1 1];
A123 = 0.5*det([one; x1 x2 x3; y1 y2 y3]);

% Formulate the A matrix

A = [1 x1 y1 0 0 0; 0 0 0 1 x1 y1; 1 x2 y2 0 0 0; 0 0 0 1 x2 y2; 1 x3 y3 0 0 0; 0 0 0 1 x3 y3];

% Compute the inverse of A matrix

Ainv = inv(A);

% Relate the strains to internal displacement using the G matrix

G = [0 1 0 0 0 0; 0 0 0 0 0 1; 0 0 1 0 1 0];

% Relate the strains to nodal displacements using the B matrix

B = G*Ainv;

% Formulate the D matrix

D = (E/(1-(nu^2)))*[1 nu 0; nu 1 0; 0 0 ((1-nu)/2)];

% Formulate the element stiffness matrix

Ke = t*(B'*D*B*A123);
