function Ke = cssr(node1, node2, node3, node4, Ey, nu, t)

x1 = node1(1); y1 = node1(2);
x2 = node2(1); y2 = node2(2);
x3 = node3(1); y3 = node3(2);
x4 = node4(1); y4 = node4(2);

if abs(x1 - x2) ~= abs(x3 - x4) || abs(y1 - y4) ~= abs(y3 - y2)
    error('The nodes do not form a rectangle')
end

a = abs(x1 - x2);
b = abs(y4 - y1);

ar = a/b;

A = (1 - nu)/2;
B = (1 + nu)/2;
C = (1 - (3*nu))/2;
D = (4 - (nu^2))/3;
E = (2 + (nu^2))/3;
A1 = A*ar; A2 = A/ar;
D1 = D*ar; D2 = D/ar;
E1 = E*ar; E2 = E/ar;

Ke = ((Ey*t)/(4*(1 - (nu^2))))*[A1+D2 B A1-D2 -C -A1-E2 -B -A1+E2 C; 
                               B A2+D1 C -A2+E1 -B -A2-E1 -C A2-D1;
                               A1-D2 C A1+D2 -B -A1+E2 -C -A1-E2 B;
                               -C -A2+E1 -B A2+D1 C A2-D1 B -A2-E1;
                               -A1-E2 -B -A1+E2 C A1+D2 B A1-D2 -C;
                               -B -A2-E1 -C A2-D1 B A2+D1 C -A2+E1;
                               -A1+E2 -C -A1-E2 B A1-D2 C A1+D2 -B;
                               C A2-D1 B -A2-E1 -C -A2+E1 -B A2+D1];

