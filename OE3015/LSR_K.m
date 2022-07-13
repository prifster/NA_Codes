function Ke = lsr(node1, node2, node3, node4, E, nu, t)

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

Keps = ((E*t)/(12*(1 - (nu^2))))*[4/ar 3*nu (-4/ar) 3*nu (-2/ar) -3*nu 2/ar -3*nu;
     3*nu 4*ar -3*nu 2*ar -3*nu -2*ar 3*nu -4*ar; 
     -4/ar -3*nu 4/ar -3*nu 2/ar 3*nu -2/ar 3*nu; 
     3*nu 2*ar -3*nu 4*ar -3*nu -4*ar 3*nu -2*ar; 
     -2/ar -3*nu 2/ar -3*nu 4/ar 3*nu -4/ar 3*nu; 
     -3*nu -2*ar 3*nu -4*ar 3*nu 4*ar -3*nu 2*ar; 
     2/ar 3*nu -2/ar 3*nu -4/ar -3*nu 4/ar -3*nu; 
     -3*nu -4*ar 3*nu -2*ar 3*nu 2*ar -3*nu 4*ar];

Kgam = ((E*t)/(24*(1 + nu)))*[4*ar 3 2*ar -3 -2*ar -3 -4*ar 3;
     3 4/ar 3 -4/ar -3 -2/ar -3 2/ar; 
     2*ar 3 4*ar -3 -4*ar -3 -2*ar 3;
     -3 -4/ar -3 4/ar 3 2/ar 3 -2/ar;
     -2*ar -3 -4*ar 3 4*ar 3 2*ar -3;
     -3 -2/ar -3 2/ar 3 4/ar 3 -4/ar;
     -4*ar -3 -2*ar 3 2*ar 3 4*ar -3; 
      3 2/ar 3 -2/ar -3 -4/ar -3 4/ar];
Ke = Keps + Kgam;
