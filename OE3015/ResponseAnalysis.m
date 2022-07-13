Consider a ship midship section with the neutral axis  above the keel and with the deck located at  above the keel. Given that the second moment of area of cross-section  and Young's modulus of steel is  
What is the bending moment required to result for a deck compressive stress ?
What is resulting curvature of the hull girder?
  
%% Code:
D = 12;
I = 30;
z_NA = 5;
E = 200000;
sigmaD = 150;

M = (sigmaD*I*(10^6))/7
Curvature = M/(E*I*(10^6));

R = 1/Curvature;

fprintf('The required bending moment is %.2f Nm\n',M);
fprintf('The curvature is %.6f/m and radius is %.2f m\n',Curvature,R);
