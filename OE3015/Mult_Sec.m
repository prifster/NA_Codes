set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',1.5)

B = 30;                     % Breadth of Ship
D = 15;                     % Depth of Ship
a = 5;                      % Geometric Parameter "a"
G = 80 * 1e9;               % Shear Modulus of Steel
L = 200;                    % Length of Ship 
Cb = 0.6;                   % Block coefficient of the Ship
th_des = 5;                 % Permissible twist angle in degrees
t= 3/1000;                  % Plate thickness

% Calculate the twisting moment generated at the midship section for
% ship heading bearing an angle 60 degrees to incident wave direction

Ftm = 2;
Cw = 0.165 + 0.95 * Cb; 
Cm = 0.45 * B^2 * Cw^2;
H = 8.13 - ((250 - 0.7 * L)/125)^3;
Mx = 250 * H * L * (Ftm * Cm ); 

% Formulate the C matrix

C_bar = zeros(8);

C_bar(1,1) = (D + a)/t;
C_bar(1,2) = -a/t;

C_bar(2,1) = -a/t;
C_bar(2,2) = (D + a)/t;
C_bar(2,3) = -a/t;

C_bar(3,2) = -a/t;
C_bar(3,3) = (4*a)/t;
C_bar(3,4) = -a/t;

C_bar(4,3) = -a/t;
C_bar(4,4) = B/t;
C_bar(4,5) = -a/t;

C_bar(5,4) = -a/t;
C_bar(5,5) = B/t;
C_bar(5,6) = -a/t;

C_bar(6,5) = -a/t;
C_bar(6,6) = (4*a)/t;
C_bar(6,7) = -a/t;

C_bar(7,6) = -a/t;
C_bar(7,7) = (D + a)/t;
C_bar(7,8) = -a/t;

C_bar(8,7) = -a/t;
C_bar(8,8) = (D + a)/t;


% Fromulate Area Vector

Area = zeros(8,1);

Area(1) = ((D - a)*a)/2;
Area(2) = ((D - a)*a)/2;
Area(3) = a^2;
Area(4) = ((B/2) - a)*a;
Area(5) = ((B/2) - a)*a;
Area(6) = a^2;
Area(7) = ((D - a)*a)/2;
Area(8) = ((D - a)*a)/2;

% Calculate normalized shear flow

q_norm = 2*G*(C_bar\Area);

% We can approximate that the amount of material used in building the ship
% is roughly proportional to the area of structural midship cross section.
% Calculate the area of cross section of the midship Acs

Acs = (2*t*(B + (2*D) -a)) + (7*a*t);

fprintf('The thickness of plate is %.2f mm !!\n',t*1000);
fprintf('The area of cross-section is %.2f m^2\n',Acs);


% Rate of Twist

dtheta_dx = Mx/(2*Area'*q_norm);

% Twist over a length dL

theta = dtheta_dx*L;

% Maximum shear stress in the section

tau = max(q_norm*dtheta_dx)/t;

fprintf('Maximum shear stress in the section is %.2f MPa\n',tau*1e-6);
