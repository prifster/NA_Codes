set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',1.5)

B = 30;                     % Breadth of Ship
D = 15;                     % Depth of Ship
G = 80 * 1e9;               % Shear Modulus of Steel
L = 200;                    % Length of Ship 
Cb = 0.6;                   % Block coefficient of the Ship
th_des = 5;                 % Permissible twist angle in degrees

H = 8.13 - ((250 - (0.7*L))/(125))^3;
Cw = 0.165 + (0.95*Cb);
Cm = 0.45*(B^2)*(Cw^2);
Cq = 5*D*Cb;
Km = 2;
Kq = 0;
d = ((D^2)/(2*D + B)) - (0.6*D);
 
% Calculate the twisting moment generated at the midship section for
% ship heading bearing an angle 60 degrees to incident wave direction

Mx = 250*H*L*((Km*Cm) + (Kq*Cq*d));

% Calculate the thickness such that the maximum twist over the length of
% the ship is less than 5 degrees

t = ((Mx*3*180*L)/(G*(2*D + B)*th_des*pi))^(1/3);

fprintf('The required thickness of plate is %.2f mm !!\n',t*1000);

% We can approximate that the amount of material used in building the ship
% is roughly proportional to the area of structural midship cross section.
% Calculate the area of cross section of the midship Acs

Acs = ((2*D) + B)*t;

fprintf('The area of cross-section is %.2f m^2\n',Acs);

% Now we want to calculate the warping of this section as it twists and
% visualize the twist

% Arc length

s = 0:0.1:B/2+D;

% Calculate the warping function about point D as given in the notes

wD = zeros(size(s));
wD(s > B/2) = (B/2)*(s(s>B/2) - B/2);

% Calculate the integrand of integral calculating the location of center of
% twist

integrand = zeros(size(s));
integrand(s > B/2) = (B/2)*t.*wD(s > B/2);

% Calculate the second moment of area of cross section about eta axis

zeta = s;
zeta(s > B/2) = B/2;


I_eta = (((B^2)*D*t)/4) + (((B^3)*t)/24);;

% Calculate the location of center of twist along eta-axis 

eta0 = -(1/I_eta)*(((B^2)*(D^2)*t)/8);;

% Calculate St. Venant's constant

J = (1/3)*((2*D)+B)*(t^3);

% Calculate Rate of Twist

dtheta_dx = Mx/(G*J);

% Calculate Twist over the length of the ship L

theta = (Mx*L)/(G*J);

% Calculate Maximum shear stress in the section

tau =  (Mx*t)/J;

fprintf('Maximum shear stress in the section is %.2f MPa\n',tau*1e-6);

% You need not edit any lines below this
%--------------------------------------------------------------------------

% Visualization

% Original section

x_orig = zeros(size(s));
y_orig = s;
y_orig(s > B/2) = B/2;
z_orig = s - B/2 - eta0;
z_orig(s < B/2) = 0 - eta0;

xf_orig = [x_orig x_orig];
yf_orig = [fliplr(-y_orig) y_orig];
zf_orig = [fliplr(z_orig) z_orig];

% Calculate the normalized warping function

wn = -wD + eta0*zeta;

% New Section

x1_new = x_orig + dtheta_dx*wn;
y1_new = y_orig - theta*z_orig;
z1_new = z_orig + theta*y_orig;

x2_new = x_orig - dtheta_dx*wn;
y2_new = -y_orig - theta*z_orig;
z2_new = z_orig - theta*y_orig;

xf_new = [fliplr(x2_new) x1_new];
yf_new = [fliplr(y2_new) y1_new];
zf_new = [fliplr(z2_new) z1_new];

figure
plot(yf_orig,zf_orig)
hold all
plot(yf_new,zf_new)
plot(0,0,'rp')
title('Cross Section before and after warping')
ylabel('z')
xlabel('y')

figure
plot3(xf_orig,yf_orig,zf_orig)
hold all
plot3(xf_new,yf_new,zf_new)
plot3(0,0,0,'rp')
xlim([-1 1])
title('Cross Section before and after warping')
ylabel('y')
xlabel('x')
zlabel('z')
