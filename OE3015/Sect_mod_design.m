% Script to optimize the scantlings of a midship section to minimize the
% weight while ensuring adequate section modulus to withstand bending 
% stresses due to stillwater and wave induced loads
%
% Developed by: Abhilash Somayajula 
% Course OE3015 - Ship Structural Analysis

syms t1 t2 real

tol = 1e-4;         % Error tolerance
L = 100;            % Length of vessel
B = 20;             % Breadth of vessel
D = 10;             % Depth of vessel
T = 6;              % Draft of vessel
r = 1;              % Bilge radius
Cb = 0.8;           % Block coefficient
fa = 1e6;           % Maximum still water distributed load

% Convert to symbolic variables for section modulus computation

B = sym(B);
D = sym(D);
r = sym(r);

% Initialize arrays

area = sym(zeros(1,9));
height = sym(zeros(1,9));
inertia_local = sym(zeros(1,9));

% Calculate the contribution of each element to the area of cross-section
% If any expression involves pi, use "pi" in the expression instead of "3.14"

area(1) = t2*B;
area(2) = t1*(D-r);
area(3) = t1*(D-r);
area(4) = t1*(D-r);
area(5) = t1*(D-r);
area(6) = t2*B;
area(7) = (t2*pi*r)/2;
area(8) = t2*(B-(2*r));
area(9) = (t2*pi*r)/2;

% Calculate the height of centroid of each element's area from baseline
% Remember the height will not depend on the thickness t1 or t2 
% If any expression involves pi, use "pi" in the expression instead of "3.14"

height(1) = D;
height(2) = (D+r)/2;
height(3) = (D+r)/2;
height(4) = (D+r)/2;
height(5) = (D+r)/2;
height(6) = r;
height(7) = (pi - 2)/pi *r;
height(8) = 0;
height(9) = (pi - 2)/pi *r;

% Calculate local area moment of inertia about the element's neutral axis.
% Note the we neglect the local area moment of inertia of the deck (element
% 1), tank top (element 6) and bottom shell (element 8).

inertia_local(2) = (1/12)*area(2)*(D-r)^2;
inertia_local(3) = (1/12)*area(3)*(D-r)^2;
inertia_local(4) = (1/12)*area(4)*(D-r)^2;
inertia_local(5) = (1/12)*area(5)*(D-r)^2;
inertia_local(7) = ((1/2) - (4/pi^2))*area(7)*(r^2);
inertia_local(9) = ((1/2) - (4/pi^2))*area(9)*(r^2);

% Calculate the total area and its first and second moments about baseline

area_total = sum(area);                                 % Total area of the cross-section
moment_1st = sum(area.*height);                         % First moment of area of cross-section about the baseline
moment_2nd = sum(area.*height.^2 + inertia_local);      % Second moment of area of cross section about the baseline


% Distance to neutral axis from baseline

h_NA = moment_1st/area_total;

% Area moment of inertia of cross section about neutral axis

I = moment_2nd - area_total*h_NA^2;

% Calculate section modulus for deck

Z = I/(D-h_NA);
Zfun = matlabFunction(Z);

% Calculate derivatives of section modulus

dZ_by_dt1 = matlabFunction(diff(Z,t1));
dZ_by_dt2 = matlabFunction(diff(Z,t2));

% Convert B, D and r into double variables

B = double(B);
D = double(D);
r = double(r);

% Wave parameter C

if L>=90 && L<300
    C = 10.75 - ((300 - L)/100)^(1.5);
elseif L>=300 && L<350
    C = 10.75;
elseif L>=350 && L<=500
    C = 10.75 - ((300 - L)/150)^(1.5);
end

% Distributed Stillwater Load

f = @(x) fa*cos(2*pi*x/L);

% Compute Stillwater Shear Force

Qs = @(x) ((fa*L)/(2*pi))*sin(2*pi*x/L);

% Compute Stillwater Bending Moment

Ms = @(x) -((fa*(L^2))/(2*pi)^2)*(cos(2*pi*x/L) - 1);

% Compute Stillwater Bending Moment at Midship

Ms_MS = Ms(L/2);

% Compute Wave Induced Bending Moment in Sagging at Midship

Mw_sag_MS = -110*1*C*(L^2)*B*(Cb + 0.7);

% Compute Wave Induced Bending Moment in Hogging at Midship

Mw_hog_MS = 190*1*C*(L^2)*B*Cb;

% Compute total bending moment in the Sagging and Hogging at Midship

Mt_sag_MS = Ms_MS + Mw_sag_MS;
Mt_hog_MS = Ms_MS + Mw_hog_MS;

% Maximum bending moment at midship

Mt_MS = max(abs([Mt_sag_MS Mt_hog_MS]));

%##########################################################################
% OPTIMIZATION - SEQUENTIAL LINEAR PROGRAMMING
%##########################################################################

% Cost function - Weight of a unit length of a beam with the cross-section 
% given by the midship cross-section defined above

f_cost = [double(diff(area_total,t1)); double(diff(area_total,t2));]/1000;

% Bounds and initial conditions for optimization of design variables

LB = [10; 10;];             % Lower bounds for design variables in mm
UB = [100; 100;];           % Upper bounds for design variables in mm

x0 = [30; 30;];             % Initial solution in mm
xn = x0;                    % Initial solution specification
fval = f_cost'*xn;          % Value of objective function at initial condition

fprintf('Cost evaluated at initial condition t1 = %.2f mm and t2 = %.2f mm is %.4f\n',x0(1),x0(2),fval)

% Iterate until convergence
while true
    
    % Limit the update of design variables to 5 % of their current value
    
    Dx = 0.05*xn;
    LBn = xn - Dx;
    UBn = xn + Dx;
    
    LBn(LBn < LB) = LB(LBn < LB);
    UBn(UBn > UB) = UB(UBn > UB);
    
    % Linearized constraint formulation
    A = zeros(2,2);
    b = zeros(2,1);
    
    ta = xn(1)/1000; tb = xn(2)/1000;
    
    A(1,1) = -dZ_by_dt1(ta,tb);
    A(1,2) = -dZ_by_dt2(ta,tb);
    A(2,1) = -dZ_by_dt1(ta,tb);
    A(2,2) = -dZ_by_dt2(ta,tb);
    
    %----------------------------------------------------------------------
    % Write your code in the space below to compute the minimum section 
    % modulus from ABS rules discussed in class. 
    %
    % The second part is to compute the vector b to capture the constraints. 
    % For this you may use the functions Zfun(t1,t2), dZ_by_Dt1(t1,t2) and 
    % dZ_by_dt2(t1,t2) defined above
    %----------------------------------------------------------------------
    
    Zmin1 = (Mt_MS)/(190*(10^6));
    Zmin2 = 0.9*C*(L^2)*B*(Cb+0.7)*(10^(-6));
    
    b(1) = Zfun(ta,tb) - Zmin1 - dZ_by_dt1(ta,tb)*ta - dZ_by_dt2(ta,tb)*tb;
    b(2) = Zfun(ta,tb) - Zmin2 - dZ_by_dt2(ta,tb)*tb - dZ_by_dt1(ta,tb)*ta;
    
    %----------------------------------------------------------------------
    
    b = b*1000;
    
    % Optimize using linear programming
    
    [xn1,fvaln,exitflag] = linprog(f_cost,A,b,[],[],LBn,UBn);
    
    % Check if optimization using linear programming succeeded
    
    if exitflag ~= 1
        error('Optimization failed. Recheck the problem parameters')
    end
    
    % Update the design variables
    
    xn = xn1;
    
    fprintf('Cost evaluated at new t1 = %.2f mm and new t2 = %.2f mm is %.4f\n',xn(1),xn(2),fvaln)
    
    % Check if tolerance is met for stopping the optimization routine
    
    if abs(fval - fvaln) < tol
        break
    end    
    
    fval = fvaln;
    
    
end

fprintf('The optimum solution is t1 = %.2f mm and t2 = %.2f mm\n',xn(1),xn(2));
