function out = combined_bending(tbot,inp)
%--------------------------------------------------------------------------
% This function calculates the properties of the cross-section of a ship
% including the sectional area A, second moments of inertia about both
% baseline I_NA and centerline I_CL, height of the neutral axis above the
% baseline h_NA.
%
% Developed by: Abhilash Somayajula
% Date: 25-Aug-2019
%--------------------------------------------------------------------------

if nargin == 2
    t = inp.t;
    B = inp.B;
    a = inp.a;
    b = inp.b;
    c = inp.c;
    D = inp.D;
    M = inp.M;
else
    t = 10/1000;
    B = 40;
    a = 3;
    b = 5;
    c = 4;
    D = a + b + c;
    M = 800 * 1e6;
end

% Array initialization

area = zeros(1,5);
height_base = zeros(1,5);
width_center = zeros(1,5);
local_inertia_horiz = zeros(1,5);
local_inertia_vert = zeros(1,5);

% Area of each element

area(1) = (B/4)*t;
area(2) = (B/4)*t;
area(3) = (a+b)*t;
area(4) = sqrt((a^2)+(c^2))*t;
area(5) = ((B/2)-a)*tbot;

% Height of each element's centroid from the baseline

height_base(1) = a+b+c;
height_base(2) = b+c;
height_base(3) = ((a+b)/2)+c;
height_base(4) = c/2;
height_base(5) = 0;

% Local second moment of inertia of each element about the horizontal axis
% passing through the centroid of the element.

local_inertia_horiz(1) = 0;
local_inertia_horiz(2) = 0;
local_inertia_horiz(3) = (area(3)*((a+b)^2))/12;
local_inertia_horiz(4) = (area(4)*(c^2))/12;
local_inertia_horiz(5) = 0;

% Width of each element's centroid from the centerline

width_center(1) = 3*B/8;
width_center(2) = 3*B/8;
width_center(3) = B/2;
width_center(4) = (B-a)/2;
width_center(5) = ((B/2)-a)/2;

% Local second moment of inertia of each element about the vertical axis
% passing through the centroid of the element

local_inertia_vert(1) = (area(1)*((B/4)^2))/12;
local_inertia_vert(2) = (area(2)*((B/4)^2))/12;
local_inertia_vert(3) = 0;
local_inertia_vert(4) = (area(4)*(a^2))/12;
local_inertia_vert(5) = (area(5)*((B/2)-a)^2)/12;

% First and second moments of area about the base and centerline

moment_1st_base = area .* height_base;
moment_2nd_base = area .* height_base.^2 + local_inertia_horiz;
moment_2nd_center = area .* width_center.^2 + local_inertia_vert;

% Location of the neutral axis relative to the baseline

h_NA = (sum(moment_1st_base)/sum(area));

% Second moment of inertia of the cross-section about the neutral axis

I_NA = 2*sum(moment_2nd_base-area*h_NA.^2);

% Second moment of inertia of the cross section about the centerline

I_CL = 2*sum(moment_2nd_center);

% Total area of the cross section

A = 2*sum(area);

% Distance of the fiber that would experience maximum stress within the
% section

z_max = height_base(1)-h_NA;
y_max = B/2;

% Compute stress due to combined loading as a function of heel angle th

sig = @(th) M*cos(th)/I_NA*z_max + M*sin(th)/I_CL*y_max

% Output

out = struct;
out.sig = sig;
out.A = A;
out.h_NA = h_NA;
out.I_NA = I_NA;
out.I_CL = I_CL;
out.y_max = y_max;
out.z_max = z_max;
out.area = area;
out.height_base = height_base;
out.width_center = width_center;
out.local_inertia_horiz = local_inertia_horiz;
out.local_inertia_vert = local_inertia_vert;
