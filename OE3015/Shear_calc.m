function out = shear_stress_fcn(tss,inp)

t = inp.t;
B = inp.B;
a = inp.a;
b = inp.b;
c = inp.c;
D = inp.D;
Q = inp.Q;

% Compute the location of neutral axis and second moment of inertia of the
% cross section using previously defined functions

out = combined_bending(tss,inp);
h_NA = out.h_NA;
I_NA = out.I_NA;

% Calculate the first moment upto the section at which shear stress or
% shear flow need to be computed.

g = (a+b+c)-h_NA;
m_AB = @(s1) t*(g)*s1

m_CD = @(s2) t*((b+c)-h_NA)*s2

m_BD = @(s3) (g)*t*(B/4) + (g*s3-(s3.*s3)/2)*tss 

m_DE = @(s3) (g)*t*(B/4) + (g*s3-(s3.*s3)/2)*tss + t*((b+c)-h_NA)*(B/4)

m_EF = @(s4)  g*t*B/4 + (g*(a+b)-((a+b)^2)/2)*tss + t*((b+c)-h_NA)*(B/4) - ((h_NA-c)*s4+(c/sqrt(a^2+c^2))*((s4.*s4)/2))*tss

m_GF = @(s5) -t*(h_NA)*s5

% You need not edit any of the following code.

s1 = linspace(0,B/4);
q1 = Q/I_NA * m_AB(s1) * 1e-6;

s2 = s1;
q2 = Q/I_NA * m_CD(s2) * 1e-6;

s3a = linspace(0,a);
s3b = linspace(a,a+b);
s3 = [s3a s3b];

q3a = Q/I_NA * m_BD(s3a);
q3b = Q/I_NA * m_DE(s3b);
q3 = [q3a q3b] * 1e-6;

s4 = linspace(0,sqrt(a^2 + c^2));
q4 = Q/I_NA * m_EF(s4) * 1e-6;

s5 = linspace(0,B/2 - a);
q5 = Q/I_NA * m_GF(s5) * 1e-6;

out = struct;
out.s1 = s1; out.s2 = s2; out.s3 = s3; out.s4 = s4; out.s5 = s5;
out.q1 = q1; out.q2 = q2; out.q3 = q3; out.q4 = q4; out.q5 = q5;

out.m_AB = m_AB; out.m_CD = m_CD; out.m_BD = m_BD; out.m_DE = m_DE; out.m_EF = m_EF; out.m_GF = m_GF;

end
