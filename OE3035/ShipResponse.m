In this assignment, you will write a program to generate P-M spectrum and will calcualte the heave and pitch responses of a 175m containership in head seas. The ship is moving a Froude-number of 0.15. Please follow the following instructions for calculations carefully
Calculate the wave spectrum using PIERSON-MOSKOWITZ (PM) spectrum for the following wave parameters : Hs=3m , Tp=12 sec; The frequency range () should be chosen between 0.15 and 1 rad/sec. 
Calculate the encountering wave spectrum
Calcualte the encountering wave frequency.
The heave and pitch (rad/m) RAO are given as z3 and z5 for a set of encountering frequency. Their corresponding encoutnering frequency are given as the variale 'OMEGA'. As you can see that frequency values of RAO and wave spectrum are different. So you will have to interpolate heave and pitch RAO values for the corresponding PM spectra frequencies. You can use the function 'interp1' for interpolation. Calculate interpolated values of heave and pitch RAO. The interpolated values should be called by the variable 'z3i' and 'z5i'.
Calculate the significant wave height and zero upcrossing period of the incident waves
Calculate the significant height and zero upcrossing period for the heave responses
Calculate the significant height and zero upcrossing period for the pitch responses

%% Code:
restoredefaultpath
savepath

load z3.mat
load z5.mat
load OMEGA

Fr = 0.15;
g = 9.81;
L = 175;

v= Fr*sqrt(g*L);% calcualte velocity of the ship

dw=2*pi/1000; %incremental frequency (dw)
w=0.15:dw:1;  % frequency range for wave spectrum

Tp = 12;
Hs = 3;

% Calcualte Sw

for og = 1:numel(w)
    wp = (2*pi)/Tp;
    wm = w(og)/wp;
    Sw(og) = (5/16)*(Hs^2/wp)*wm^(-5)*exp((-5/4)*wm^(-4));
end

% Calcualte Swe

for i = 1:numel(w)
    Swe(i) = Sw(i)/(1- (2*w(i)*v*cos(pi)/g));
end


% Calcualte we

for d = 1:numel(w)
    we(d) = w(d)*(1 - (v*(w(d)/g)*cos(pi)));
end

figure
plot(we,Swe,'-r',w,Sw,'-b')

Hs_wave= 4*sqrt(trapz(w,Swe)); % calculate significant wave height



Tz_wave=12; % calcualte zero upcrossing wave period

z3i=interp1(OMEGA,abs(z3),we);% interpolate heave RAO values
z5i=interp1(OMEGA,abs(z5),we);% interpolate pitch RAO values


Sheave=z3i.^2.*Swe;% calcualte heave encountering spectrum
Spitch=z5i.^2.*Swe;% calculate pitch encoutnering spectrum



Hs_heave=4*sqrt(trapz(we,Sheave));% calcualte significant heave height
Tz_heave=2*pi*sqrt(trapz(we,Sheave)/trapz(we,we.^2.*Sheave)); % calcualte heave zero upcrossing period

Hs_pitch=4*sqrt(trapz(we,Spitch));% calculate significant pitch angle
Tz_pitch=2*pi*sqrt(trapz(we,Spitch)/trapz(we,we.^2.*Spitch)); % calcualte zero upcrossing period

figure;subplot(2,1,1),plot(we,Sheave)
subplot(2,1,2),plot(we,Spitch)

