set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',14)

% Dimensions

B = 40;
a = 3;
b = 5;
c = 4;
D = a + b + c;

% Array initialization

tbot = (10:30)/1000;
sig0 = zeros(1,numel(tbot));
sig30 = zeros(1,numel(tbot));
sig90 = zeros(1,numel(tbot));

dsig0_by_da = zeros(1,numel(tbot));
dsig30_by_da = zeros(1,numel(tbot));

% Compute stress derivative d/da_m(sigma) for different ship heel angles

for i = 1:numel(tbot)
    [sig, A, h_NA, I_NA, I_CL, y_max, z_max] = combined_bending(tbot(i));
    
    sig0(i) = sig(0);
    sig30(i) = sig(pi/6);
    sig90(i) = sig(pi/2);
    
    zm = -h_NA;                                   % Location of member with respect to neutral axis
    zl = D-h_NA;                                   % Location of point at which stress is computed with respect to neutral axis
    ym = ((B/2)-a)/2;
    dsig0_by_da(i) = -sig0(i)*cos(0/6)*( (zm/(A*zl)) + ((zm^2)/I_NA) ) - sig90(i)*sin(0/6)*(((ym^2) + (((B/2)-a)^2)/12)/I_CL);                        % Stress derivative for  0 deg heel
    dsig30_by_da(i) = -sig0(i)*cos(pi/6)*(zm/(A*zl) + zm^2/I_NA)-sig90(i)*sin(pi/6)*((ym^2+ (B/2-a)^2/12)/I_CL);                      % Stress derivative for 30 deg heel (This is more complicated than simply using a formula from class)
    
end

% Choose one value of bottom thickness and compare actual stress variation
% with the linear approximation at this point. You can try changing ind and 
% see how the results look at a different value of tbot.

ind = 11;
am = tbot*(B - 2*a);

% Slope and intercept of linear approximation for 0 deg heel

m0 = dsig0_by_da(ind);
c0 = sig0(ind) - dsig0_by_da(ind)*am(ind);

% Slope and intercept of linear approximation for 30 deg heel

m30 = dsig30_by_da(ind);
c30 = sig30(ind) - dsig30_by_da(ind)*am(ind);

% Generate plot for visual comparison

plot(tbot*1000,sig0*1e-6,'r',tbot*1000,(m0*am + c0)*1e-6,'r--')
hold all
plot(tbot*1000,sig30*1e-6,'k',tbot*1000,(m30*am + c30)*1e-6,'k--')
plot(tbot(ind)*1000,sig0(ind)*1e-6,'bo',tbot(ind)*1000,sig30(ind)*1e-6,'bo')
grid on
title('Combined bending stress vs t_{bot}')
xlabel('t_{bot} in mm')
ylabel('\sigma in MPa')
legend('Heel  0^0','Heel  0^0 - Linear Approx.','Heel 30^0', 'Heel 30^0 - Linear Approx.')
legend boxoff
box off
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
