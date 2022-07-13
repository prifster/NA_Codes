The given ship is a 175m containership. The frame of reference (co-ordinate system) is located at the mean water level and at midship. The positive x-axis is toward bow and positive y axis toward port and positive z-axis is upward. You will use 2D-strip theory for the calculation. 

 The ship has 25 stations  and the station numbering starts from the aft most station. 'x_loc' shows the  x-location of these staitons
 'offsetN' shows the offset of the sections. Each station is subdivided into 15 segments. The first 2 rows of 'offsetN' corrspond to the y and z locations of the aft most station i.e.station no 1. The next 2 rows (i.e. row 3 and 4) correspong to the y and z location of the 2nd station etc. 
For each station, you will start from left most (port side) point and follw anti-clock wise direction until you reach the starboard end. This will help you to be maintain consistancy with the solution during the estimation of unit normal to the surface.
Unit normals in y and z needs to be calculated for each line segment on the station.
Calcualte Ship volume, longitudinal and vertical center of bouyancy.
Sectional Froude krylov heave and pitch forces/moments need to be calculated. First you need to calculate the sectional Froude-Krylov forces  where Nz is the unit normal vecotor acting in z direction and  is the wave heading, m is the number of stations and n is the number of frequencies. Here the summation is over a particular section or station. dl is the length of the small segment along the length of the cross section. Once you have sectional forces, you can calculate the Froude-Krylov force acting on the whole ship for each given frequencies. You can use the command 'trapz' for integration along the length. There are 10 frequencies. So you will have 10 values of gloabl/final F-K force. The FK force will be in complex number. The size of the global heave force will be n x 1 , where n is the number of frequencies.
You can calculate the sectional Froude-Krylov moment using the formulation ,where x is longitudinal distance to the section measured from the frame of reference.Then You have to calculate the Froude-Krylov moment acting on the whole ship by integrating these sectional forces for each given frequencies. There are 10 frequney that I have given you. So you will have 10 values of final F-K moments for the whole ship.The FK moment will be in complex number. The size of the global pitch moment will be nx1 , where n is the number of frequencies
The diffraction forces are given.
The addedmass and damping forces are given
You need to calculate the heave and pitch restoring coefficients
Then solve the 2DoF equation of motion in frequency domaain


%% Code:
restoredefaultpath
savepath

hangle=180*pi/180; % WAVE HEADING ANGLE
load('offsetN.mat') % LOADING OFFSET - OFFSET START FROM AFT SECTION AND PROCEEDS TOWARD FORWARD SECTION.
load('x_loc') % LOADING STATION(SECTION) X LOCATION- FIRST STATION SHOWS THE AFT MOST STATION

load('MASS_PROPERTIES') % LOADING MASS PROPERTY
load('DIFFRACTION_FORCE') % LOADING HEAVE AND PITCH DIFFRACTION FORCES
load('ADDEDMASS') % LOADING HEAVE AND PITCH ADDEDMASS
load('DAMPING') % LOADING HEAVE AND PITCH DAMPING COEFFICIENTS
load('OMEGA') % LOADING WAVE FREQUENCY- VARIABLE NAME 'OMEGA'

Mass=MASS_PROPERTIES(1); % MASS OF THE SHIP
I5=MASS_PROPERTIES(2); % MASS MOMENT OF INERTIA OF THE SHIP FOR PITCH MOTION
lcg=MASS_PROPERTIES(3); % LCG OF THE SHIP WITH REFERENCE TO THE CO-ORDINATE SYSTEM
zcg=MASS_PROPERTIES(4); % VCG OF THE SHIP FROM THE MEAN WATER LEVEL WHERE THE CO-ORDINATE SYSTEM IS LOCATION

FD3=DIFFRACTION_FORCE(:,1); % HEAVE DIFFRACTION FORCE
FD5=DIFFRACTION_FORCE(:,2); % PITCH DIFFRACTION FORCE

A33=ADDEDMASS(:,1); % HEAVE ADDEDMASS
A35=ADDEDMASS(:,2); % HEAVE ADDEDMASS DUE TO PITCH MOTION
A53=ADDEDMASS(:,3); % PITCH ADDEDMASS DUE TO HEAEV MOTION
A55=ADDEDMASS(:,4); % PITCH ADDEDMASS DUE TO PITCH MOTION

B33=DAMPING(:,1); % HEAVE DAMPING COEFFICIENT
B35=DAMPING(:,2); % HEAVE DAMPING COEFFICIENT DUE TO PITCH MOTION
B53=DAMPING(:,3);% PITCH DAMPING DUE TO HEAEV MOTION
B55=DAMPING(:,4);% PITCH DAMPING DUE TO PITCH MOTION


lbp=175; % LBP
g=9.81;
density=1000; % DENSITY OF WATER

%% calculate sectional areas
for p=1:length(x_loc)
     areas(p) = trapz(offsetN(2*p-1,1:15), offsetN(2*p,1:15));
     fmzarea(p) = trapz(offsetN(2*p-1,:),((offsetN(2*p,:).^2)/2))/areas(p);

end

 
%% CCCCCC UNIT NORMAL CALCULATION CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

[row col]=size(offsetN);

      for p=1:length(x_loc)
                for l=1:col-1
                        
%                     Write your program for estimation of the unit normal  in y and z directions for each line segment in a station
                    y1=offsetN(2*p-1,l);
                    y2=offsetN(2*p-1,l+1);
                    z1=offsetN(2*p,l);
                    z2=offsetN(2*p,l+1);
                    Normaly(p,l)=(z1-z2)/sqrt((y2-y1)^2+(z2-z1)^2);
                    Normalz(p,l)=(y2-y1)/sqrt((y2-y1)^2+(z2-z1)^2);
                    midz(p,l)=(z1+z2)/2;
                    midy(p,l)=(y1+y2)/2;
                    dl(p,l)=sqrt((y2-y1)^2+(z2-z1)^2);  


                 end
      end
            
%% CCCCCCCCCCC VOLUME AND CENTRE OF BUOYANCY  OF THE SHIP CCCCCCCCCCCCCCCCCCCCCCCCCC

% Write your program to calculate the volume and the longitudinal and verticcal centre of bouyancy of the ship
for h = 4:20
    if h == 4 || h == 20
        s(h) = areas(h);
    elseif h == 6:2:18
        s(h) = 2*areas(h);
    else
        s(h) = 4*areas(h);
    end
end
volume = trapz(x_loc,areas);

COBX = trapz(x_loc,x_loc.*areas)/trapz(x_loc,areas);
COBZ = trapz(x_loc,fmzarea.*areas)/trapz(x_loc,areas);

%% CCCCCC FROUDE KRYLOV FORCE - FOR SOLVING FREQUENCY DOMAIN PROBLEM% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

fk3=zeros(length(OMEGA),length(x_loc));
fk5=zeros(length(OMEGA),length(x_loc));

for n=1:length(OMEGA)
     for p=1:length(x_loc)
        
       
       fk3(n,p)= density*g*sum(Normalz(p,:).*(exp((OMEGA(1,n).^2)/g*midz(p,:)).*exp(-i*((OMEGA(1,n).^2)/g*x_loc(p)*cos(hangle)-sin(hangle)*((OMEGA(1,n).^2)/g*midy(p,:))))).*dl(p,:));  % sectional Froude-Krylov heave force
       fk5(n,p)= -x_loc(p)*fk3(n,p);% sectional Froude-Krylov pitch moment
       
         end
    
 FK3(n,1)=trapz(x_loc,fk3(n,:)); % global F-K heave force. You can use the command trapz to calculate this force
 FK5(n,1)=trapz(x_loc,fk5(n,:));% global F-K pitch moment. You can use the command trapz to calculate this force
end



%% CCCCCC TOTAL EXCITING FORCE  - FROUDE-KRYLOV+DIFFRACTION FORCE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(

F3=(FD3+FK3); % THIS SHOULD BE A nX1 MATRIX, where n is the number of frequencies
F5=(FD5+FK5); % THIS SHOULD BE A nX1 MATRIX


% STIFFNESS MATRICES 

for p = 1:length(x_loc)
    y(p) = offsetN(2*p-1,1);
end
areaw = 2*trapz(x_loc,y);

C33 = density*g*areaw;% CALCULATE THE HEAVE-HEAVE RESTORING COEFFICIENT
C35 = -density*g*2*trapz(x_loc,y.*x_loc);% CALCUALTE THE HEAVE-PITCH RESTORING COEFFICIENT
C53 = C35;% CALCULATE THE PITCH-HEAVE RESTORING COEFFICIENT
C55 = density*g*volume*2*trapz(x_loc,(y.*(x_loc).^2)/volume);% CALCCUALTHE PITCH-PITCH RESTORING COEFFICIENT

 MASSMATRIX = [Mass -Mass*lcg;-Mass*lcg I5];% WRITE MASS MATRIX HERE . THIS WILL BE 2X2 MATRIX
 STIFFNESS = [C33,C35;C53,C55];% WRITE STIFFNESS MATRIX. THIS WILL BE 2X2 MATRIX



%% CCCCCC FREQUENCY DOMAIN SOLUTION OF EQUATION OF MOTION CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

% SOLVE THE EQUATION OF MOTION HERE BASED ON FREQUENCY DOMAIN SOLUTION
for m=1:length(OMEGA)
    M = [A33(m),A35(m);A53(m),A55(m)]+MASSMATRIX;
    B = [B33(m),B35(m);B53(m),B55(m)];
    F = [F3(m);F5(m)];
    req = inv([-OMEGA(m)^2*M + 1i*OMEGA(m)*B + STIFFNESS])*F;
    
    z3(m)=req(1);
    z5(m)=req(2);
end
z3= z3;% Heave RAO - THIS SHOULD HAVE VALUES FOR EACH FREQUENCY
z5= z5;% PITCH RAO- THIS SHOULD HAVE VALUES FOR EACH FREQUENCY

wl=g*(2*pi./OMEGA).^2/2/pi/lbp;
k0=OMEGA.^2/g;
figure;subplot(2,1,1),plot(wl,abs(z3)),xlim([ 0. 2.5])
subplot(2,1,2),plot(wl,abs(z5)./k0.'),xlim([ 0. 2.5])

