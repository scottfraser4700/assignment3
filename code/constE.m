%% Section 1: Constant Electric Field
% In this section, a constant electric field is applied in the X-direction
% of the simulation plane. Each electron experiences an equal force from
% the field, since the magnitude of the field is equal everywhere, and the
% electric force is F = eE, which is of course also constant for all
% involved electrons in the system. The acceleration for each particle is
% then A = eE/m, which is also constant for every particle. Using a value
% of 1V shows a very defined curving to the right, which leaves an electric
% field of 1e7V/m. The magnitude of the force experienced by each particle
% is then 1.602e-12 N, which in turn leads to an acceleration of
% 1.7588m/s^2.
% The current density component from electron drift is given by J = qnuE,
% where u is the mobility of electrons, E is the electric field, n is the
% number of free electrons, and q is the elementary charge. 
clear
C.q_0 = 1.60217653e-19;
C.m_0 = 9.10938215e-31;
C.kb = 1.3806504e-23;
C.T = 300;
frameWidth = 200e-9;
frameHeight = 100e-9;
nAtoms = 1000;
Vth = sqrt(2*C.kb*C.T /(0.26*C.m_0));
dt = frameHeight/Vth/100;
Tstop = 750*dt;
t = 0;
freepath = 0.2e-12;
Pscatter = 1 - exp(-dt/freepath);
Voltage = 1;
Efield = Voltage / frameWidth;
Force = Efield * C.q_0;
Accel = Force / C.m_0;
J = C.q_0*nAtoms*Efield;

%initializing vectors
Xnext = zeros(1,nAtoms);
Ynext = zeros(1,nAtoms);
VX = Vth * randn(1,nAtoms);
VY = Vth * randn(1,nAtoms);
V = sqrt(VY.*VY+VX.*VX);
X = frameWidth * rand(1, nAtoms);
Y = frameHeight * rand(1, nAtoms);
R = zeros(1, nAtoms);
Temperature = zeros(1, 100);
iteration = 1;


while t < Tstop
    %determines which particles scatter and performs calculations on them
    %to determine mean free path and time between collisions
    R = rand(1,nAtoms);
    VX(R<Pscatter) = Vth*randn(1);
    VY(R<Pscatter) = Vth*randn(1);
    VX = VX + Accel*dt;
    V = sqrt(VY.*VY+VX.*VX);
    
    
    Xnext = X + VX*dt;
    Ynext = Y + VY*dt;
    %X boundary conditions set
    right = Xnext>frameWidth;
    left = Xnext<0;
    Xnext(right) = Xnext(right)-frameWidth;
    Xnext(left) = Xnext(left) + frameWidth;
    %Y boundary conditions set
    top = Ynext > frameHeight;
    bottom = Ynext < 0;
    VY(top | bottom) = VY(top | bottom) * -1;
    %calculations for temperature
    Temperature(iteration) = 0.26*C.m_0*mean(V.^2)/4/C.kb;
    figure(1)
    xlim([0 frameWidth])
    ylim([0 frameHeight])
    hold on
    %plotting, but avoid plotting the full horizontal jump
    if abs(Xnext(1) - X(1)) < 2*abs(VX(1))*dt
        figure(1)
        plot([Xnext(1) X(1)], [Ynext(1) Y(1)], 'blue')
    end
    if abs(Xnext(2) - X(2)) < 2*abs(VX(2))*dt
        figure(1)
        plot([Xnext(2) X(2)], [Ynext(2) Y(2)], 'red')
    end
    if abs(Xnext(3) - X(3)) < 2*abs(VX(3))*dt
        figure(1)
        plot([Xnext(3) X(3)], [Ynext(3) Y(3)], 'green')
    end
    if abs(Xnext(4) - X(4)) < 2*abs(VX(4))*dt
        figure(1)
        plot([Xnext(4) X(4)], [Ynext(4) Y(4)], 'black')
    end
    
    %updating positions, and advancing time a step forward so the while
    %loop works
    X = Xnext;
    Y = Ynext;
    t = t+dt;
    iteration = iteration + 1;
    pause(0.0001);
end


%electron density map
figure(2)
EDM = hist3([X',Y'],[30,30]);
pcolor(EDM')
title('Electron Density Map')
view(2)

%temperature map
xLim = linspace(0,frameWidth,100);
yLim = linspace(0,frameHeight,100);
xTempReg = discretize(X,xLim);
yTempReg = discretize(Y,yLim);
for q=1:1:100
    for w=1:1:100
        %Temperature contained in defined region
        tempReg = (q == xTempReg) & (w == yTempReg); 
        
        %Total velocities in region
        vxTot=sum(VX(tempReg));
        vyTot=sum(VY(tempReg));
        vTot = sqrt((vxTot)^2+(vyTot)^2);
        
        %Calculate Temperature
        tempMap(q,w) = C.m_0*0.26*(vTot)^2/(2*C.kb);
    end
end
figure(3)
surf(tempMap)
view(2)
title('Temperature Map')