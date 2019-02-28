%% Section 3 - Coupled Simulations
%first solve the potential and field, then have the particles interact with
%it
nx = 200;
ny = 100;
G = sparse(nx*ny, nx*ny);
B = zeros(1,nx*ny);
%setting conductivity map for 3rd part
cMap = ones(nx,ny);
 for q = 1:nx
     for w = 1:ny
         if ((q<(0.6*nx)&&q>(0.4*nx)&&w>(0.6*ny)) || (q<(0.6*nx)&&q>(0.4*nx)&&w<(0.4*ny)))
             cMap(q,w) = 1e-2;
         end
     end
 end

for i = 1:nx
    for j = 1:ny
        %Setting up the G-Matrix
        n = j + (i-1)*ny;
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,n) = 1;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G(n,n) = -(rxm + rxp + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        elseif j == ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            G(n,nxm) = rxm;
            G(n,nym) = rym;
            G(n,n) = -(rxm + rxp + rym + ryp);
        end
    end
end
%solving for matrix of potentials
V = G\B';
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end
[Ex,Ey] = gradient(-Vmap);


%tighten up my boundary conditions on the bottleneck so there is no
%penetration
%also, in order to use the EX,EY, and VX,VY vectors together, I need to
%index things in a different way. each VX and VY corresponds to a
%position in the X,Y vectors. at each step, I need to use the find which E
%vectors I'm accelerating each particle with.
C.q_0 = 1.60217653e-19;
C.m_0 = 9.10938215e-31;
C.kb = 1.3806504e-23;
C.T = 300;
frameWidth = 200e-9;
frameHeight = 100e-9;
nAtoms = 1000;
Vth = sqrt(2*C.kb*C.T /(0.26*C.m_0));
dt = frameHeight/Vth/100;
Tstop = 1000*dt;
t = 0;
freepath = 0.2e-12;
Pscatter = 1 - exp(-dt/freepath);
AccelX = Ex .* C.q_0 ./ C.m_0;
AccelY = Ey .* C.q_0 ./ C.m_0;

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
    Xindex = ceil(X * 1e9);
    Yindex = ceil(Y * 1e9);
    R = rand(1,nAtoms);
    VX(R<Pscatter) = Vth*randn(1);
    VY(R<Pscatter) = Vth*randn(1);
    VX = VX + AccelX*dt;
    VY = VY + AccelY*dt;
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
