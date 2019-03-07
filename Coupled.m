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
            B(n) = 3;
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
%from here nothing above really matters, the field is solved and its all
%good


C.q_0 = 1.60217653e-19;
C.m_0 = 9.10938215e-31;
C.kb = 1.3806504e-23;
C.T = 300;
frameWidth = 200e-9;
frameHeight = 100e-9;
nAtoms = 10;
Vth = sqrt(2*C.kb*C.T /(0.26*C.m_0));
dt = frameHeight/Vth/100;
Tstop = 1000*dt;
t = 0;
freepath = 0.2e-12;
Pscatter = 1 - exp(-dt/freepath);
AccelX = 2e9* Ey .* C.q_0 ./ C.m_0;
AccelY = 1e9* Ex .* C.q_0 ./ C.m_0;

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

%defining box boundaries
boxleft = X>0.8e-7;
boxright = X<1.2e-7;
box1bottom = Y>0.6e-7;
box2top = Y<0.4e-7;
boxedin = (boxleft&boxright&box1bottom)|(boxleft&boxright&box2top);
%removing particles from the boxes
while sum(boxedin)>0
X(boxedin) = rand*frameWidth;
Y(boxedin) = rand*frameHeight;
boxleft = X>0.8e-7;
boxright = X<1.2e-7;
box1bottom = Y>0.6e-7;
box2top = Y<0.4e-7;
boxedin = (boxleft&boxright&box1bottom)|(boxleft&boxright&box2top);
end
%all particles are out of the boxes to begin with

while t < Tstop
    
    Xindex = ceil(X * 1e9);
    Yindex = ceil(Y * 1e9);
    %scatter
    R = rand(1,nAtoms);
    VX(R<Pscatter) = Vth*randn(1);
    VY(R<Pscatter) = Vth*randn(1);
    %accelerate
    ain = sub2ind(size(AccelX),Xindex, Yindex);
    aiin = sub2ind(size(AccelY),Xindex, Yindex);
    VX = VX + AccelX(ain)*dt;
    VY = VY + AccelY(aiin)*dt;
    V = sqrt(VY.*VY+VX.*VX);
    
    %set next
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
    Ynext = Y + VY*dt;
    %set boundary for sides of boxes
    box1sides = (Ynext>0.6e-7)&(Xnext>0.78e-7)&(Xnext<1.22e-7);
    box2sides = (Ynext<0.4e-7)&(Xnext>0.78e-7)&(Xnext<1.22e-7);
    topbox = (Y<0.57e-7)&(Ynext>0.57e-7)&(Xnext>0.8e-7)&(Xnext<1.2e-7);
    bottombox = (Y>0.43e-7)&(Ynext<0.43e-7)&(Xnext>0.8e-7)&(Xnext<1.2e-7);
    VY(topbox|bottombox) = VY(topbox|bottombox) * -1;
    VX(topbox|bottombox|box1sides|box2sides) = VX(topbox|bottombox|box1sides|box2sides) * -1;
    
    %calculations for temperature
    Temperature(iteration) = 0.26*C.m_0*mean(V.^2)/4/C.kb;
    figure(1)
    xlim([0 frameWidth])
    ylim([0 frameHeight])
    plot([0.8e-7 0.8e-7],[0 0.4e-7], 'black')
    plot([1.2e-7 1.2e-7],[0 0.4e-7], 'black')
    plot([0.8e-7 1.2e-7],[0.4e-7 0.4e-7], 'black')
    plot([0.8e-7 0.8e-7],[1e-7 0.6e-7], 'black')
    plot([1.2e-7 1.2e-7],[1e-7 0.6e-7], 'black')
    plot([0.8e-7 1.2e-7],[0.6e-7 0.6e-7], 'black')
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
