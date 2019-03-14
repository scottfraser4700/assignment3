%% Finite Difference Solution of Electric Field
% This section prepares the electric field solution for use in the next
% section of the assignment. In the previous section, we assumed a constant
% field across the entirety of the simulation window. Here, we have a
% non-constant conductivity across the field of view. The different
% boundary conditions produce different results for the potential and
% electric field, both of which are plotted. The finite difference matrix method
% is used, in exactly the same way as in assignment 2.
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

 
%plotting potential
figure(4)
surf(Vmap)

title('Potential')
ylabel('X dimension (L=200)')
xlabel('Y dimension (W=100)')

%plotting electric field
[Ex,Ey] = gradient(-Vmap);
figure(5)
quiver(Ex,Ey)
title('Electric field')
ylabel('X dimension (L=200)')
xlabel('Y dimension (W=100)')
