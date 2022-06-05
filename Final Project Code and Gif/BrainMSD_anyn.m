% Final project
% Lily Clay
% interconnected mass spring dampers
clear;clc;close;

r = 4; %rows/cols of nodes
n = r*r; %number of nodes


%initial pos
x = zeros(n,2);

j = 0; i = r-1;
for k = 1:n
    x(k,1) = j; x(k,2) = i;
    j = j+1;
    if j == r
        j = 0; i = i - 1;
    end
end


%initial dof vector
q0 = zeros(2*n,1);
for i = 1:n
    q0(2*i-1) = x(i,1);
    q0(2*i) = x(i,2);
end

%radii
r0 = 0.01; %rod
s0 = 0.01; %sphere

%mass
rho = 3000;
M = zeros(2*n);
for k = 1:2*n
    M(k,k) = 4/3*pi*s0^3*rho;
end

%damper
b = zeros(2*n,1);
for k = 1:2*n
    b(k) = 10;
end

%viscosity
visc = 50;
C = zeros(2*n);
for k = 1:2*n
    C(k,k) = 6*pi*s0*visc;
end

%external force
P = zeros(2*n,1);
P(2*r-1) = -30;
P(4*r-1) = -30;

%gravity
g = 9.81;
W = zeros(2*n,1);
for k = 1:n
    W(2*k) = g;
end


%spring quantities
Y = 1e5;
EA = Y*pi*r0^2;

%time
Ttot = 5;
dt = 0.01;
nsteps = Ttot/dt;

tol = 1e-3;

%initial vel
u = zeros(size(q0));

%calculation
for i = 2:nsteps

    q = q0;
    %Newton Raphson
    err = 10*tol;
    while err > tol
        %inertia
        f = M/dt * ((q-q0)/dt-u);
        J = M/dt^2;
    
        %elastic forces
        %linear spring 1 between node k and k+1
        for k = 1:n-1
            if mod(k,r) == 0
                dF = zeros(6,1); dJ = zeros(6);
            else
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
             
            l_k = 1;
            dF = gradEs(xk,yk,xkp1,ykp1,l_k,EA);
            dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EA);
            ind = 2*k-1:2*k+2;
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
            end
        end

        %linear spring 1 between vertical nodes
        for k = 1:n-r

            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*(k+r)-1);
            ykp1 = q(2*(k+r));
     
            l_k = 1;
            dF = gradEs(xk,yk,xkp1,ykp1,l_k,EA);
            dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EA);
            ind = [2*k-1 2*k 2*(k+r)-1 2*(k+r)];
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
            
        end

        %linear spring diagonals
        for k = 1:n-r
            if mod(k,r) == 0
                dF = zeros(6,1); dJ = zeros(6);
            else
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*(k+r+1)-1);
                ykp1 = q(2*(k+r+1));
                
                l_k = sqrt(2);
                dF = gradEs(xk,yk,xkp1,ykp1,l_k,EA);
                dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EA);
                ind = [2*k-1 2*k 2*(k+r+1)-1 2*(k+r+1)];
                f(ind) = f(ind) + dF;
                J(ind,ind) = J(ind,ind) + dJ;
            end
        end

        %damping between nodes
        fd = zeros(2*n,1);
        for j = 1:r
            for k = 2*r*j-2*r+3:2*r*j-2
                fd(k) = (b(k)+b(k+2))*(q(k)-q0(k))/dt - ...
                b(k)*(q(k-2)-q0(k-2))/dt - ...
                b(k+2)*(q(k+2)-q0(k+2))/dt;
            end
        end

        Jd = zeros(2*n);

        for k = 1:2*n-2
        for j = 1:2*n
            bi = b(k);

            if k == 2*n-1
                bip1 = b(1);
            elseif k == 2*n
                bip1 = b(2);
            else
                bip1 = b(k+2);
            end
    
            if k == j
                Jd(k,j) = (bi+bip1)/dt;
            elseif j == k-1
                Jd(k,j) = bi/dt;
            elseif j == k+1
                Jd(k,j) = bip1/dt;
            else
                Jd(k,j) = 0;
            end
        end
        end

        f = f + fd;
        J = J + Jd;

%         %damping between nodes vertical
%         fd = zeros(2*n,1);
%         for j = 2:2
%             for k = 2*r*j-2*r+1:2*r*j
%                 fd(k) = (b(k)+b(k+2*r))*(q(k)-q0(k))/dt - ...
%                 b(k)*(q(k-2*r)-q0(k-2*r))/dt - ...
%                 b(k+2*r)*(q(k+2*r)-q0(k+2*r))/dt;
%             end
%         end
% 
%         f = f + fd;
    
        %viscous force
        f = f + C*(q-q0)/dt;
        J = J + C/dt;

        %external force for first second
        if i <= 3/dt
            f = f + P;
        end

        %update
        q = q-J\f;
        err = sum(abs(f));
    end

    %update
    u = (q - q0)/dt;
    q0 = q;

    figure(1);
    plot(q(1:2:end),q(2:2:end),'ro-');
    axis equal
    drawnow
    frame = getframe(figure(1));
    im{i} = frame2im(frame);
   
    
end

% filename = 'BrainMSD.gif';
% for idx = 1:299
%     [A,map] = rgb2ind(im{idx+1},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
%     end
% end
