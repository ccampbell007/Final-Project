function   [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,Pb,dl,EAb,EIb,L0,r,n,b,impact] = InitializeBrain(Rb,brain_base,N,step,q_stem);
% Runs Initializes Matrices to Calculate Brain Newton's Method
global dt RunTime Nsteps g midpoint Ro Ri Ro_h Ri_h L accel ang_accel mu rho rho_s rho_l Y EIb EAb Ro_s Ri_s Ro_b Ri_b Pb scale
L0 = L;                 % Starting height
Lb = 2*pi*Rb;            % Beam Length (m)
dl = Lb/(N-1);           % Beam Length Between Nodes (m)
I = pi/4*(Ro_b^4 - Ri_b^4); % Mass Moment of Intertia (m^4)
A = pi*(Ro_b^2 - Ri_b^2);   % Beam Cross Sectional Area (m^2)
r = 4;                      % Rows/cols of nodes
n = r*r;                    % Number of nodes
N=n;
impact = 0;
% Head Properties

head_base = L0;                            % Head Base Index
%% Sum Radii
% Nodes
R(1:N) = .005;                             % Sphere Radii (m)
R(N) = .005;                               % Head base Radius (m)

%% Assign mass, weight and viscous damping arrays
w=zeros(N);         % Initialize sphere weight array
m=zeros(N);          % Initialize sphere mass array
c=zeros(N);         % Initialize Viscous Damping array
for i = 1:N
    m(i) = pi*(Ro_b^2 - Ri_b^2)*Lb*rho/(N-1)*2;  % Sphere Mass
    w(i) = -2000;                               % Sphere Weight 
    c(i) = 6*pi*mu*R(i);                        % Damping Coefficient
end

%% Define brain indices

%initial pos
x = zeros(n,2);
scale = 30;         % Variable to scale brain size **(Set to brain size later)**
j = 0; i = r-3;
for k = 1:n
    x(k,1) = j/scale; x(k,2) = L0 + i/scale;
    j = j+1;
    if j == r
        j = 0; i = i + 1;
    end
end

%initial dof vector
q0 = zeros(2*n,1);
for i = 1:n
    q0(2*i-1) = x(i,1);
    q0(2*i) = x(i,2);
end

%% Nodes Definition
for i = 1:N
head_index(i) = i;                         % Head Indices
end

% Free and Fixed Index
for ii = 1:2*N
    if ii <= 2 || (isreal((ii-1)/8) == 1 && rem((ii-1)/8,1) == 0) %|| (isreal((ii)/2) && rem((ii)/2,1)==0)        % Ignore brain stem connected nodes    
    else
    free_index(ii) = ii;                                   % Free beam indices
    end
end
[zero_index] = find(~free_index);
free_index = nonzeros(free_index)';                        % Remove zero terms
fixed_index = [];
%% Define Starting Node Locations
nodes = zeros(N,2);                                       % Initialize Nodes (x,y) Matrix
for i = 1:N
    nodes(i,1) = q0(2*i-1);                               % X Locations
    nodes(i,2) = q0(2*i);                  % Y = 0 Initial Condition
end

%% Equations of Motion

% Mass Matrix
M = zeros(2*N,2*N);         % Initialize Mass Matrix

for i = 1:N
M(2*i-1,2*i-1) = m(i);
M(2*i,2*i) = m(i);
end

% Weight Matrix
W = zeros(2*N,1);           % Initialize Weight Matrix

for i = 1:N
W(2*i,1) = w(i);
end

%damper
b = zeros(2*n,1);
for k = 1:2*n
    b(k) = 10;
end

% Damping Matrix
C = zeros(2*N,2*N);        % Initialize Damping Matrix

for i = 1:N
C(2*i-1,2*i-1) = c(i);
C(2*i,2*i) = c(i);
end

% Initial DOF Vector
q0 = zeros(2*N,1);          % Initialize Vector

for i = 1:N
q0(2*i - 1) = nodes(i,1);   % X Coordinate
q0(2*i) = nodes(i,2);       % Y Coordinate
end

% Initial Position and Velocity
q = q0;                     % DOF
u = (q - q0)/dt;            % Velocity

%  Initialize Store Vectors
y_mid_store = zeros(Nsteps,1);
v_mid_store = zeros(Nsteps,1);

y_mid_store(1) = q(2*midpoint);     % Store Position at First Timestep
v_mid_store(1) = u(2*midpoint);     % Store Velocity at First Timestep

tol = Y*I/Lb^2*1e-3;                 % Tolerance for N-R Method

%% Netwon's Raphson Method
    i = step;
    err = 10 * tol;
    fprintf('Time = %f\n',(i-1)*dt);  % Print Time Step
    q_fixed = 0;
    q = q0;                           % Initial Guess
    q(fixed_index) = q_fixed;
    q_free = q(free_index);
    
    % Newton Raphson Method to Calculate Position Q
    while err > tol
        
    % Interia Term
    f = M/dt * ((q - q0)/dt - u);
    J = M/dt^2;                     
    
   for k = 1:n-1
            if mod(k,r) == 0
                dF = zeros(6,1); dJ = zeros(6);
            else
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
             
            l_k = 1;
            dF = gradEs(xk,yk,xkp1,ykp1,l_k,EAb);
            dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EAb);
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
            dF = gradEs(xk,yk,xkp1,ykp1,l_k,EAb);
            dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EAb);
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
                dF = gradEs(xk,yk,xkp1,ykp1,l_k,EAb);
                dJ = hessEs(xk,yk,xkp1,ykp1,l_k,EAb);
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
    
  % Viscous Force Term
  %f = f + C*(q - q0)/dt;
  %J = J + C/dt;                            
  
 % Weight Term
 % f = f - W;
  
  % Acceleration Term
  %f = f - Pb;
   
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
    end
   
y_offset = q(2) - q_stem(end);           % Y offset Subtract 
x_offset = q(1) - q_stem(end-1);         % X offset for head
    
  % Update position and velocity
  for i = 1:N
   q(2*i-1) = q(2*i-1) - x_offset;                                  % Offset for neck movement
   q(2*i) = q(2*i) - y_offset;                                      % Offset for neck movement
  end
   
  % Update position and velocity
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  


end