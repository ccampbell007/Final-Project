function [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,P,dl,EAh,EIh,L0] = InitializeHead(Rh,head_base,N,step,q_neck)
% Runs Initializes Matrices to Calculate Head Newton's Method
%% Global variables
global dt RunTime Nsteps g midpoint Ro Ri Ro_h Ri_h L accel ang_accel mu rho rho_s rho_l Y EI EA Yh EAh EIh
%% Properties (**Make global**)
L0 = L;                     % Starting height
Lh = 2*pi*Rh;                % Beam Length (m)
dl = Lh/(N-1);               % Beam Length Between Nodes (m)
I = pi/4*(Ro_h^4 - Ri_h^4); % Mass Moment of Intertia (m^4)
A = pi*(Ro_h^2 - Ri_h^2);   % Beam Cross Sectional Area (m^2)
%% Head
% Head Properties
head_base = L0;                            % Head Base Index
%% Sum Radii
% Nodes
R(1:N) = .005;          % Sphere Radii (m)
R(N) = .025;            % Head base Radius (m)

%% Assign mass, weight and viscous damping arrays
w=zeros(N);         % Initialize sphere weight array
m=zeros(N);          % Initialize sphere mass array
c=zeros(N);         % Initialize Viscous Damping array
for i = 1:N
    m(i) = pi*(Ro_h^2 - Ri_h^2)*Lh*rho/(N-1);  % Sphere Mass
    w(i) = -2000;                              % Sphere Weight 
    c(i) = 6*pi*mu*R(i);                       % Damping Coefficient
end

%% Define head indices
[head_nodes] = circularcoordinates(Rh,head_base,N);
%[brain_index] = circularcoordinates(Rh,brain_base,Nb);
%% Nodes Definition
for i = 1:N
head_index(i) = i;                         % Head Indices
end

%% Applied Accelerations/Forces (*Pull from Main Function Later*)
accel_index = [head_index];     % Acceleration Node Index

for i = 1:length(accel_index)              % Adjust indices to match x coord
    if accel_index(i) == 1
    elseif accel_index(i) == 2
        accel_index(i) = accel_index(i)+1;
    else accel_index(i) = accel_index(i)+2;
    end
end

p = zeros(2*N,1);                          % Iniialize force matrix  
p(accel_index) = accel.*m(accel_index);    % Applied forces (N) F=ma

% Free and Fixed Index
for i = 1:2*N-4
free_index(i) = i+2;                      % Free beam indices
end
fixed_index = [];                      % Fixed beam indices

%% Define Starting Node Locations
nodes = zeros(N,2);                % Initialize Nodes (x,y) Matrix
for i = 1:N
    nodes(i,1) = head_nodes(2*i-1);                % X Locations
    nodes(i,2) = head_nodes(2*i);         % Y = 0 Initial Condition
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

% Damping Matrix
C = zeros(2*N,2*N);        % Initialize Damping Matrix

for i = 1:N
C(2*i-1,2*i-1) = c(i);
C(2*i,2*i) = c(i);
end

% Force Matrix
P = zeros(2*N,1);           % Initialize Vector

for i = 1:2*N
P(i,1) = p(i);
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

tol = Yh*I/Lh^2*1e-3;                 % Tolerance for N-R Method

%% Netwon's Raphson Method
    i = step;
    err = 10 * tol;
    fprintf('Time = %f\n',(i-1)*dt);  % Print Time Step
    q_fixed =0;
    q = q0;                           % Initial Guess
    q(fixed_index) = q_fixed;
    q_free = q(free_index);
    
    % Newton Raphson Method to Calculate Position Q
    while err > tol
        
    % Interia Term
    f = M/dt * ((q - q0)/dt - u);
    J = M/dt^2;                     
    
    % Elastic Forces
    % Use GradEs Function (Derivatives of Stretching Elastic Energy)
    % Linear Springs
    for n = 1:N-1
        xk = q(2*n-1);
        yk = q(2*n);
        xkpl = q(2*n + 1);
        ykpl = q(2*n +2);
        dF = gradEs(xk,yk,xkpl,ykpl,dl,EAh);             % Run GradEs Function[4]
        dJ = hessEs(xk,yk,xkpl,ykpl,dl,EAh);             % Run hessEs Function[5]
        f(2*n-1:2*n + 2) = f(2*n-1:2*n + 2) + dF;       % Add dF to the first 4 terms
        J(2*n-1:2*n + 2,2*n-1:2*n + 2) = ...
            J(2*n-1:2*n + 2,2*n-1:2*n + 2) + dJ;        % Add to Jacobian Matrix
    end
   
    % Use GradEb Function (Derivatives of Bending Elastic Energy)
    % Bending Springs
    for n = 2:N-1
    xkm1 = q(2*n -3);
    ykm1 = q(2*n - 2);
    xk = q(2*n - 1);
    yk = q(2*n);
    xkp1 = q(2*n + 1);
    ykp1 = q(2*n + 2);
    curvature0 = 0;
    dF = gradEb(xkm1,ykm1,xk,yk,xkp1,...
        ykp1,curvature0,dl,EIh);                     % Run gradEb Function[6] 
    dJ = hessEb(xkm1,ykm1,xk,yk,xkp1,...
        ykp1,curvature0,dl,EIh);                     % Run HessEb Function[7] 
    f(2*n - 3:2*n + 2) = f(2*n - 3:2*n + 2) + dF;   % Add dF to the first 4 terms
    J(2*n - 3:2*n + 2,2*n - 3:2*n + 2) = ...
        J(2*n - 3:2*n + 2,2*n - 3:2*n + 2) + dJ;    % Add to Jacobian Matrix
    end
    
  % Viscous Force Term
  % f = f + C*(q - q0)/dt;
  % J = J + C/dt;                            
   
  % Viscous Force
  %Explicit f = f + C*v;
  
 % Weight Term
 % f = f - W;
  
  % Acceleration Term
  f = f - P;
   
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
   end
   
  % Update position and velocity
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  
%   figure(1);
%   hold off 
%   plot(q(1:2:end), q(2:2:end), 'ko-');  % Plot x and y coordinates
%   hold on
%   axis equal
%   drawnow update
    
%   % Store Position and Velocity
%   y_mid_store(i) = q(2*midpoint);        % Store Position for Radius 2
%   v_mid_store(i) = u(2*midpoint);        % Store Velocity for Radius 2
%   y_max(i) = min(q);                     % Max vertical displacement
    
% figure(2)                                % Plot velocity vs time
% time = (1:Nsteps)*dt;
% plot(time, v_mid_store,'k-')
% xlabel('Time (s)')
% ylabel('Velocity(m/s)')
% title('Center Node Velocity vs Time')

end