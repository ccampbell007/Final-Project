% MAE 259B Project
% Concussion Simulation

% REFERENCES
% 
% [1] 	K. M. Jawed, "Conservative Force and Potential Energy," in Discrete Simulation of Slender Structures. 
% [2] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 3, Los Angeles: University of California Los Angeles, 2022. 
% [3] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 4, Los Angeles: University of California Los Angeles, 2022.
% [4] 	K. M. Jawed, gradEs, University of California Los Angeles, 2022. 
% [5] 	K. M. Jawed, gradEb, University of California Los Angeles, 2022. 
% [6] 	K. M. Jawed, hessEs, University of California Los Angeles, 2022. 
% [7] 	K. M. Jawed, hessEb, University of California Los Angeles, 2022. 
clc; clear all; close all;

%writerObj = VideoWriter('rebound.avi');
%writerObj.FrameRate = 15;
%open(writerObj);

%% General Properties
N = 5;                          % Nodes
dt = .005;                       % Timestep (s)
RunTime = 0.15;                    % Simulation Runtime (s)
Nsteps = round (RunTime/dt);    % Timesteps 
g = 9.8;                        % Gravity (m/s^2)
midpoint = round((N-1)/2 + 1);  % Middle Node

%% Neck
% Neck Beam
L = 0.3;                  % Beam Length (m)
dl = L/(N-1);           % Beam Length Between Nodes (m)
Ro = .013;              % Beam Radius (m)
Ri = .011;              % Beam Inner Radius (m)
I = pi/4*(Ro^4 - Ri^4); % Mass Moment of Intertia (m^4)
A = pi*(Ro^2 - Ri^2);   % Beam Cross Sectional Area (m^2)

% Material Properties
Y = 10e6;           % Elastic Modulus (Pa)
EI = Y*I;           
EA = Y*A;
rho_s = 3000;       % Sphere density (kg/m^3)
rho_l = 1000;       % Liquid Density (kg/m^3)
rho = rho_s-rho_l;  % Density Difference (kg/m^3)
mu = 1000;          % Viscosity (Pa-s)

%% Head
% Head Properties
Rh = .075;                                % Head Radius
Rb = .03;                                 % Brain Radius

% Head Properties
Nh = 5;                                   % Head Nodes
Nb = 5;                                   % Brain Nodes
head_base = L;                            % Head Base Starting Point (y)

% Brain Stem Properies
Ls = .1;                                   % Brain Stem Length (m)
Ns = 5;                                    % Brain Stem Nodes
dls = Ls/(Ns-1);                           % Beam Length Between Nodes (m)
brain_base = L+Ls;                         % Brain Base
%% Sum Number of Nodes
Nnet = N + Nh + Nb;
% Nodes
R(1:Nnet) = .005;       % Sphere Radii (m)
R(N) = .025;            % Head base Radius (m)

%% Assign mass, weight and viscous damping arrays
w=zeros(N);         % Initialize sphere weight array
m=zeros(N);          % Initialize sphere mass array
c=zeros(N);         % Initialize Viscous Damping array
for i = 1:N
    m(i) = pi*(Ro^2 - Ri^2)*L*rho/(N-1);  % Sphere Mass
    w(i) = -2000;                         % Sphere Weight 
    c(i) = 6*pi*mu*R(i);                  % Damping Coefficient
end

%% Nodes Definition
for i = 1:N-1
    neck_index(i) = i;                         % Neck Indices
end
head_base = N;                             % Head Index
brain_base = N;

%% Wall Definition (Turn into function)
x_wall = 0.3;                             % Wall location (m)
N_wall = 2;                               % Wall Nodes


%% Applied Accelerations/Forces
accel = 2200;                              % Acceleration (m/s^2)
accel_index = [neck_index];     % Acceleration Node Index

for i = 1:length(accel_index)              % Adjust indices to match x coord
    if accel_index(i) == 1
    elseif accel_index(i) == 2
        accel_index(i) = accel_index(i)+1;
    else accel_index(i) = accel_index(i)+2;
    end
end
ang_accel = 100;                           % Rad/s^2

p = zeros(2*N,1);                          % Iniialize force matrix  
p(accel_index) = accel.*m(accel_index);    % Applied forces (N) F=ma

% Free and Fixed Index
for i = 1:2*N-2
    free_index(i) = i+2;                      % Free beam indices
end
fixed_index = [1,2];                      % Fixed beam indices

%% Define Starting Node Locations
nodes = zeros(N,2);                % Initialize Nodes (x,y) Matrix
for i = 1:N
    nodes(i,1) = 0;                % X Locations
    nodes(i,2) = (i-1)*dl;         % Y = 0 Initial Condition
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

tol = Y*I/L^2*1e-3;                 % Tolerance for N-R Method

global dl EA EI



%% Netwon's Raphson Method
for i=2:Nsteps
    err = 10 * tol;
    fprintf('Time = %f\n',(i-1)*dt);  % Print Time Step
    q_fixed =0;
    q = q0;                           % Initial Guess
    q(fixed_index) = q_fixed;
    free_index = [3:2*N];
    q_free = q(free_index);
    
    % Newton Raphson Method to Calculate Position Q
    while err > tol

    %% Normal Method for calculing next q
    fprintf('err = %f\n',err);

    % Defining Radius of the head
    head_rad = 0.17;

    % Calculating Forces
    [J,f] = calculateForce(q,q0,u,N,M,dt,i,P);

    % Update position
    f_free = f(free_index);
    J_free = J(free_index, free_index);
    q_free = q_free - J_free \ f_free;  
  
    q(free_index) = q_free;
   
    err = sum(abs(f_free));  

   %% Deviation at Wall
   % If the err is below tolerance (simulation has finished calculating the
   % next q value) and the x/y values of the last node are past the imaginary
   % circle, then the x value of that node is set to the wall value, it is no
   % longer a free index, and the q/err is recalculated
   if err < tol && ((q(end)-head_rad)^2 + q(end-1)^2) > head_rad^2
      bad_rad = sqrt((q(end)-head_rad)^2 + q(end-1)^2);
      scale = head_rad/bad_rad;
      q(end-1) = q(end-1)*scale;
      q(end) = (q(end)-head_rad)*scale + head_rad;
      [J,f] = calculateForce(q,q0,u,N,M,dt,i,P);

      % Redefine the free index
      free_index = [3:2*N-2 2*N];

      % Update values
      f_free = f(free_index);
      J_free = J(free_index, free_index);
      q_free = q(free_index);
      q_free = q_free - J_free \ f_free;  
      q(free_index) = q_free;
      err = sum(abs(f_free));
      % The new q vector calculated like this anytime the last node crosses
      % the imaginary wall. THE FREE INDEX IS RESET AFTER THIS. That way,
      % once the last node starts moving to the left again (at err < 0),
      % the simulation continues as normal
   else
       continue;
   end
 
   
   end



   
  % Update position and velocity
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  
  % Head Function
  if i == 2
  [q_h,q0_h,q_fixed_h,tol_h,free_index_h,fixed_index_h,M_h,dt_h,u_h,N_h,P_h,dl_h,EA_h,EI_h,L0_h] = InitializeHead(Rh,head_base,N,i,q);
  else
  [q_h,q0_h,q_fixed_h,tol_h,free_index_h,fixed_index_h,M_h,dt_h,u_h,N_h,P_h,dl_h,EA_h,EI_h] = HeadNewtonMethod(Rh,head_base,N_h,i,q_h,q0_h,q_fixed_h,tol_h,free_index_h,fixed_index_h,M_h,dt_h,u_h,P_h,dl_h,EA_h,EI_h,q,L0_h);
  end
  
  % Brain Stem Function
  if i == 2
  [q_s,q0_s,q_fixed_s,tol_s,free_index_s,fixed_index_s,M_s,dt_s,u_s,N_s,P_s,dl_s,EA_s,EI_s,L0_s] = InitializeStem(Ls,N,i,q,dls);
  else
  [q_s,q0_s,q_fixed_s,tol_s,free_index_s,fixed_index_s,M_s,dt_s,u_s,N_s,P_s,dl_s,EA_s,EI_s] = StemNewtonMethod(Ls,N_s,i,q_s,q0_s,q_fixed_s,tol_s,free_index_s,fixed_index_s,M_s,dt_s,u_s,P_s,dl_s,EA_s,EI_s,q,L0_s,dls);
  end
  
  % Brain Function
  if i == 2
  [q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,N_b,P_b,dl_b,EA_b,EI_b,L0_b] = InitializeBrain(Rb,brain_base,Nb,i,q);
  else
  [q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,N_b,P_b,dl_b,EA_b,EI_b] = BrainNewtonMethod(Rb,brain_base,N_b,i,q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,P_b,dl_b,EA_b,EI_b,q_s,L0_b);
  end

  % Plotting the skull
  head_rad = 0.17;
  num_of_angles = 20;
  angles = 0:2*pi/num_of_angles:2*pi;
  x_skull = zeros(length(angles),1);
  y_skull = zeros(length(angles),1);  
  for j = 1:length(angles)
    x_skull(j) = head_rad*cos(angles(j));
    y_skull(j) = head_rad*sin(angles(j));
  end
  
  fig = figure(1);
  hold off 
  plot(q(1:2:end), q(2:2:end), 'ko-');  % Plot neck
  hold on
  plot(x_skull+q(1),y_skull+q(2)+head_rad,'b');
  axis equal
%  Frame = getframe(fig);
%  writeVideo(writerObj, Frame);
  drawnow update
  

  % Store Position and Velocity
  y_mid_store(i) = q(2*midpoint);        % Store Position for Radius 2
  v_mid_store(i) = u(2*midpoint);        % Store Velocity for Radius 2
  y_max(i) = min(q);                     % Max vertical displacement
end 


% set the seconds per image
% open the video writer


%close(writerObj);
