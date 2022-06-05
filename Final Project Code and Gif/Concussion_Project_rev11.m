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
% [8]   J. G. Beckwith, W. Zhao, S. Ji, A. G. Ajamil, R. P. Bolander, J. J. Chu, T. W. McAllister, J. J. Crisco, S. M. Duma, S. Rowson, S. P. Broglio, K. M. Guskiewicz, J. P. Mihalik, S. Anderson, B. Schnebel, P. Gunnar Brolinson, M. W. Collins, and R. M. Greenwald, “Estimated brain tissue response following impacts associated with and without diagnosed concussion,” Annals of Biomedical Engineering, vol. 46, no. 6, pp. 819–830, 2018. 

%writerObj = VideoWriter('rebound.avi');
%writerObj.FrameRate = 15;
%open(writerObj);

clc; clear all; close all;
%% Global variables
global dt RunTime Nsteps g midpoint Ro Ri Ro_h Ri_h L accel ang_accel mu rho rho_s rho_l Y EI EA Ro_s Ri_s Ls Ns Ro_b Ri_b C W Ys EIs EAs Yh EAh EIh dl x_wall head_rad Rb Pb EIb EAb EIl EAl Cb head_x_offset backimpact

%% General Properties
N = 16;                         % Nodes
dt = .01;                       % Timestep (s)
RunTime = .15;                  % Simulation Runtime (s)
Nsteps = round (RunTime/dt);    % Timesteps 
g = 9.8;                        % Gravity (m/s^2)
midpoint = round((N-1)/2 + 1);  % Middle Node
impact = 0;
backimpact = 0;
oldhead_angle = 0;              % Initiate Old Head Angle
head_angle = 0;                 % Initiate Head Angle
%% Neck
% Neck Beam
L = .25;                        % Neck Length (m)
dl = L/(N-1);                   % Neck Length Between Nodes (m)
Ro = .021;                      % Neck Radius (m)
Ri = .013;                      % Neck Inner Radius (m)
Ro_h = .0875;                   % Head Node Radius (m)
Ri_h = .0875-.0065;             % Head Node Radius (m)
Ro_s = .013;                    % Stem Node Radius (m)
Ri_s = .001;                    % Stem Node Radius (m)
Ro_b = .133;                    % Stem Node Radius (m)
Ri_b = .001;                    % Stem Node Radius (m)
Rl = .15;                       % Helmet Outer diameter (m)
I = pi/4*(Ro^4 - Ri^4);         % Neck Mass Moment of Intertia (m^4)
A = pi*(Ro^2 - Ri^2);           % Neck Cross Sectional Area (m^2)
Ib = pi/4*(Ro_b^4 - Ri_b^4);    % Neck Mass Moment of Intertia (m^4)
Ab = pi*(Ro_b^2 - Ri_b^2);      % Brain Cross Sectional Area (m^2)
Il = pi/4*(Ro_b^4 - Ri_b^4);    % Neck Mass Moment of Intertia (m^4)**Need to updates radii**
Al = pi*(Ro_b^2 - Ri_b^2);      % Neck Cross Sectional Area (m^2)

% Neck Material Properties
Y = 2.5e6;                      % Elastic Modulus (Pa)
EI = Y*I;           
EA = Y*A;
% for Yl = 1.78e2:100:1.78e4 
Yl = 1.78e2;                    % Helmet Elastic Modulus 
EIl = Yl*Il;
EAl = Yl*Al;
Yb = 1.39e2;                    % Brain Elastic Modulus 
EIb = Yb*Ib;
EAb = Yb*Ab;
Yh = 10e7;                      % Brain Stem Elastic Modulus
EIh = Yh*I;
EAh = Yh*A;
Ys = 10e7;                      % Brain Stem Elastic Modulus
EIs = Ys*I;
EAs = Ys*A;
rho_s = 1900;                   % Sphere density (kg/m^3)
rho_l = 1000;                   % Liquid Density (kg/m^3)
rho = rho_s-rho_l;              % Density Difference (kg/m^3)
mu = 1e-3;                      % Viscosity (Pa-s)

%% Head
% Head Properties
Rh = .0285;                     % Head Radius
Rb = .066;                      % Brain Radius

% Head Properties
Nh = 16;                        % Head Nodes
Nb = 16;                        % Brain Nodes
Nl = 16;                        % Helmet Nodes
head_base = L;                  % Head Base Starting Point (y)

% Brain Stem Properies
Ls = .01;                       % Brain Stem Length (m)
Ns = 10;                        % Brain Stem Nodes
dls = Ls/(Ns-1);                % Beam Length Between Nodes (m)
brain_base = L+Ls;              % Brain Base

% Head Offsets
head_y_offset = Ro_h/2 + Ls;
head_x_offset = Ro_h/2 + Ls;
%% Sum Number of Nodes
Nnet = N + Nh + Nb;

% Nodes
R(1:Nnet) = .005;               % Sphere Radii (m)
R(N) = .005;                    % Head base Radius (m)

%% Assign mass, weight and viscous damping arrays
w=zeros(N);                     % Initialize sphere weight array
m=zeros(N);                     % Initialize sphere mass array
c=zeros(N);                     % Initialize Viscous Damping array
for i = 1:N
    m(i) = pi*(Ro^2 - Ri^2)*L*rho/(N-1);                % Sphere Mass
    mb(i) = pi*(Ro_b^2 - Ri_b^2)*(2*pi*Rb)*rho/(Nb-1);  % Brain Node Mass
    w(i) = -1000;                                       % Sphere Weight 
    c(i) = 6*pi*mu*R(i);                                % Damping Coefficient
end

for i = 1:Nb
    mb(i) = pi*(Ro_b^2 - Ri_b^2)*(2*pi*Ro_b)*rho/(Nb-1);  % Brain Node Mass
    cb(i) = 6*pi*mu*R(i);                                 % Brain Damping Coefficient
end
%% Nodes Definition
for i = 1:N
neck_index(2*i-1) = 2*i-1;                   % Neck Indices
neck_index(2*i) = 2*i;                          
end

%% Wall Definition (Turn into function)
x_wall = .3;                                 % Wall location (m)
N_wall = 2;                                  % Wall Nodes

% Establish points for plotting
for i = 1:2
    xw(i) = x_wall;                          % Wall x coordinates
end
yw(1) = 0;                                   % Wall y coordinates
yw(2) = .5;                                  % Wall y coordinates

%% Applied Accelerations/Forces
% for accel = 250:100:2500
accel = 250;                                 % Acceleration (m/s^2)
accelb = 0;
accel_index = [neck_index];                  % Acceleration Node Index

% Angular Acceleration
ang_accel = 100;                             % Rad/s^2

% Forces
p = zeros(2*N,1);                            % Initialize force matrix  
p(accel_index) = accel.*m(1);                % Applied forces (N) F=ma
pb(accel_index) = accelb.*mb(1);             % Brain Applied Forces
% Free and Fixed Index
for i = 1:2*N-5
free_index(i) = i+4;                         % Free beam indices
end
fixed_index = [1,2];                         % Fixed beam indices

%% Define Starting Node Locations
nodes = zeros(N,2);                    % Initialize Nodes (x,y) Matrix
for i = 1:N
    nodes(i,1) = 0;                    % X Locations
    nodes(i,2) = (i-1)*dl;             % Y Initial Condition
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
W(2*i,1) = w(i);            % Assign weight matrix
end

% Damping Matrix
C = zeros(2*N,2*N);        % Initialize Damping Matrix

for i = 1:N
C(2*i-1,2*i-1) = c(i);     % Assign Damping Matrix
C(2*i,2*i) = c(i);
end

for i = 1:Nb
Cb(2*i-1,2*i-1) = c(i);    % Assign Brain Damping Matrix
Cb(2*i,2*i) = c(i);
end

% Force Matrix
P = zeros(2*N,1);          % Initialize Vector

for i = 1:2*N
P(i,1) = p(i);             % Assign force vector
end

for i = 1:2*Nb
Pb(i,1) = pb(1);           % Brain Force Vector
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
q_b = zeros(2*N,1);                 % Initialize Brain Vector
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
    head_rad = Ro_h;

    % Calculating Forces
    [J,f] = calculateForce(q,q0,u,N,M,dt,i,P);

    % Update position
    f_free = f(free_index);
    J_free = J(free_index, free_index);
    q_free = q_free - J_free \ f_free;  
  
    q(free_index) = q_free;
   
    err = sum(abs(f_free));  

        %% Deviation at Wall (Helmet)
   % If the err is below tolerance (simulation has finished calculating the
   % next q value) and the x value of the last node is over the imaginary
   % wall, then the x value of that node is set to the wall value, it is no
   % longer a free index, and the q/err is recalculated
   if err < tol && q(end-1) > x_wall-Rl-head_x_offset
      q(end-1) = x_wall-Rl-head_x_offset; 
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
      impact = 1;
   else
       continue;
   end
   
    
   %% Deviation at Wall
   % If the err is below tolerance (simulation has finished calculating the
   % next q value) and the x value of the last node is over the imaginary
   % wall, then the x value of that node is set to the wall value, it is no
   % longer a free index, and the q/err is recalculated
   if err < tol && q(end-1) > x_wall-head_rad-head_x_offser
      q(end-1) = x_wall-head_rad-head_x_offset; 
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
      impact = 1;
   else
       continue;
   end
    end

   % Brain Impact
   for k = 1:2:length(q_b)-1
    if err < tol && (q_b(k) <= (q(end-1) - .08))
    backimpact = 1;
    else
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
   
  ang = 1;                                          % Determines if the angle is pos or negative
  % Head Angle for brain orientation
  if (q(end-1) >= q(end-3)) && ang == 1             % Positive head angle
  head_angle = tan(abs(q(end) - q(end-2))/abs(q(end-1) - q(end-3)));
  dtheta = head_angle - oldhead_angle;
  ang = 1;
  elseif (q(end-1) < q(end-3)) && ang == 0          % Negative head angle
  head_angle = tan(abs(q(end) - q(end-2))/abs(q(end-3) - q(end-1)));
  dtheta = head_angle - oldhead_angle;
  ang = 0;
  elseif (q(end-1) < q(end-3)) && ang == 1         % Head moves from positive to negative
  head_angle = tan(abs(q(end) - q(end-2))/abs(q(end-1) - q(end-3)));
  dtheta = head_angle - oldhead_angle;
  ang = 0;
  elseif (q(end-1) > q(end-3)) && ang == 0         % Head moves from negative to positive
  head_angle = tan(abs(q(end) - q(end-2))/abs(q(end-1) - q(end-3)));
  dtheta = head_angle - oldhead_angle;
  ang = 1;
  end
  
  % Brain Function
  if i == 2
  [q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,N_b,P_b,dl_b,EA_b,EI_b,L0_b,rb,nb,bb,impact_b] = InitializeBrain(Rb,brain_base,Nb,i,q);
  else
  [q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,N_b,P_b,dl_b,EA_b,EI_b,impact] = BrainNewtonMethod(Rb,brain_base,N_b,i,q_b,q0_b,q_fixed_b,tol_b,free_index_b,fixed_index_b,M_b,dt_b,u_b,P,dl_b,EA_b,EI_b,q,L0_b,P_b,rb,nb,bb,impact,dtheta);
  end
  
  % Redefine Head Angle
  oldhead_angle = head_angle;
  
  % Static Head
  num_of_angles = 20;
  angles = 0:2*pi/num_of_angles:2*pi;
  x_head = zeros(length(angles),1);
  y_head = zeros(length(angles),1);  
  for j = 1:length(angles)
    x_head(j) = head_rad*cos(angles(j))+head_x_offset;
    y_head(j) = head_rad*sin(angles(j))+head_y_offset;
  end
  
  % Helmet
  helmet_y_offset = Rl/2;
  head_x_offset = head_rad/2;
  num_of_angles = 20;
  angles = 0:2*pi/num_of_angles:pi;
  x_helmet = zeros(length(angles),1);
  y_helmet = zeros(length(angles),1);  
  for j = 1:length(angles)
    x_helmet(j) = Rl*cos(angles(j))+head_x_offset;
    y_helmet(j) = Rl*sin(angles(j))+helmet_y_offset;
  end
    
  
  figure(1);
  hold off 
  plot(q(1:2:end), q(2:2:end), 'ko-');  % Plot neck
  hold on
  plot(x_head+q(end-1),y_head+q(end),'b'); % Plot Head
  plot(x_helmet+q(end-1),y_helmet+q(end),'k--'); % Plot Helmet
  plot(xw,yw,'ko-')
  plot(q_b(1:2:end), q_b(2:2:end), 'ko-');  % Plot brain
  axis equal
  drawnow update
    
  % Store Position and Velocity
  y_mid_store(i) = q(2*midpoint);        % Store Position for Radius 2
  v_mid_store(i) = u_b(2*midpoint);      % Store Velocity for Radius 2
  y_max(i) = min(q);                     % Max vertical displacement
  if i == 2
  else
  a_mid_store(i-2) = (v_mid_store(i)-v_mid_store(i-1))/dt;  % Brain Acceleration
  end
  
%   frame = getframe(figure(1));
%   im{i} = frame2im(frame);
%    
 end 
% 
% % Create gif
% filename = 'BrainMSD.gif';
% for idx = 1:14
%     [A,map] = rgb2ind(im{idx+1},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
%     end
% end

%     
figure(2)                                % Plot velocity vs time
time = (1:Nsteps)*dt;
plot(time, v_mid_store,'k-')
xlabel('Time (s)')
ylabel('Velocity(m/s)')
title('Brain Velocity vs Time')

figure(3)                                % Plot acceleration vs time
time = (2:Nsteps-1)*dt;
plot(time, a_mid_store,'k-')
xlabel('Time (s)')
ylabel('Acceleration(m/s^2)')
title('Brain Acceleration vs Time')


% set the seconds per image
% open the video writer


%close(writerObj);