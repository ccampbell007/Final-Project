function [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,P,dl,EA,EI] = HeadNewtonMethod(Rh,head_base,N,i,q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,P,dl,EA,EI,q_neck,L0)
%% Netwon's Raphson Method
% Connect to neck by matching indices
match_index = [1 2];    % Indices connected to the neck

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
        dF = gradEs(xk,yk,xkpl,ykpl,dl,EA);             % Run GradEs Function[4]
        dJ = hessEs(xk,yk,xkpl,ykpl,dl,EA);             % Run hessEs Function[5]
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
        ykp1,curvature0,dl,EI);                     % Run gradEb Function[6] 
    dJ = hessEb(xkm1,ykm1,xk,yk,xkp1,...
        ykp1,curvature0,dl,EI);                     % Run HessEb Function[7] 
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
   
% Starting point is (0,N). Whatever the neck end has offset, offset all
% head points by

y_offset = q(2) - q_neck(end);  % Y offset Subtract 
x_offset = q(1) - q_neck(end-1);         % X offset for head
    
  % Update position and velocity
  for i = 1:N
  q(2*i-1) = q(2*i-1) - x_offset;        % Offset for neck movement
  q(2*i) = q(2*i) - y_offset;          % Offset for neck movement
  end
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  
%   figure(1);
%   hold off 
%   plot(q(1:2:end), q(2:2:end), 'ko-');  % Plot x and y coordinates
%   hold on
%   axis equal
%   drawnow update
    
  % Store Position and Velocity
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