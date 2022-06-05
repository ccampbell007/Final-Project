function  [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,P,dl,EA,EI] = StemNewtonMethod(Ls,N,i,q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,P,dl,EA,EI,q_neck,L0,dls,impact);
% Runs Newton's Method for the head
global Cb W

%% Netwon's Raphson Method
    for ii = 1:length(u)
        u(ii) = 0;
    end
if impact == 1 || impact == 2 || impact || 3
    P(1:2:end) = -50;
    P(2:2:end) = 50; 
else
end
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
        
 [J,f] = calculateForces(q,q0,u,N,M,dt,i,P);

   
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
    end
   
% Starting point is (0,N). Whatever the neck end has offset, offset all
% head points by

y_offset = q(2) - q_neck(end);           % Y offset Subtract 
x_offset = q(1) - q_neck(end-1);         % X offset for head
    
  % Update position and velocity (Update for just the first node)
  for i = 1:4
  q(2*i-1) = q(2*i-1) - x_offset;        % Offset for neck movement
  q(2*i) = q(2*i) - y_offset;            % Offset for neck movement
  q0(2*i-1) = q0(2*i-1) - x_offset;        % Offset for neck movement
  q0(2*i) = q0(2*i) - y_offset;            % Offset for neck movement
  end
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position

end