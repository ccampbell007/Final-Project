function   [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,P,dl,EA,EI] = HelmetNewtonMethod(Rb,brain_base,N,i,q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,P,dl,EA,EI,q_stem,L0,Pb);
% Runs Newton's Method for the head
global C W x_wall head_rad Rb Pb

%% Netwon's Raphson Method
% Connect to neck by matching indices
match_index = [1 2];    % Indices connected to the neck
    ii = 0;             % Counting Variable
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
    
    [J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb);

        
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
  
  % Calculate Impact
  for k = 1:2:length(q)-1 
   if err < tol && q(k) > x_wall-.03 && ii < 10
       ii = ii +1;
     % Wall Impact **update with impact**
           Pb(1:2:end) = -2000;
           Pb(2:2:end) = 2000; 
      [J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb);
      
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
    end
  Pb(1:end) = 0; 
% Starting point is (0,N). Whatever the neck end has offset, offset all
% head points by

y_offset = q(2) - q_stem(end);  % Y offset Subtract 
x_offset = q(1) - q_stem(end-1);         % X offset for head
    
  % Update position and velocity
  for i = 1:N
  q(2*i-1) = q(2*i-1) - x_offset;        % Offset for neck movement
  q(2*i) = q(2*i) - y_offset;          % Offset for neck movement
  end
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  
end