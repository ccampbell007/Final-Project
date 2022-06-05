function   [q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,N,P,dl,EA,EI,impact] = BrainNewtonMethod(Rb,brain_base,N,i,q,q0,q_fixed,tol,free_index,fixed_index,M,dt,u,P,dl,EA,EI,q_stem,L0,Pb,r,n,b,impact,dtheta);
% Runs Newton's Method for the head
global C W x_wall head_rad Rb Pb EAb EIb Cb scale dl Ro_h backimpact

%% Netwon's Raphson Method
% Connect to neck by matching indices
    for ii = 1:length(u)
        u(ii) = 0;
    end
    
    match_index = [1 2];    % Indices connected to the neck
    ii = 0;                 % Counting Variable
    err = 10 * tol;
    fprintf('Time = %f\n',(i-1)*dt);  % Print Time Step
    q_fixed = 0;
    q = q0;                           % Initial Guess
    q(fixed_index) = q_fixed;
    q_free = q(free_index);
    
    %% Newton Raphson Method to Calculate Position Q
    % Free and Fixed Index (Distortion only occurs during impact)
    if impact == 1 || impact == 2 || impact == 3
    while err > tol

    % Interia Term
    f = M/dt * ((q - q0)/dt - u);
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
  f = f + Cb*(q - q0)/dt;
  J = J + Cb/dt;                            

  
 % Weight Term
 % f = f - W;
 

 % Acceleration Term
   f = f - Pb;
    
  %[J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb);

        
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
    end
    
 impact = impact + 1;
    
  % Calculate Impact
  for k = 1:2:length(q)-1 
   if err < tol && impact >= 1
       ii = ii +1;

     % Wall Impact **update with impact**
           %Pb = -M*((q-q0)/dt-u)/dt; %change in momentum mv over dt is force

           Pb(1:2:end) = -25;
           Pb(2:2:end) = 25; 
      [J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb,n,r,b);
      
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
      % Brain Compression/Extension during impact
   end
  end
 else
    end
  Pb(1:end) = 0; 
% Starting point is (0,N). Whatever the neck end has offset, offset all
% head points by

y_offset = q(2) - q_stem(end);           % Y offset Subtract 
x_offset = q(1) - q_stem(end-1);         % X offset for head

if impact == 1 || impact == 2 || impact == 3
% Update position and velocity
  for i = 1:N
   q(2*i-1) = q(2*i-1) - x_offset;                                  % Offset for neck movement
   q(2*i) = q(2*i) - y_offset;                                      % Offset for neck movement
  end    
else
  % Update position and velocity
  for i = 1:N
      if i == 1
   q(2*i-1) = q(2*i-1) - x_offset;                                  % Offset for neck movement
   q(2*i) = q(2*i) - y_offset;                                      % Offset for neck movement
%  q0(2*i-1) = q0(2*i-1) - x_offset;                                % Offset for neck movement
%  q0(2*i) = q0(2*i) - y_offset;                                    % Offset for neck movement
      elseif i > 1 && i <= 4
   q(2*i-1) = q(2*i-1) - x_offset - 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset + 1/scale*cos(dtheta)/2;              % Offset for neck movement
      elseif i > 4 && i <= 8
   q(2*i-1) = q(2*i-1) - x_offset - 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      elseif i > 8 && i <= 12
   q(2*i-1) = q(2*i-1) - x_offset - 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      elseif i > 12 && i <= 16
   q(2*i-1) = q(2*i-1) - x_offset - 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      end
  end
end

    %% Newton Raphson Method to Calculate Position Q
    % Free and Fixed Index (Distortion only occurs during impact)
    if backimpact == 1 || backimpact == 2 || backimpact == 3
    while err > tol

    % Interia Term
    f = M/dt * ((q - q0)/dt - u);
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
  f = f + Cb*(q - q0)/dt;
  J = J + Cb/dt;                            

  
 % Weight Term
 % f = f - W;
 

 % Acceleration Term
   f = f - Pb;
    
  %[J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb);

        
  % Update position
  f_free = f(free_index);
  J_free = J(free_index, free_index);
  q_free = q_free - J_free \ f_free;                    
    
  err = sum(abs(f_free));         % Evaluate error to continue N-R Method
    
  q(free_index) = q_free;
    end
    
 backimpact = backimpact + 1;
    
  % Calculate Impact
  for k = 1:2:length(q)-1 
   if err < tol && backimpact >= 1
       ii = ii +1;

     % Wall Impact **update with impact**
           Pb(1:2:end) = 25;
           Pb(2:2:end) = -25; 
     % Pb = -M*((q-q0)/dt-u)/dt; %change in momentum mv over dt is force

      [J,f] = calculateForceb(q,q0,u,N,M,dt,i,Pb,n,r,b);
      
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
      % Brain Compression/Extension during impact
   end
  end
 else
    end
  Pb(1:end) = 0; 
% Starting point is (0,N). Whatever the neck end has offset, offset all
% head points by

y_offset = q(2) - q_stem(end);           % Y offset Subtract 
x_offset = q(1) - q_stem(end-1);         % X offset for head

   
if impact == 1 || impact == 2 || impact == 3
    
  % Update position and velocity
  for i = 1:N
      if i == 1
   q(2*i-1) = q(2*i-1) - x_offset;                                  % Offset for neck movement
   q(2*i) = q(2*i) - y_offset;                                      % Offset for neck movement
%  q0(2*i-1) = q0(2*i-1) - x_offset;                                % Offset for neck movement
%  q0(2*i) = q0(2*i) - y_offset;                                    % Offset for neck movement
      elseif i > 1 && i <= 4
   q(2*i-1) = q(2*i-1) - x_offset% + 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset% + 1/scale*cos(dtheta)/2;              % Offset for neck movement
      elseif i > 4 && i <= 8
   q(2*i-1) = q(2*i-1) - x_offset% + 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset% + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      elseif i > 8 && i <= 12
   q(2*i-1) = q(2*i-1) - x_offset% + 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset% + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      elseif i > 12 && i <= 16
   q(2*i-1) = q(2*i-1) - x_offset% + 1/scale*sin(dtheta);            % Offset for neck movement
   q(2*i) = q(2*i) - y_offset% + 1/scale*cos(dtheta)/2;               % Offset for neck movement
      end
  end
        
  u = (q - q0)/dt;                % Velocity
  q0 = q;                         % Position
  
end