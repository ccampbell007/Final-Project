function [J,f] = calculateForces(q,q0,u,N,M,dt,i,P)
    global EA dl EI

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
   
  
 % Weight Term
 % f = f - W;
 

  % Acceleration Term
  if i == 1 | i == 2 | i == 3
    f = f - P;
  end

end