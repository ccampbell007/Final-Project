function [J,f] = calculateForceb(q,q0,u,N,M,dt,i,P,n,r,b)
    global EAb dl EIb Cb

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
   f = f - P;
end