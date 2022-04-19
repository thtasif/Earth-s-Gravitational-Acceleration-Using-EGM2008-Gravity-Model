function [P,scaleFactor] = Associted_Legendre(N,R)

x = R(1,3)/norm(R); %x = sin(phi) of Associated Legendre function (not position x)

%%%%%%%%%%%%%%%%%%%% Constant Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi = zeros(N+1,N+1); %(1,1) actually means (n,m) = (0,0)
kappa = zeros(N+1,N+1);
zi = zeros(N+1,1);

zi(2,1) = sqrt(3); %n = 0
for n = 2 : N     
    zi(n+1,1) = sqrt((2*n+1)/(2*n));  %n = 2 to N    
end


for n = 0:N
     for m = 0:n-1 %as m < n    what if m = n?       
        psi(n+1,m+1) = sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)));      
    end
end


for n = 2 : N
    for m = 0:n-2        
        kappa(n+1,m+1) = psi(n+1,m+1)/psi(n-1+1,m+1);     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(N+1,N+1); %(1,1) actually means (n,m) = (0,0)
P(1,1) = 1; %% n = 0, m = 0, P_00

for n = 1:N %P_n,n
   P(n+1,n+1) = zi(n+1,1)*sqrt(1-x^2)*P(n-1+1,n-1+1); 
end


for n = 1:N %P_n,n-1
   P(n+1,n-1+1) = psi(n+1,n-1+1)*x*P(n-1+1,n-1+1);
end


for n = 2:N %P_n,m    
    for m = 0:n-2       
        P(n+1,m+1) = psi(n+1,m+1)*x*P(n-1+1,m+1) - kappa(n+1,m+1)*P(n-2+1,m+1);               
    end
end

%%%%%%%%%%%%%%%%%%%%%%% scaleFactor Calculation %%%%%%%%%%%%%%%%%%%%%%%%

scaleFactor = zeros(N+1,N+1);
scaleFactor(1,1) = 0;
scaleFactor(2,1) = 1;
scaleFactor(2,2) = 0;

for n = 2:N
    k = n + 1;
    for m = 0:n
        p = m + 1;               
        if (n == m)            
            scaleFactor(k,k,:) = 0;            
        elseif (m == 0)           
            scaleFactor(k,p,:) = sqrt((n+1)*(n)/2);            
        else            
            scaleFactor(k,p,:) = sqrt((n+m+1)*(n-m));            
        end
    end
end


end







