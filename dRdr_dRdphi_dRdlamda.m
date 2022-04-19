function [dRdr,dRdphi,dRdlamda] = dRdr_dRdphi_dRdlamda(N_ii,M_ii,r_eq,R,meu,C_nm,S_nm)

norm_R = sqrt(R(1,:)*R(1,:)');
rho = sqrt(R(1,1)^2 + R(1,2)^2);
lamda = atan2(R(1,2),R(1,1));
phi = asin(R(1,3)/norm_R);

F(1) = (R(1,:)*R(1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,scaleFactor] = Associted_Legendre(N_ii+2,R(1,:)); 

drdR = [R(1,1)/norm_R, R(1,2)/norm_R, R(1,3)/norm_R];

sum_dRdr = 0;

for n = 2:N_ii       
    G_r = F(1)^(-(n+2)/2);
    
    sum_dRdr1 = 0;

     for m = 0:M_ii              
        cos_sin_f = C_nm(n+1,m+1)*cos(m*lamda) + S_nm(n+1,m+1)*sin(m*lamda);
                
        sum_dRdr1 = sum_dRdr1 + P(n+1,m+1)*cos_sin_f;      
    end
    
    sum_dRdr = sum_dRdr + (n+1)*r_eq^n*G_r*sum_dRdr1;
       
end

dRdr = -meu*sum_dRdr*drdR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dphidR = (1/rho)*[-(R(1,1)/norm_R)*(R(1,3)/norm_R), -(R(1,2)/norm_R)*(R(1,3)/norm_R), rho^2/norm_R^2];

sum_dRdphi = 0;

for n = 2:N_ii
   for m = 0:M_ii
       
       g_n1 = F(1)^(-(n+1)/2);
       tan_phi = R(1,3)/sqrt(R(1,1)^2 + R(1,2)^2);
       
       sum_dRdphi = sum_dRdphi + r_eq^n*g_n1*(P(n+1,m+1+1)*scaleFactor(n+1,m+1) - m*tan_phi*P(n+1,m+1))*(C_nm(n+1,m+1)*cos(m*lamda) + S_nm(n+1,m+1)*sin(m*lamda));
    
   end
end

dRdphi = meu*sum_dRdphi*dphidR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dlamdadR = (1/rho)*[-R(1,2)/rho, R(1,1)/rho, 0];

sum_dRdlamda = 0;

for n = 2:N_ii
        
     G_lamda = F(1)^(-(n+1)/2);
    
    sum_dRdlamda1 = 0;

     for m = 0:M_ii
                
        del_cos_sin_f = m*(-C_nm(n+1,m+1)*sin(m*lamda) + S_nm(n+1,m+1)*cos(m*lamda));
        
        sum_dRdlamda1 = sum_dRdlamda1 + P(n+1,m+1)*del_cos_sin_f;
        
    end
    
       sum_dRdlamda = sum_dRdlamda + r_eq^n*G_lamda*sum_dRdlamda1;

end

dRdlamda = meu*sum_dRdlamda*dlamdadR;

end


