var Y R Pi D Z Eps_R;
varexo Shock_D Shock_Z Shock_R ;
 
parameters Beta Gamma Kappa  Phi_D Phi_Z  Rho_D Rho_Z Rho_R Alpha_Y;
 
 
Beta = 0.99;
Gamma = 1;
Kappa = 0.006;
Phi_D = 0.971;
Phi_Z = -0.463;
Rho_D = 0.992;
Rho_Z = 0.37;
Rho_R = 0.92;
Alpha_Y = 0.9999;

model(linear);


// Euler

Y = Alpha_Y*Y(+1) - (1/Gamma)*(R-Pi(+1)) + D;

// PC

Pi = Beta*Pi(+1) + Kappa*Y + Z;

// Real Rate Rule

R - Pi(+1) = Phi_D*D + Phi_Z*Z + Eps_R;


// Exogenous Processes

D = Rho_D*D(-1) + Shock_D;

Z = Rho_Z*Z(-1) + Shock_Z;

Eps_R = Rho_R*Eps_R(-1) + Shock_R;

end;
 

shocks;
var Shock_D; stderr 0.102;
var Shock_Z; stderr 0.324;
var Shock_R; stderr 0;
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid(1);
steady;
check;

stoch_simul(order = 1,irf=15,nograph, periods = 10000) Y Pi R;

Ymat = R;
Xmat = [Pi Y];
Beta_hat = regress(Ymat,Xmat);
Phi_Pi = Beta_hat(1)
Phi_Y = Beta_hat(2)

BBt(3,1)=(1-(1/Gamma)*Phi_D)/(1-Alpha_Y*Rho_D);
BBt(3,2)=(-(1/Gamma)*Phi_Z)/(1-Alpha_Y*Rho_Z);
BBt(3,3)=-(1/Gamma)*1/(1-Alpha_Y*Rho_R);
BBt(1,1)=1/(1-Beta*Rho_D)*(Kappa*BBt(3,1));
BBt(1,2)=1/(1-Beta*Rho_Z)*(1+Kappa*BBt(3,2));
BBt(1,3)=1/(1-Beta*Rho_R)*(Kappa*BBt(3,3));
BBt(2,1)=Phi_D+Rho_D*BBt(1,1);
BBt(2,2)=Phi_Z+Rho_Z*BBt(1,2);
BBt(2,3)=1+Rho_R*BBt(1,3);

A=zeros(3,3);
B=BBt;

a1  =   B(1,1);
a2  =   B(1,2);
a3  =   B(1,3);
b1  =   B(2,1);
b2  =   B(2,2);
b3  =   B(2,3);
c1  =   B(3,1);
c2  =   B(3,2);
c3  =   B(3,3);

% Implied Taylor Rule

phipi       =   (b2-c2*b1/c1)/(a2-a1*c2/c1)
phiy        =   (b1-a1*phipi)/c1

Determinacy = phiy*((1-Beta)/Kappa) + phipi







