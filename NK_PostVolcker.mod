var Y R Pi D Z Eps_R;
varexo Shock_D Shock_Z Shock_R ;
 
parameters 
Beta ${\beta}$
Gamma ${\gamma}$
Kappa ${\kappa}$ 
Phi_Pi ${\phi_{\Pi}}$
Phi_Y ${\phi_{Y}}$
Rho_D ${\rho_{D}}$
Rho_Z ${\rho_{Z}}$
Rho_R ${\rho_{M}}$
;
 
 
Beta = 0.99;
Gamma = 1;
Kappa = 0.678;
Phi_Pi = 1.769;
Phi_Y = -0.014;
Rho_D = 0.962;
Rho_Z = 0.994;
Rho_R = 0.457;

model(linear);


// Euler

Y = Y(+1) - (1/Gamma)*(R-Pi(+1)) + D;

// PC

Pi = Beta*Pi(+1) + Kappa*Y + Z;

// Taylor Rule

R = Phi_Pi*Pi + Phi_Y*Y + Eps_R;


// Exogenous Processes

D = Rho_D*D(-1) + Shock_D;

Z = Rho_Z*Z(-1) + Shock_Z;

Eps_R = Rho_R*Eps_R(-1) + Shock_R;

end;
 

shocks;
var Shock_D; stderr 0.092;
var Shock_Z; stderr 0.575;
var Shock_R; stderr 0.811;
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid(1);
steady;
check;

stoch_simul(order = 1,irf=15);

 

estimated_params;
Kappa,gamma_pdf,0.1,0.05;
// Gamma, normal_pdf, 1.5,0.375;
Phi_Pi,normal_pdf,1.5,0.25;
Phi_Y,normal_pdf,0.125,0.025;
Rho_D,beta_pdf,0.5,0.2;
Rho_Z,beta_pdf,0.5,0.2;
Rho_R,beta_pdf,0.5,0.2;
stderr Shock_D,inv_gamma_pdf,0.1,2; 
stderr Shock_Z,inv_gamma_pdf,0.1,2; 
stderr Shock_R,inv_gamma_pdf,0.1,2; 

end;

estimated_params_init(use_calibration) ;
end;
 
varobs Pi R Y;

identification;
 
estimation(mh_tune_jscale,datafile=Post_Volcker_Data,mode_compute=6,cova_compute=1,mode_check,mh_replic=1000000,tex,mh_nblocks = 1) Y R Pi;

write_latex_prior_table;  


