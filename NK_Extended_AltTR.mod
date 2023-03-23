var Y R Pi D Z Eps_R;
varexo Shock_D Shock_Z Shock_R ;
  

parameters 
Beta ${\beta}$
Gamma ${\gamma}$
Kappa ${\kappa}$ 
Zeta ${\zeta}$ 
Iota ${\iota}$ 
Phi_Y ${\phi_{Y}}$
Phi_LagPi ${\phi_{Pi_{-1}}}$
Phi_LagY ${\phi_{Y_{-1}}}$
Phi_Pi ${\phi_{Pi}}$
Rho_D ${\rho_{D}}$
Rho_Z ${\rho_{Z}}$
Rho_R ${\rho_{M}}$
;
 
Beta = 0.99;
Gamma = 1;
Kappa = 0.7;
Phi_Pi = 1.5;
Phi_Y = 0.125;
Phi_LagPi = 0;
Phi_LagY = 0;
Rho_D = 0.992;
Rho_Z = 0.371;
Rho_R = 0.924;
Zeta = 0.2;
Iota = 0.2;

model(linear);
#C1 = Zeta/(1+Zeta);

// Euler


Y = (1-C1)*Y(+1) + C1*Y(-1) - (1/Gamma)*(R-Pi(+1)) + D;


// PC

Pi = (Beta/(1+Beta*Iota))*Pi(+1) + (Iota/(1+Beta*Iota))*Pi(-1) + Kappa*Y   + Z;

// Taylor Rule

R  = Phi_Y*Y + Phi_Pi*Pi + Phi_LagY*Y(-1) + Phi_LagPi*Pi(-1) + Eps_R;


// Exogenous Processes

D = Rho_D*D(-1) + Shock_D;

Z = Rho_Z*Z(-1) + Shock_Z;

Eps_R = Rho_R*Eps_R(-1) + Shock_R;

end;
 

shocks;
var Shock_D; stderr 0.102;
var Shock_Z; stderr 0.324;
var Shock_R; stderr 0.057;
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid(1);
steady;
check;

stoch_simul(order = 1,irf=15,nograph);

 

estimated_params;
Kappa,gamma_pdf,0.1,0.05;
// Gamma, normal_pdf, 1.5,0.375;
Phi_Pi,normal_pdf,1.5,0.25;
Phi_Y,normal_pdf,0.125,0.05;
Phi_LagPi,normal_pdf,0,0.25;
Phi_LagY,normal_pdf,0,0.25;
Iota,beta_pdf,0.5,0.2;
Zeta,beta_pdf,0.5,0.2;
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

dynare_sensitivity;

identification;
 
estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=9,mode_check,mh_replic=500000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;

// estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=6,optim=('ncov-mh',500000),mode_check,mh_replic=500000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;


write_latex_prior_table;  

stoch_simul(order = 1,irf=15,nograph,periods = 10000);



