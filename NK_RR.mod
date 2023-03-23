var Y R Pi D Z Eps_R;
varexo Shock_D Shock_Z Shock_R ;
  

parameters 
Beta ${\beta}$
Gamma ${\gamma}$
Kappa ${\kappa}$ 
Phi_D ${\phi_{D}}$
Phi_Z ${\phi_{Z}}$
Rho_D ${\rho_{D}}$
Rho_Z ${\rho_{Z}}$
Rho_R ${\rho_{M}}$
Alpha_Y ${\alpha_{Y}}$
;
 
Beta = 0.99;
Gamma = 1;
Kappa = 0.006;
Phi_D = 0.971;
Phi_Z = -0.463;
Rho_D = 0.992;
Rho_Z = 0.371;
Rho_R = 0.924;
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
Phi_D,uniform_pdf,,,-2,2;
Phi_Z,uniform_pdf,,,-2,2;
Rho_D,beta_pdf,0.5,0.2;
Rho_Z,beta_pdf,0.5,0.2;
Rho_R,beta_pdf,0.5,0.2;
stderr Shock_D,inv_gamma_pdf,0.1,2; 
stderr Shock_Z,inv_gamma_pdf,0.1,2; 
stderr Shock_R,inv_gamma_pdf,0.1,2; 

end;

// estimated_params_init(use_calibration) ;
// end;
 
varobs Pi R Y;

dynare_sensitivity;

identification;
 
// estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=9,mode_check,mh_replic=300000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;
estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=6,mode_check,mh_replic=300000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;


write_latex_prior_table;  



model_comparison NK_RR(0.5) NK_Uniform(0.5);



