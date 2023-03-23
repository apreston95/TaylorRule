var Y R Pi D Z Eps_R;
varexo Shock_D Shock_Z Shock_R ;
  

parameters 
Beta ${\beta}$
Gamma ${\gamma}$
Kappa ${\kappa}$ 
Zeta ${\zeta}$ 
Iota ${\iota}$ 
Phi_D ${\phi_{D}}$
Phi_Z ${\phi_{Z}}$
Phi_Y ${\phi_{Y}}$
Phi_Pi ${\phi_{Pi}}$
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
Phi_Pi = 0;
Phi_Y = 0;
Rho_D = 0.992;
Rho_Z = 0.371;
Rho_R = 0.924;
Alpha_Y = 0.9999;
Zeta = 0;
Iota = 0;

model(linear);
#C1 = Zeta/(1+Zeta);

// Euler

// Y = Alpha_Y*((1-Zeta)*Y(+1) + Zeta*Y(-1)) - (1/Gamma)*(R-Pi(+1)) + D;

Y = Alpha_Y*((1-C1)*Y(+1) + C1*Y(-1)) - (1/Gamma)*(R-Pi(+1)) + D;


// PC

Pi = (Beta/(1+Beta*Iota))*Pi(+1) + (Iota/(1+Beta*Iota))*Pi(-1) + Kappa*Y   + Z;

// Real Rate Rule

R - Pi(+1) = Phi_Y*Y(-1) + Phi_Pi*Pi(-1) + Phi_D*D + Phi_Z*Z + Eps_R;


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
Phi_D,uniform_pdf,,,-2,2;
Phi_Z,uniform_pdf,,,-2,2;
Phi_Y,uniform_pdf,,,-2,2;
Phi_Pi,uniform_pdf,,,-2,2;
Iota,beta_pdf,0.5,0.2;
Zeta,beta_pdf,0.5,0.2;
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
 
// estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=8,optim = ('MaxFunEvals',20000),mode_check,mh_replic=500000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;

estimation(mh_tune_jscale,datafile=NK_Data,mode_compute=6,optim=('ncov-mh',500000),mode_check,mh_replic=500000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;


write_latex_prior_table;  

stoch_simul(order = 1,irf=15,nograph,periods = 10000);


