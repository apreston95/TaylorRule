var C dC N W dW R Pi Pi_W D Z Z_W Eps_R;
varexo Shock_D Shock_Z Shock_R Shock_W;
 
parameters 
Beta ${\beta}$
Psi ${\psi}$
Gamma ${\gamma}$
Kappa ${\kappa}$ 
Kappa_W ${\kappa_{W}}$ 
Phi_Pi ${\phi_{\Pi}}$
Rho_D ${\rho_{D}}$
Rho_Z ${\rho_{Z}}$
Rho_R ${\rho_{M}}$
Rho_W ${\rho_{W}}$
;
 
 
Beta = 0.99;
Psi = 1;
Gamma = 1;
Kappa = 0.678;
Kappa_W = 0.7;
Phi_Pi = 1.769;
Rho_D = 0.962;
Rho_Z = 0.994;
Rho_R = 0.457;
Rho_W = 0.962;


model(linear);


// Euler

C = C(+1) - (1/Gamma)*(R-Pi(+1)) + D;

// PC

Pi = Beta*Pi(+1) + Kappa*W + Z;

// Wage PC

Pi_W = Beta*Pi_W - Kappa_W*(W - (C + Psi*N)) + Z_W;

// Wage Accounting

W = W(-1) + Pi_W - Pi;


// Goods Market clearing

C = W + N;

// Taylor Rule

R = Phi_Pi*Pi + Eps_R;

// Wage Growth

dW = W - W(-1);

// Cons Growth

dC = C - C(-1);


// Exogenous Processes

D = Rho_D*D(-1) + Shock_D;

Z = Rho_Z*Z(-1) + Shock_Z;

Z_W = Rho_W*Z_W(-1) + Shock_W;


Eps_R = Rho_R*Eps_R(-1) + Shock_R;

end;
 

shocks;
var Shock_D; stderr 0.092;
var Shock_Z; stderr 0.575;
var Shock_R; stderr 0.811;
var Shock_W; stderr 0.1;
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
Kappa_W,gamma_pdf,0.1,0.05;
// Gamma, normal_pdf, 1.5,0.375;
Phi_Pi,normal_pdf,1.5,0.25;
Rho_D,beta_pdf,0.5,0.2;
Rho_Z,beta_pdf,0.5,0.2;
Rho_R,beta_pdf,0.5,0.2;
Rho_W,beta_pdf,0.5,0.2;
stderr Shock_D,inv_gamma_pdf,0.1,2; 
stderr Shock_Z,inv_gamma_pdf,0.1,2; 
stderr Shock_R,inv_gamma_pdf,0.1,2; 
stderr Shock_W,inv_gamma_pdf,0.1,2; 


end;

estimated_params_init(use_calibration) ;
end;
 
varobs Pi R dC dW;

identification;
 
estimation(mh_tune_jscale,datafile=HANK_Data,mode_compute=9,cova_compute=1,mode_check,mh_replic=100000,tex) C R Pi;

write_latex_prior_table;  

// model_comparison HANK_Sticky_Wage_RR(0.5) HANK_Sticky_Wage(0.5);


