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
Kappa = 0.1;
Phi_D = 0;
Phi_Z = 0;
Phi_R = 0;
Rho_D = 0.9;
Rho_Z = 0.9;
Rho_R = 0.9;
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
var Shock_D; stderr 1;
var Shock_Z; stderr 1;
var Shock_R; stderr 1;
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

identification;
 
estimation(datafile=Pre_Volcker_Data,mode_compute=6,optim=('ncov-mh',500000),cova_compute=1,mode_check,mh_replic=3000000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;

// estimation(mh_tune_jscale,datafile=Pre_Volcker_Data,mode_compute=9,cova_compute=1,mode_check,mh_replic=500000,mh_nblocks=1,mh_drop=0.25,tex) Y R Pi;

write_latex_prior_table;  


posterior_function(function='Determinacy_Condition',sampling_draws = 5000);
Determinacy_Pre=cell2mat(oo_.posterior_function_results(:,1));
[f_Pre,xi_Pre] = ksdensity(Determinacy_Pre);
figure('Name','Determinacy Condition Posterior Distribution')
plot(xi_Pre,f_Pre);

mean_Pre = mean(Determinacy_Pre);
save('Determinacy_Cond_Pre.mat','f_Pre','xi_Pre','Determinacy_Pre','mean_Pre')


