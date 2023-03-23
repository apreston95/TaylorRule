/*
 * modified by Franck Portier and Andy Preston (May 2021) to add Euler discouting and a real rate rule
 * 
 * This file provides replication files for 
 * Smets, Frank and Wouters, Rafael (2007): "Shocks and Frictions in US Business Cycles: A Bayesian
 * DSGE Approach", American Economic Review, 97(3), 586-606, that are compatible with Dynare 4.5 onwards
 *
 * To replicate the full results, you have to get back to the original replication files available at
 * https://www.aeaweb.org/articles.php?doi=10.1257/aer.97.3.586 and include the respective estimation commands and mode-files.
 *
 * Notes:
 *  - The consumption Euler equation in the paper, equation (2), premultiplies the risk premium process \varepsilon_t^b,
 *      denoted by b in this code, by the coefficient c_3. In the code this prefactor is omitted by setting the 
 *      coefficient to 1. As a consequence, b in this code actually is b:=c_3*\varepsilon_t^b. As a consequence, in 
 *      the arbitrage equation for the value of capital in the paper, equation (4), the term 1*\varepsilon_t^b
 *      is replaced by 1/c_3*b, which is equal to \varepsilon_t^b given the above redefinition. This rescaling also explains why the 
 *      standard deviation of the risk premium shock in the AR(1)-process for b has a different standard deviation than reported
 *      in the paper. However, the results are unaffected by this scaling factor (except for the fact that the posterior distribution
 *      reported in the paper cannot be directly translated to the present mod-file due to parameter correlation in the posterior.  
 *  - As pointed out in Del Negro/Schorfheide (2012): "Notes on New-Keynesian Models"
 *      in the code implementation of equation (8) for both the flex price and the sticky price/wage economy, 
 *      there is a (1/(1+cbetabar*cgamma)) missing in the i_2 in front of q_t (denoted qs in the code). 
 *      Equation (8) in the paper reads:  
 *          (1-(1-delta)/gamma)*(1+beta*gamma^(1-sigma))*gamma^2*varphi
 *      which translates to the code snippet:
 *          (1-(1-ctou)/cgamma)*(1+cbetabar*cgamma)*cgamma^2*csadjcost
 *      But the code implements
 *          (1-(1-ctou)/cgamma)*cgamma^2*csadjcost
 *      which corresponds to an equation reading
 *          (1-(1-delta)/gamma)*gamma^2*varphi
 *  - Chib/Ramamurthy (2010): "Tailored randomized block MCMC methods with application to DSGE models", Journal of Econometrics, 155, pp. 19-38
 *      have pointed out that the mode reported in the original Smets/Wouters (2007) paper is not actually the mode. \bar \pi (constepinf) is estimated lower
 *      while \bar \l (constelab) is higher.
 *  - Note that at the prior mean, [cmap,crhopinf] and [cmaw,crhow] are pairwise collinear. Thus, running identification at the prior
 *      mean will return a warning. But this is only a local issue. These parameters are only indistinguishable at the prior mean, but not 
 *      at different points.
 *  - In the prior Table 1A in the paper, the 
 *          - habit parameter $\lambda$ is erroneously labeled h
 *          - the fixed cost parameter $\phi_p$ is labeled $\Phi$ 
 *  - Table 1B claims that $\rho_{ga}$ follows a beta prior with B(0.5,0.2^2), but the code shows that it actually
 *      follows a normal distribution with N(0.5,0.25^2)
 *
 * This file was originally written by Frank Smets and Rafeal Wouters and has been updated by
 * Johannes Pfeifer. 
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 2007-2013 Frank Smets and Raf Wouters
 * Copyright (C) 2013-15 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */

var labobs      ${lHOURS}$      (long_name='log hours worked') 
    robs        ${FEDFUNDS}$    (long_name='Federal funds rate') 
    pinfobs     ${dlP}$         (long_name='Inflation') 
    dy          ${dlGDP}$       (long_name='Output growth rate') 
    dc          ${dlCONS}$      (long_name='Consumption growth rate') 
    dinve       ${dlINV}$       (long_name='Investment growth rate') 
    dw          ${dlWAG}$       (long_name='Wage growth rate') 
    ewma        ${\eta^{w,aux}}$ (long_name='Auxiliary wage markup moving average variable')  
    epinfma     ${\eta^{p,aux}}$ (long_name='Auxiliary price markup moving average variable')  
    zcapf       ${z^{flex}}$    (long_name='Capital utilization rate flex price economy') 
    rkf         ${r^{k,flex}}$  (long_name='rental rate of capital flex price economy') 
    kf          ${k^{s,flex}}$  (long_name='Capital services flex price economy') 
    pkf         ${q^{flex}}$    (long_name='real value of existing capital stock flex price economy')  
    cf          ${c^{flex}}$    (long_name='Consumption flex price economy') 
    invef       ${i^{flex}}$    (long_name='Investment flex price economy') 
    yf          ${y^{flex}}$    (long_name='Output flex price economy') 
    labf        ${l^{flex}}$    (long_name='hours worked flex price economy') 
    wf          ${w^{flex}}$    (long_name='real wage flex price economy') 
    rrf         ${r^{flex}}$    (long_name='real interest rate flex price economy')
    mc          ${\mu_p}$       (long_name='gross price markup') 
    zcap        ${z}$           (long_name='Capital utilization rate') 
    rk          ${r^{k}}$       (long_name='rental rate of capital') 
    k           ${k^{s}}$       (long_name='Capital services') 
    pk          ${q}$           (long_name='real value of existing capital stock') 
    c           ${c}$           (long_name='Consumption')
    inve        ${i}$           (long_name='Investment')
    y           ${y}$           (long_name='Output')
    lab         ${l}$           (long_name='hours worked')
    pinf        ${\pi}$         (long_name='Inflation')
    w           ${w}$           (long_name='real wage')
    r           ${r}$           (long_name='nominal interest rate')
    a           ${\varepsilon_a}$       (long_name='productivity process')
    b           ${c_2*\varepsilon_t^b}$ (long_name='Scaled risk premium shock')
    g           ${\varepsilon^g}$       (long_name='Exogenous spending')
    qs          ${\varepsilon^i}$       (long_name='Investment-specific technology')
    ms          ${\varepsilon^r}$       (long_name='Monetary policy shock process') 
    spinf       ${\varepsilon^p}$       (long_name='Price markup shock process')
    sw          ${\varepsilon^w}$       (long_name='Wage markup shock process')
    kpf         ${k^{flex}}$            (long_name='Capital stock flex price economy') 
    kp          ${k}$           (long_name='Capital stock') 
    ;    
 
varexo ea       ${\eta^a}$      (long_name='productivity shock')
    eb          ${\eta^b}$      (long_name='risk premium shock')
    eg          ${\eta^g}$      (long_name='Spending shock')
    eqs         ${\eta^i}$      (long_name='Investment-specific technology shock')
    em          ${\eta^m}$      (long_name='Monetary policy shock')
    epinf       ${\eta^{p}}$    (long_name='Price markup shock')  
    ew          ${\eta^{w}}$    (long_name='Wage markup shock')  
        ;  
 
parameters curvw ${\varepsilon_w}$  (long_name='Curvature Kimball aggregator wages')  
    cgy         ${\rho_{ga}}$       (long_name='Feedback technology on exogenous spending')  
    curvp       ${\varepsilon_p}$   (long_name='Curvature Kimball aggregator prices')  
    constelab   ${\bar l}$          (long_name='steady state hours')  
    constepinf  ${\bar \pi}$        (long_name='steady state inflation rate')  
    constebeta  ${100(\beta^{-1}-1)}$ (long_name='time preference rate in percent')  
    cmaw        ${\mu_w}$           (long_name='coefficient on MA term wage markup')  
    cmap        ${\mu_p}$           (long_name='coefficient on MA term price markup')  
    calfa       ${\alpha}$          (long_name='capital share')  
    czcap       ${\psi}$            (long_name='capacity utilization cost')  
    csadjcost   ${\varphi}$         (long_name='investment adjustment cost')  
    ctou        ${\delta}$          (long_name='depreciation rate')  
    csigma      ${\sigma_c}$        (long_name='risk aversion')  
    chabb       ${\lambda}$         (long_name='external habit degree')  
    ccs         ${d_4}$             (long_name='Unused parameter')  
    cinvs       ${d_3}$             (long_name='Unused parameter')  
    cfc         ${\phi_p}$          (long_name='fixed cost share')  
    cindw       ${\iota_w}$         (long_name='Indexation to past wages')  
    cprobw      ${\xi_w}$           (long_name='Calvo parameter wages')   
    cindp       ${\iota_p}$         (long_name='Indexation to past prices')  
    cprobp      ${\xi_p}$           (long_name='Calvo parameter prices')   
    csigl       ${\sigma_l}$        (long_name='Frisch elasticity')   
    clandaw     ${\phi_w}$          (long_name='Gross markup wages')   
    crdpi       ${r_{\Delta \pi}}$  (long_name='Unused parameter')  
    
    crhoa       ${\rho_a}$          (long_name='persistence productivity shock')  
    crhoas      ${d_2}$             (long_name='Unused parameter')  
    crhob       ${\rho_b}$          (long_name='persistence risk premium shock')  
    crhog       ${\rho_g}$          (long_name='persistence spending shock')  
    crhols      ${d_1}$             (long_name='Unused parameter')  
    crhoqs      ${\rho_i}$          (long_name='persistence risk premium shock')  
    crhoms      ${\rho_r}$          (long_name='persistence monetary policy shock')  
    crhopinf    ${\rho_p}$          (long_name='persistence price markup shock')  
    crhow       ${\rho_w}$          (long_name='persistence wage markup shock')  
    ctrend      ${\bar \gamma}$     (long_name='net growth rate in percent')  
    cg          ${\frac{\bar g}{\bar y}}$     (long_name='steady state exogenous spending share')  
    
    phia        $\phi_a$            (long_name='real rate loading on a')
    phib        $\phi_b$            (long_name='real rate loading on b')
    phig        $\phi_g$            (long_name='real rate loading on g')
    phiqs       $\phi_{i^s}$           (long_name='real rate loading on qs')
    phispinf    $\phi_{p^inf}$         (long_name='real rate loading on spinf')
    phisw       $\phi_{w^s}$            (long_name='real rate loading on sw')  
    
    phiewma     $\phi_{\eta^{w}}$
    phiepinfma  $\phi_{\eta^{p}}$
    phiyf       $\phi_{y^{p}}$
    phiy        $\phi_y$
    phir        $\phi_r$
    phikpf      $\phi_{k^{p,s}}$
    phikp       $\phi_{k^{s}}$
    phicf       $\phi_{c^{p}}$
    phiinvef    $\phi_{i^{p}}$
    phic        $\phi_c$
    phiinve     $\phi_i$
    phipinf     $\phi_{\pi}$ 
    phiw        $\phi_w$

    alphac     $\alpha_c$           (long_name='Euler discounting parameter')
    
    ;

               
// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

% FPJ/AP
alphac=.999;
% alphac=1;

// estimated parameters initialisation
calfa=0.211;
cbeta=.997944;
csigma=1.607;
cfc=1.613;
cgy=0.575;

csadjcost= 5.856;
chabb=    0.710;    
cprobw=   0.713;
csigl=    1.155;
cprobp=   0.678;
cindw=    0.608;
cindp=    0.178;
czcap=    0.737;


crhoa=    0.982;
crhob=    0.275;
crhog=    0.986;
crhols=   0.9;
crhoqs=   0.809;
crhoas=1; 
crhoms=0.169;
crhopinf=0.965;
crhow=0.970;
cmap = 0.807;
cmaw  = 0.901;

constelab= -1.316;

% FPJ/AP
phia    =   -0.023;
phib    =   0.576;
phig    =   0.067;
phiqs   =   0.200;
phispinf =   -0.426;
phisw    =   -0.234;

phiewma = -0.027;
phiepinfma = -3.055;
phiyf = 2.289;
phiy = -2.408;
phir = 2.840;
phikpf = -0.923;
phikp = 0.034;
phicf = 1.986;
phiinvef = -1.593;
phic = 0.559;
phiinve = 1.094;
phipinf = 1.580;
phiw  = -0.396;

%% Added by JP to provide full calibration of model before estimation
constepinf=0.609;
constebeta=0.160;
ctrend=0.385;

model(linear); 
//deal with parameter dependencies; taken from usmodel_stst.mod 
#cpie=1+constepinf/100;         %gross inflation rate
#cgamma=1+ctrend/100 ;          %gross growth rate
#cbeta=1/(1+constebeta/100);    %discount factor

#clandap=cfc;                   %fixed cost share/gross price markup
#cbetabar=cbeta*cgamma^(-csigma);   %growth-adjusted discount factor in Euler equation
#cr=cpie/(cbeta*cgamma^(-csigma));  %steady state net real interest rate
#crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou); %R^k_{*}: steady state rental rate
#cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));      %steady state real wage
//cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
#cikbar=(1-(1-ctou)/cgamma);        %k_1 in equation LOM capital, equation (8)
#cik=(1-(1-ctou)/cgamma)*cgamma;    %i_k: investment-capital ratio
#clk=((1-calfa)/calfa)*(crk/cw);    %labor to capital ratio
#cky=cfc*(clk)^(calfa-1);           %k_y: steady state output ratio
#ciy=cik*cky;                       %investment-output ratio
#ccy=1-cg-cik*cky;                  %consumption-output ratio
#crkky=crk*cky;                     %z_y=R_{*}^k*k_y
#cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy; %W^{h}_{*}*L_{*}/C_{*} used in c_2 in equation (2)
#cwly=1-crk*cky;                    %unused parameter

#conster=(cr-1)*100;                %steady state federal funds rate ($\bar r$)

// flexible economy

          [name='FOC labor with mpl expressed as function of rk and w, flex price economy']
	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      [name='FOC capacity utilization, flex price economy']
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
          [name='Firm FOC capital, flex price economy']
	      rkf =  (wf)+labf-kf ;
          [name='Definition capital services, flex price economy']
	      kf =  kpf(-1)+zcapf ;
% FPJ
          [name='Investment Euler Equation, flex price economy']
% 	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
	      invef = alphac*((1/(1+cbetabar*cgamma))* (invef(-1) + cbetabar*cgamma*invef(1))+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          [name='Arbitrage equation value of capital, flex price economy']
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
% FPJ
	      [name='Consumption Euler Equation, flex price economy']
%         cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
          cf = alphac*((chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1)) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      [name='Aggregate Resource Constraint, flex price economy']
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      [name='Aggregate Production Function, flex price economy']
          yf = cfc*( calfa*kf+(1-calfa)*labf +a );
          [name='Wage equation, flex price economy']
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
          [name='Law of motion for capital, flex price economy (see header notes)']              
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy
          [name='FOC labor with mpl expressed as function of rk and w, SW Equation (9)']
	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      [name='FOC capacity utilization, SW Equation (7)']
          zcap =  (1/(czcap/(1-czcap)))* rk ;
          [name='Firm FOC capital, SW Equation (11)']
	      rk =  w+lab-k ;
          [name='Definition capital services, SW Equation (6)']
	      k =  kp(-1)+zcap ;
% FPJ/AP
          [name='Investment Euler Equation, SW Equation (3)']
% 	      inve = (1/(1+cbetabar*cgamma))* (inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
	      inve = alphac*(1/(1+cbetabar*cgamma))* (inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          [name='Arbitrage equation value of capital, SW Equation (4)']
          pk = -r+pinf(1)-0*b 
                + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b 
                + (crk/(crk+(1-ctou)))*rk(1) 
                + ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
% FPJ/AP
	      [name='Consumption Euler Equation, SW Equation (2)']
% 	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) 
%             + (1/(1+chabb/cgamma))*c(+1) 
%             +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) 
%            - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
          c = alphac*((chabb/cgamma)/(1+chabb/cgamma)*c(-1) 
                + (1/(1+chabb/cgamma))*c(+1)) 
                +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) 
                - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;          
	      [name='Aggregate Resource Constraint, SW Equation (1)']
          y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
          [name='Aggregate Production Function, SW Equation (5)']
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
          [name='New Keynesian Phillips Curve, SW Equation (10)']
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
	      [name='Wage Phillips Curve, SW Equation (13), with (12) plugged for mu_w']
          w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;
% FPJ/AP
%           [name='Taylor rule, SW Equation (14)']
% 	      r =  crpi*(1-crr)*pinf
%                +cry*(1-crr)*(y-yf)     
%                +crdy*(y-yf-y(-1)+yf(-1))
%                +crr*r(-1)
%                +ms  ;
          [name='Real Rate rule ']
	      r -  pinf(+1) = 
               phia*a     
               +phib*b     
               +phig*g     
               +phiqs*qs     
               +phispinf*spinf     
               +phisw*sw

               +phiewma*ewma(-1)
               +phiepinfma*epinfma(-1)
               +phiyf*yf(-1)
               +phiy*y(-1)
               +phir*r(-1) 
               +phikpf*kpf
               +phikp*kp
               +phicf*cf(-1)
               +phiinvef*invef(-1)
               +phic*c(-1)
               +phiinve*inve(-1)
               +phipinf*pinf(-1)
               +phiw*w(-1) 
             
               +ms ;
              
          [name='Law of motion for productivity']              
	      a = crhoa*a(-1)  + ea;
          [name='Law of motion for risk premium']              
	      b = crhob*b(-1) + eb;
          [name='Law of motion for spending process']              
	      g = crhog*(g(-1)) + eg + cgy*ea;
	      [name='Law of motion for investment specific technology shock process']              
          qs = crhoqs*qs(-1) + eqs;
          [name='Law of motion for monetary policy shock process']              
	      ms = crhoms*ms(-1) + em;
          [name='Law of motion for price markup shock process']              
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
          [name='Law of motion for wage markup shock process']              
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew; 
          [name='Law of motion for capital, SW Equation (8) (see header notes)']              
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

// measurement equations
[name='Observation equation output']              
dy=y-y(-1)+ctrend;
[name='Observation equation consumption']              
dc=c-c(-1)+ctrend;
[name='Observation equation investment']              
dinve=inve-inve(-1)+ctrend;
[name='Observation equation real wage']              
dw=w-w(-1)+ctrend;
[name='Observation equation inflation']              
pinfobs = 1*(pinf) + constepinf;
[name='Observation equation interest rate']              
robs =    1*(r) + conster;
[name='Observation equation hours worked']              
labobs = lab + constelab;

end; 

steady_state_model;
dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs = constepinf;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
labobs = constelab;
end;

steady;
check;

% return

shocks;
var ea;
stderr 0.461;
var eb;
stderr 0.217;
var eg;
stderr 0.475;
var eqs;
stderr 0.338;
var em;
stderr 0;
var epinf;
stderr 0.133;
var ew;
stderr 0.340;
end;




stoch_simul(order = 1, irf=31,nograph,periods=100000) y r pinf yf g;


L = 4;
y_yf = y - yf;

Y4 = r(2+L:end);
X4 = [r(1+L:end-1) y_yf(2+L:end) (y_yf(2+L:end)-y_yf(1+L:end-1)) pinf(2+L:end)];
Z4 = [r(1:end-(L+1)) y_yf(2:end-L) (y_yf(2:end-L)-y_yf(1:end-(L+1))) pinf(2:end-L)];

Beta_hat =  inv(X4'*X4)*X4'*Y4;
Phi_r = Beta_hat(1)
Phi_y_yf = Beta_hat(2)/(1-Phi_r)
Phi_y_yf_diff = Beta_hat(3)
Phi_Pinf = Beta_hat(4)/(1-Phi_r)


save SW2007_RR_Full_Rule_sim y r pinf yf g;



