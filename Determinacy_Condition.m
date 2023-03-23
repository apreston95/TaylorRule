function output_cell =posterior_function_demo(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info)
% output_cell =posterior_function_demo(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info);
% This is an example file computing statistics on the prior/posterior draws. The
% function allows read-only access to all Dynare structures. However, those
% structures are local to this function.  Changing them will not affect
% other Dynare functions and you cannot use them to pass results to other
% Dynare functions.
% The function takes one and only one output argument: an 1 by n cell.
% Using functions like cell2mat, the contents of the cell can be easily
% transformed back to matrices. See the fs2000_posterior_function.mod for
% an example

% INPUTS
%   xparam1                      Current parameter draw
%   M_           [structure]     Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%   options_     [structure]     Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_          [structure]     Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%   estim_params_[structure]     Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   bayestopt_   [structure]     Matlab's structure describing the parameter options (initialized by dynare, see @ref{bayestopt_}).
%   dataset_     [structure]     Matlab's structure storing the dataset
%   dataset_info [structure]     Matlab's structure storing the information about the dataset

% Output
%   output_cell  [1 by n cell]   1 by n Matlab cell allowing to store any
%                                desired computation or result (strings, matrices, structures, etc.)

% Copyright (C) 2015 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


%% store the slope based on the parameter draw
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names{ii,:});
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

sigmad = 1;
alphar = 1/Gamma;
alphay = Alpha_Y;
rhod = Rho_D;
phimu = Phi_Z;
rhomu = Rho_Z;
sigmanu = 1;
sigmamu = 1;
rhonu = Rho_R;
beta = Beta;
gammay = Kappa;
gammar = 0;
phid = Phi_D;



BBt(3,1)=(sigmad-alphar*phid)/(1-alphay*rhod);
BBt(3,2)=(-alphar*phimu)/(1-alphay*rhomu);
BBt(3,3)=-alphar*sigmanu/(1-alphay*rhonu);
BBt(1,1)=1/(1-beta*rhod)*(gammay*BBt(3,1)+gammar*phid);
BBt(1,2)=1/(1-beta*rhomu)*(sigmamu+gammay*BBt(3,2)+gammar*phimu);
BBt(1,3)=1/(1-beta*rhonu)*(gammay*BBt(3,3)+gammar*sigmanu);
BBt(2,1)=phid+rhod*BBt(1,1);
BBt(2,2)=phimu+rhomu*BBt(1,2);
BBt(2,3)=sigmanu+rhonu*BBt(1,3);

A=zeros(3,3);
B=BBt;
R=[rhod 0 0;0 rhomu 0 ;0 0 rhonu];

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

phipi       =   (b2-c2*b1/c1)/(a2-a1*c2/c1);
phiy        =   (b1-a1*phipi)/c1;

output_cell{1,1}=phiy*((1-Beta)/Kappa) + phipi;

end