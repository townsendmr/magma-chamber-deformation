function [rho_g,drho_g_dP,drho_g_dT] = eos_g(P,T)

% % ideal gas law
% gasConst   = 462;       % specific gas constant for water (J kg^-1 K^-1)
% 
% rho_g       = P/gasConst/T;
% drho_g_dP   = 1/gasConst/T;
% drho_g_dT   = -P/gasConst/T^2;

% rho_g       = 2400;%P/gasConst/T;
% drho_g_dP   = 0;%1/gasConst/T;
% drho_g_dT   = 0;%-P/gasConst/T^2;

% parametrization of redlich kwong taken from Huber et al. 2010
rho_g        = -112.528*(T-273.15)^-0.381 + 127.811*(P*1e-5)^-1.135 + 112.04*(T-273.15)^-0.411*(P*1e-5)^0.033; 
drho_g_dP    = (-1.135)*127.811*(P*1e-5)^-2.135 + 0.033*112.04*(T-273.15)^-0.411*(P*1e-5)^-0.967;
drho_g_dT    = (-0.381)*(-112.528)*(T-273.15)^-1.381 + (-0.411)*112.04*(T-273.15)^-1.411*(P*1e-5)^0.033;

rho_g        = rho_g*1e3; %display(rho_g)
drho_g_dP    = drho_g_dP*1e-2;%display(drho_gdP)
drho_g_dT    = drho_g_dT*1e3;%display(drho_gdT)