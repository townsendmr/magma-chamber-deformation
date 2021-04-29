function [eps_x, deps_x_dT,deps_x_deps_g]= crystal_fraction(T,eps_g)

% the b value controls shape of melting curve (0.5 approximates a dacite, 
% increase b for more mafic, decrease for more silicic, but don't 
% forget to also change the solidus and liquidus)
b              = 0.1; 

T_sol          = 700 + 273;
T_liq          = 950 + 273;


phi_x          = 1 - ((T-T_sol)/(T_liq-T_sol))^b;
dphi_x_dT      = - b*(T-T_sol)^(b-1)/(T_liq-T_sol)^b;
if T<T_sol
    phi_x=1;
    dphi_x_dT=0;
elseif T>T_liq
    phi_x=0;
    dphi_x_dT=0;
end
eps_x          = (1-eps_g)*phi_x;
deps_x_dT      = (1-eps_g)*dphi_x_dT;
deps_x_deps_g  = -phi_x;
% eps_x          = 0.2;
% deps_x_dT      = 0;
% deps_x_deps_g  = 0;

