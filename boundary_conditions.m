function [Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out, Hdot_in, Hdot_out, P_loss] ...
    = boundary_conditions(sw,param,P,T,eps_g,eps_x,V,rho_m,rho_x,m_eq,rho_g,rho,...
    c, m_eq_in, T_in, mdot_in,time)



 global  Mdot_out_pass


 
 
% set inflow conditions
rho_m_in       = rho_m; %2400
rho_x_in       = rho_x; %2600;

P_in           = P;
%T_in           = 1200;
eps_g_in       = 0.0;
[eps_x_in, ~,~] = crystal_fraction(T_in,eps_g_in);
[rho_g_in,~,~]  = eos_g(P_in,T_in);

rho_in         = (1-eps_g_in-eps_x_in)*rho_m_in + eps_g_in*rho_g_in + eps_x_in*rho_x_in;
c_in           = ((1-eps_g_in-eps_x_in)*rho_m_in*param.c_m + eps_g_in*rho_g_in*param.c_g + eps_x_in*rho_x_in*param.c_x)/rho_in;%c;

    
Mdot_in        = mdot_in;
Mdot_v_in      = m_eq_in*rho_m_in*(1-eps_g_in-eps_x_in)*Mdot_in/rho_in + rho_g_in*eps_g_in*Mdot_in/rho_in;
Hdot_in        = c_in*T_in*Mdot_in;

% would have been better if this was passed from the main file
% global P_0_pass % MRT
% P_lit=P_0_pass; % MRT
  %P_lit = 100e6;
% set outflow conditions
if sw.eruption ==0
    Mdot_out       = 0;
elseif sw.eruption==1
    Mdot_out       = Mdot_out_pass;
else
    disp('eruption not specified')
end

% outer shell conditions
a                    = (V/(4*pi/3))^(1/3);   % chamber radius (m)
cc                   = 10*a; %outer shell radius (m)
dr                   = 0.1*a;

% % material properties
% kappa                = 1e-6;  % thermal conductivity (W m^-1 K^-1)
% rho_r                = 2500;  % density crust (kg/m3)
% cp                   = 1200;  % heat capacity crust (J kg^-1 K-1)
%Tb                   = 500;   %
% Tb = Tb;   % MRT


% precision
maxn  = 10000; % number of grid points for temeparture profile in the surrounding shell
%, increases precision, slows down code

Mdot_v_out     = m_eq*rho_m*(1-eps_g-eps_x)*Mdot_out/rho + rho_g*eps_g*Mdot_out/rho;
if sw.heat_cond==1

    % simplified version
%     dr                   = 10000;
%     small_q              = -kappa*rho_r*cp*(500-T)/(2*a);
%     surface_area_chamber = 4*pi*a^2;
%     Q_out                = small_q*surface_area_chamber;
    
    
    % heat loss
    Q_out          = heat_conduction_chamber(maxn,a,cc,dr,param.kappa,param.rho_r,param.c_r,param.Tb); 

    %Q_out          = heat_conduction_chamber(maxn,a,cc,dr,param.kappa,param.rho_r,param.c_r,param.Tb); 
elseif sw.heat_cond==0
    Q_out          = 0;
else
    disp('heat_cond not specified')
end

if isnan(Q_out)
    Q_out=0;
    disp('Q_out is NaN')
end

Hdot_out       = c*T*Mdot_out + Q_out;


% viscous relaxation
% nn   = 1.9;   % power law exponent
% AA   = 2e-4*(1e6)^-nn;  % material dependent constant (Pa^-n s^-1)
% G    = 141e3; % activation energy for creep (J/mol)
% M    = 8.314472;   % Molar gas constant (which gas?) (J/mol/K = kg m^2 s^-2 mol^-1 K^-1)
    GLQ_n =64;
    [quadpts,weights] = GLQ_points_weights_hard(GLQ_n);
    b     = a +cc;
    quadpts_r = (b-a)/2.*quadpts + (a+b)/2;
    
  if sw.heat_cond==1  
    
    Trt = heat_conduction_chamber_profile(maxn,a,cc,quadpts_r,param.kappa,param.Tb);
         A = 4.25e7; % material-dependent constant for viscosity law (Pa s)
B = 8.31; % molar gas constant (J/mol/K)
G = 141e3; % activation energy for creep (J/mol)
%dev_stress = 20e6;
%eta_rt     = (dev_stress.^(1-nn)./AA).*exp(G./M./Trt);
eta_rt     = A.*exp(G./B./Trt);
     I          = (b-a)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
     eta_r      = 3*a^3*I;
 
  else
      eta_rt = param.eta_rt0;
      Trt = param.Trt0;
      quadpts_r = param.quadpts_r0;
     I          = (b-a)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
     eta_r      = 3*a^3*I;
  end
      


     
eta_rt = sort(eta_rt);
quadpts_r = sort(quadpts_r);
Trt = sort(Trt, 'descend');




if sw.visc_relax ==1
    
    %P_loss =  (P - param.P_lit)/eta_r    
    if strcmp(param.viscModel, 'shell')    
    P_loss = viscous_P_loss_shell(quadpts_r,Trt,param);
    elseif strcmp(param.viscModel, 'full') 
        P_loss = viscous_P_loss_full(eta_r,param);
    end
elseif sw.visc_relax==0
    P_loss = 0;
else
    disp('visc_relax not specified')
end
