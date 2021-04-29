%%% This code links the magma chamber box model of 
%%% Degruyter and Huber (2014) with models for viscoelastic crust and
%%% surface deformation. The two viscoelastic models are from 
%%% Bonafede and Ferrari (2009) for a fully viscoelastic half space
%%% and Segall (2010, ch 7) for a viscoelastic shell in a half space
%%% Model inputs include: chamber size, depth, temperature, water content,
%%% recharge conditions, and viscoelastic parameters.
%%% Model outputs include: time series of pressure, temperature, volume,
%%% mass, crystal fraction, gas fraction, and surface deformation
%%% Please cite both the original magma chamber model paper (Degruyter and
%%% Huber, 2014, EPSL) as well as the paper presenting the link with 
%%% surface deformation (Townsend, 2021, EPSL)

clear all
close all
clc
tic

% create structure for switches
sw.heat_cond  = 1; % switch cooling module on/off
sw.visc_relax = 1; % switch viscous relaxation on/off
sw.plot_recharge = 0; % switch for plotting recharge rate over time

%% Choose Viscoelastic Model
% 'shell' uses the viscoelastic shell model of Segall (2010) Chapter 7 and
% lets the user choose either to prescribe the shell radius and viscosity
% or set a cutoff viscosity/timescale
% 'full' treats the entire crust as viscoelastic using model of Bonafede
% and Ferrari (2009) and calculates effective viscosity using model of
% Lensky et al., (2001)

% param.viscModel = 'shell';
 param.viscModel = 'full';

% % For shell model, determine shell radius and viscosity
% param.shellOverride=1; % override calculation of shell parameters?
% param.R_shell = 1500; % Radius of viscoelastic shell (m)
% param.eta_shell = 4.1995e+17; % Viscosity of viscoelastic shell (Pa s)
% cutoff to determine shell radius and viscosity
% param.shellOverride = 0; % override calculation of shell parameters?
% param.tObs = 3600*24*365*30; % timescale of observation (seconds), to calculate thickness of viscoelastic shell;

% simulation time
end_time_yr   = 10000;
end_time      = 3600*24*365*end_time_yr;   % maximum simulation time in years

%% Set magma recharge style and parameters
% 'constant' uses a constant rate of recharge
% 'pulse' uses a Gaussian-shaped distribtuion of recharge
% '2chamber' uses recharge from a connected deeper reservoir

param.recharge = 'constant';
rate_km3_per_yr = 0.0023; % volumetric flow rate in km^3 / year
param.rate     = (rate_km3_per_yr*1e9*2400)/(365*24*3600); % Recharge rate (kg/s)

% param.recharge  = 'pulse';
% param.peak_time = 3600*24*365*2.5; % time of first peak recharge (s)
% param.sig_time  = 3600*24*365*(5/7); % parameter that controls duration (s)
% param.mass_injected  = 1.2*10^11; %2.4*10^10; % total mass injected during single recharge event (kg)
% param.repeat    = 0; % repeat recharge events?
% param.periodicity = 3600*24*365*25; % periodicity of repeated recharge events (s)

% disp('peak recharge rate (kg/s)')
% disp(param.mass_injected/(param.sig_time*sqrt(2*pi))) % max inflow rate (kg/s)
% disp('volume injected per recharge event (m^3)')
% disp(param.mass_injected/2400)
% disp('duration of event (approximate, years)')
% disp(param.sig_time*7./(3600*24*365))
% disp('ave volumetric flow rate of single event (km^3/year)')
% disp((param.mass_injected/2400/1e9)/(param.sig_time*7./(3600*24*365)))

% param.recharge  = '2chamber';
% param.cond      = 1e-4; % Hydraulic conductivity between deep and shallow reservoirs (m s)
% param.P_deep    = 300e6; % Constant pressure of the deep reservoir (Pa)

%% Set other parameters

% error tolerances used in ode method
rel_tol        = 1e-6;
abs_tol        = 1e-6;

% eruptions
max_eruptions = 65;  % maximum number of eruptions
P_crit        = 20e6; % critical overpressure for eruption (Pa)

% magma constants
param.beta_m         = 1e10;      % bulk moduluis melt (Pa)
param.beta_x         = 1e10;      % bulk moduluis crystals (Pa)
param.alpha_m        = 1e-5;      % thermal expansion melt (K^-1)
param.alpha_x        = 1e-5;      % thermal expansion crystals (K^-1)
param.L_e            = 610e3;     % latent heat of exsolution (J kg^-1)
param.L_m            = 290e3;     % latent heat of melting (J kg^-1)
param.c_m            = 1200;      % specific heat of melt (J kg^-1 K-1)
param.c_g            = 3880;      % specific heat of gas (J kg^-1 K-1)
param.c_x            = 1200;      % specific heat of crystals (J kg^-1 K-1)

% crustal properties
param.alpha_r        = 1e-5;      % thermal expansion country rock (K^-1)
param.beta_r         = 1e10;      % bulk modulus country rock (Pa)
param.kappa          = 1e-6;      % thermal conductivity (W m^-1 K^-1)
param.rho_r          = 2550;      % density crust (kg/m3)
param.c_r            = 1200;      % heat capacity crust (J kg^-1 K-1)
param.pr_r           = 0.25;      % Poisson's Ratio
param.mu_r           = 3*param.beta_r*(1-2*param.pr_r)/(2*(1+param.pr_r)); % Shear modulus of crust (Pa)

% bm = 1/param.beta_m;
% bc = 1/param.beta_r;
% B_s = bc/(bm+bc)

% chamber radius (m)
log_volume_km3 = log10(20);  % log volume in km3
volume_km3 = 10.^log_volume_km3; % volume in km3
a     = 1000.*(volume_km3./(4*pi/3)).^(1/3); % radius in m

% water content
param.water      = 0.06;%0.04:0.01:0.07;%:0.01:0.05;

% storage depth (m)
param.depth      = 8e3; %6e3:1e3:12e3;

%% Run chamber code
    
    % initial conditions

    % time
    begin_time    = 0;       % initialize time
    %end_time      = param.tObs;      % maximum simulation time in years
    %end_time      = end_time*365*24*3600;      % maximum simulation time in seconds
    
    % thermal gradient
    T_surface     = 0+273; % surface temperature (K)
    T_gradient    = 32/1e3; % thermal gradient (K/m)
    param.Tb      = T_surface+T_gradient*param.depth; % background temperature crust (K)
    
    %lithostatic pressure
    global P_lit
    grav_acc      = 9.81; % gravitational acceleration (m/s2)
    P_lit         = param.rho_r*grav_acc*param.depth; % Lithostatic pressure (Pa)
    param.P_lit   = P_lit;
    P_0           = P_lit; % initial chamber pressure (Pa)
    
    % magma chamber initial conditions
    rho_m0        = 2400;       % initial melt density (kg/m^3)
    rho_x0        = 2600;       % initial crystal density (kg/m^3)
    V_0           = 4*pi/3*a^3; % initial volume of the chamber (m^3)
    T_0           = 1100;       % initial chamber temperature (K)
    
    % recharging magma properties
    T_in          = 1200;       % Temperature of recharging magma (K)
    m_eq_in       = param.water;
    
    % eruption criteria   
    mdot_out      = 1e5;        % outflow during eruption (kg/s)    
    
    mainChamber % run mainChamber.m
    surf_disp   % run surf_disp.m
toc
%% Plotting (only for one case)

%Pressure
figure
plot(time./(3600*24*365),P./1e6,'LineWidth',2), hold on
xlabel('time (yr)','FontSize',16), ylabel('pressure (MPa)','FontSize',16)
set(gca,'FontSize',16)

% Gas volume fraction
figure
plot(time./(3600*24*365),eps_g,'LineWidth',2),hold on
xlabel('time (yr)','FontSize',16), ylabel('gas volume fraction','FontSize',16)
set(gca,'FontSize',16)

% Crystal volume fraction
figure
plot(time./(3600*24*365),eps_x,'LineWidth',2),hold on
xlabel('time (yr)','FontSize',16), ylabel('crystal volume fraction','FontSize',16)
set(gca,'FontSize',16)

% Temperature
figure
plot(time./(3600*24*365),T,'LineWidth',2),hold on
xlabel('time (yr)','FontSize',16), ylabel('temperature (K)','FontSize',16)
set(gca,'FontSize',16)

% Volume change
figure
plot(time./(3600*24*365),V./V(1),'LineWidth',2), hold on
xlabel('time (yr)','FontSize',16), ylabel('volume/V_0','FontSize',16)
set(gca,'FontSize',16)

% Mass change
figure
plot(time./(3600*24*365),(rho.*V)./(rho(1)*V(1)),'LineWidth',2), hold on
xlabel('time (yr)','FontSize',16), ylabel('mass/M_0','FontSize',16)
set(gca,'FontSize',16)

% Surface uplift
figure
plot(time(1:length(uz_s))./(365*24*3600),uz_s,'LineWidth', 2)
xlabel('time in years');  ylabel('displacement (m)')
set(gca,'FontSize',16)

