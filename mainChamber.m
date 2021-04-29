% !!!!!!!!!!!!!
% USE runCode.m
% !!!!!!!!!!!!!

% calculate additional initial conditions
[m_eq0,~,~]   = exsolve(P_0,T_0); % initial dissolved water mass fraction
[rho_g0,~,~]  = eos_g(P_0,T_0); % initial gas density

% Gas fraction
% initial estimate of crystal fraction (no gas)
[eps_x0, ~,~] = crystal_fraction(T_0,0);
% initial estimate of melt fraction (no gas)
eps_m0 = 1-eps_x0;
% initial estimate of bulk density (kg/m^3) (no gas)
rho_0 = (1-eps_x0)*rho_m0 + eps_x0*rho_x0;
% initial estimate of water concentration in melt
wt_g = m_eq_in*rho_0/(rho_m0*eps_m0);

% get a better estimate of gas volume fraction
if wt_g > m_eq0
    eps_g_guess = 0.0001;
    [eps_x0, ~,~] = crystal_fraction(T_0,eps_g_guess);
    wt_g_est = (rho_g0/rho_m0)*(eps_g_guess/(1-eps_g_guess-eps_x0)) + m_eq0;
    while wt_g_est < wt_g
        eps_g_guess = eps_g_guess + 0.0001;
        [eps_x0, ~,~] = crystal_fraction(T_0,eps_g_guess);
        wt_g_est = (rho_g0/rho_m0)*(eps_g_guess/(1-eps_g_guess-eps_x0)) + m_eq0;
    end
    eps_g0 = eps_g_guess;
else
    eps_g0 = 0;
end
% update initial crystal volume fraction
[eps_x0, ~,~] = crystal_fraction(T_0,eps_g0);
% update initial melt volume fraction
eps_m0 = 1-eps_x0-eps_g0;
% update initial bulk density (kg/m^3)
rho_0  = (1-eps_g0-eps_x0)*rho_m0 + eps_g0*rho_g0 + eps_x0*rho_x0;

% track total volatile mass
if eps_g0 > 0
    mf_ex0 = eps_g0*V_0*rho_g0;
    mf_dis0 = m_eq0*eps_m0*V_0*rho_m0;
    param.mf_0 = mf_ex0 + mf_dis0;
    param.conc_0 = m_eq0;
else
    param.mf_0 = wt_g*(1-eps_x0)*V_0*rho_m0;
    param.conc_0 = param.mf_0/((1-eps_x0)*V_0*rho_m0);
end

% Declare and initialize global variables
global DP_crit switch_Tprofile
global storeSumk storeSumk_2 storeSumk_old storeSumk_oldold storeSumk_2_old storeSumk_2_oldold lengthTime maxTime

lengthTime     = 0;
maxTime        = 0;
switch_Tprofile= 0;

global storeTime storeTemp Mdot_out_pass storeP storeV storeMF storeCONC
global ii
ii      = 1;   % initial iteration for while loop
ii_max  = max_eruptions +2;
phase0  = 0;   % number that indicates why the code stopped

% Initial temperature and viscosity profile around chamber
maxn = 10000;
GLQ_n =64;
[quadpts,weights] = GLQ_points_weights_hard(GLQ_n);
cc    = 10*a; %outer radius for heat conduction (m)
b     = a +cc;
quadpts_r = (b-a)/2.*quadpts + (a+b)/2;
storeTime =0;
storeTemp = T_0;
Trt = heat_conduction_chamber_profile(maxn,a,cc,quadpts_r,param.kappa,param.Tb);
A = 4.25e7; % material-dependent constant for viscosity law (Pa s)
B = 8.31; % molar gas constant (J/mol/K)
G = 141e3; % activation energy for creep (J/mol)
eta_rt     = A.*exp(G./B./Trt);

eta_rt = sort(eta_rt);
quadpts_r = sort(quadpts_r);
Trt = sort(Trt, 'descend');

param.Trt0 = Trt;
param.quadpts_r0 = quadpts_r;
param.eta_rt0 = eta_rt;

% figure
% plot((quadpts_r-a)./1e3,Trt,'r')
% xlabel('distance from chamber (km)')
% ylabel('temperature (K)')
% 
% figure
% plot((quadpts_r-a)./1e3,eta_rt,'k')
% xlabel('distance from chamber (km)')
% ylabel('log viscosity (Pa s)')
% set(gca, 'YScale', 'log')

t_R_vec = 0;
for i = 1:length(Trt)
    R2 = quadpts_r(i);
    T_S = Trt(i);
    eta_shell   = calc_viscosity(a,R2, T_0, T_S);
    t_R_vec(i) = (3*eta_shell*(1-param.pr_r)*(R2/a)^3)/(param.mu_r*(1+param.pr_r));
end

% figure
% plot((quadpts_r-a)./1e3,t_R_vec./(3600*24*365),'k')
% xlabel('distance from chamber (km)')
% ylabel('log relaxation time (years)')
% set(gca, 'YScale', 'log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate effective viscosity
if sw.visc_relax == 1
    if strcmp(param.viscModel,'shell')
        if param.shellOverride == 0
            % given an observation timescale, calculate the thickness
            % and effective viscosity of the wall rocks that
            % participate in viscoelastic deformation
            t_R =0; i = 2;
            while t_R<param.tObs
                R2 = quadpts_r(i);
                T_S = Trt(i);
                eta_shell   = calc_viscosity(a,R2, T_0, T_S);
                t_R = (3*eta_shell*(1-param.pr_r)*(R2/a)^3)/(param.mu_r*(1+param.pr_r));
                i = i+1;
            end
            
            param.eta_shell = eta_shell;
            param.R_shell   = R2;
            disp('shell parameters calculated based on timescale')
            disp('viscoelastic shell radius (m)')
            disp(param.R_shell)
            disp('shell effective viscosity (Pa s)')
            disp(param.eta_shell)
            disp('shell relaxation timescale (years)')
            disp(t_R/(3600*24*365))
            
        elseif param.shellOverride == 1
            tau = (3*param.eta_shell*(1-param.pr_r)*(param.R_shell/a)^3)/(param.mu_r*(1+param.pr_r));
            disp('shell parameters prescribed by user')
            disp('crust effective viscosity (Pa s)')
            disp(param.eta_shell)
            disp('shell relaxation timescale (years)')
            disp(tau/(3600*24*365))
        else
            disp('error: shell override must be 0 or 1')
        end
        
        
        
    elseif strcmp(param.viscModel,'full')
        b = 11*a;
        param.eta_shell   = calc_viscosity(a, b, T_0, param.Tb);
        disp('crust effective viscosity (Pa s)')
        disp(param.eta_shell)
    else
        disp('error: must determine viscoelastic model')
        
    end
end


% time scales (s)
% tau_inj            = rho_0*V_0/(strength*(1/(sig_time*sqrt(2*pi)))); % injection timescale
% tau_cooling        = a^2/kappa;        % cooling timescale
% tau_visco_elastic  = eta_r0/P_crit;     % relaxation timescale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % write file for input parameters
% dlmwrite(['input/input' num2str(run_n,'%0.5i') '.dat'],[P_0 ;T_0; eps_g0;...
%     a; rho_m0; rho_x0; peak_time; rho_0; beta_r; ...
%     tau_inj; tau_cooling; tau_visco_elastic; ...
%     P_lit; Tb; P_crit; rho_r; mdot_out; wt_g; m_eq_in; eta_r0;param.depth; sig_time],'precision',10)
%
%


% use global variable
Mdot_out_pass  = mdot_out;

% initialize vector to store quantities
storeTime = [];
storeTemp = [];
storeP = [];
storeV = [];
storeMF = [];
storeCONC = [];
time  = [];
P     = [];
T     = [];
eps_g = [];
V     = [];
rho_m = [];
rho_x = [];
time_erupt = [];
dur_erupt = [];
mass_erupt = [];
vol_erupt = [];

%figure
while ii<ii_max && begin_time < end_time && phase0 ==0
    % set up ode
    tspan       = [begin_time end_time]; % solve from t=0s to t=end_time
    sw.eruption = 0; % switch to turn eruption on/off
    IC          = [P_0 T_0 eps_g0 V_0 rho_m0 rho_x0];   % store initial conditions
    options     = odeset('RelTol',rel_tol,'AbsTol',abs_tol,...
        'Events',@(z,y) stopChamber(y,sw.eruption, param.P_lit, P_crit),'Refine',2);
    
    % solve ode
    [X,Y,TE,YE,IE] = ode15s(@(t,y) odeChamber(t,y,        ...
        sw, ...
        param, ...
        m_eq_in, ...
        begin_time, ...
        T_in),  ...
        tspan,    ...
        IC,       ...
        options);
    
    
    % store output, not the most proper way of coding but it works
    time  = [time;  X];
    P     = [P;     Y(:,1)];
    T     = [T;     Y(:,2)];
    eps_g = [eps_g; Y(:,3)];
    V     = [V;     Y(:,4)];
    rho_m = [rho_m; Y(:,5)];
    rho_x = [rho_x; Y(:,6)];
    
    eps_g(eps_g < 1e-8) = 0;
    
    % re-initialize
    P_0        = P(end);
    T_0        = T(end);
    eps_g0     = eps_g(end);
    V_0        = V(end);
    rho_m0     = rho_m(end);
    rho_x0     = rho_x(end);
    
    [eps_x0,~,~] =  crystal_fraction(T_0,eps_g0);
    rho_er_start = (1-eps_g0-eps_x0).*rho_m0 + eps_g0.*rho_g0 + eps_x0.*rho_x0;
    
    % re-initialize time
    begin_time = time(end);
    storeTime = storeTime(storeTime<time(end));
    storeTemp = storeTemp(storeTime<time(end));
    
    
    if isempty(IE)
        disp('you reached the end of time')
        phase0=5;
        
    elseif IE==4 && eps_x0<0.5
        
        
        % set up ode
        tspan       = [begin_time end_time]; % solve from t=0s to t=1e12s
        IC          = [P_0 T_0 eps_g0 V_0 rho_m0 rho_x0];   % store initial conditions
        sw.eruption = 1;
        options     = odeset('RelTol',rel_tol,'AbsTol',abs_tol,...
            'Events',@(z,y) stopChamber(y,sw.eruption, param.P_lit, P_crit),'Refine',2);
        
        [X,Y,TE,YE,IE] = ode15s(@(t,y) odeChamber(t,y,        ...
            sw, ...
            param, ...
            m_eq_in, ...
            begin_time, ...
            T_in),  ...
            tspan,    ...
            IC,       ...
            options);
        
        % store output
        time  = [time;  X];
        P     = [P;     Y(:,1)];
        T     = [T;     Y(:,2)];
        eps_g = [eps_g; Y(:,3)];
        V     = [V;     Y(:,4)];
        rho_m = [rho_m; Y(:,5)];
        rho_x = [rho_x; Y(:,6)];
        
        eps_g(eps_g < 1e-8) = 0;
        
        
        % store eruption times and durations
        dur = time(end)-begin_time;
        mass_er = dur*mdot_out;
        vol_er = mass_er/rho_er_start;
        
        time_erupt = [time_erupt begin_time];
        dur_erupt   = [dur_erupt dur];
        mass_erupt  = [mass_erupt mass_er];
        vol_erupt   = [vol_erupt vol_er];
        
        %re-initialize
        begin_time=time(end);
        P_0       = P(end);
        T_0       = T(end);
        eps_g0    = eps_g(end);
        V_0       = V(end);
        rho_m0    = rho_m(end);
        rho_x0    = rho_x(end);
        
        storeTime = storeTime(storeTime<time(end));
        storeTemp = storeTemp(storeTime<time(end));
        
    elseif IE==1
        disp('eps_g became 0.')
        phase0=1;
    elseif IE==2
        disp('eps_x became 0.')
        phase0=2;
    elseif IE==3
        disp('eps_x/(1-eps_g) became 0.8')
        phase0=3;
    elseif IE==6
        disp('eps_x became 0.5')
        phase0=4;
        %toc
    elseif IE==4 && eps_x0>=0.5
        disp('critical pressure reached but eps_x>0.5.')
        phase0=6;
    else
        disp('IE')
        disp(IE)
        disp('time')
        disp(time(end))
        disp('values')
        disp(eps_x0)
        disp(eps_g0)
        disp(P_0)
        disp(T_0)
        
        if ismember(1,IE)
            disp('eps_g became 0.')
            phase0=1;
        elseif ismember(2,IE)
            disp('eps_x became 0.')
            phase0=2;
        elseif ismember(3,IE)
            disp('eps_x/(1-eps_g) became 0.8')
            phase0=3;
        elseif ismember(6,IE)
            disp('eps_x became 0.5')
            phase0=4;
        else
            disp('keep calm and find out why the code stopped')
        end
    end
    
    ii=ii+1;
end



%crystal volume fraction
eps_x = zeros(size(P));
for i = 1:length(eps_x)
    [eps_x(i),~,~] =  crystal_fraction(T(i),eps_g(i));
end

% dissolved water mass fraction
m_eq = zeros(size(P));
for i = 1:length(m_eq)
    [m_eq(i),~,~] = exsolve(P(i),T(i));
end

% gas density
rho_g = zeros(size(P));
for i = 1:length(rho_g)
    [rho_g(i),~,~] = eos_g(P(i),T(i));
end

% bulk density
% global rho % MRT
rho            = (1-eps_g-eps_x).*rho_m + eps_g.*rho_g + eps_x.*rho_x;

% bulk heat capacity
c  = ((1-eps_g-eps_x).*rho_m.*param.c_m + eps_g.*rho_g.*param.c_g + eps_x.*rho_x.*param.c_x)./rho;
%
% % write file for output
% dlmwrite(['output/output' num2str(run_n,'%0.5i') '.dat'],...
%     [time P T eps_g V rho_m rho_x eps_x m_eq rho_g rho c],'precision',10)
%
% % write file for output
% dlmwrite(['output/stats' num2str(run_n,'%0.5i') '.dat'],...
%     [time_erupt' dur_erupt' mass_erupt' vol_erupt'],'precision',10)
%
% % write file for output
% dlmwrite(['output/eta' num2str(run_n,'%0.5i') '.dat'],...
%     [eta_r_vec' eta_t_vec'],'precision',10)


