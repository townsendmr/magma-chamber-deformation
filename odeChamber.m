function dydz = odeChamber(time,y,        ...
    sw, ...
    param, ...
    m_eq_in, ...
    begin_time, ...
    T_in)


global storeTime storeTemp storeP storeV storeMF storeCONC

P              = y(1);
T              = y(2);
eps_g          = y(3);
V              = y(4);

% Set recharge rate
if strcmp(param.recharge, 'constant')
    mdot_in = param.rate;
    
elseif strcmp(param.recharge, 'pulse')
    
    if param.repeat == 1
        cycle = floor(time(end)/param.periodicity);
        peak_time2 = param.peak_time+param.periodicity*cycle;
        mdot_in = param.mass_injected*(1/(param.sig_time*sqrt(2*pi)))*exp(-0.5*((time(end)-peak_time2)./param.sig_time).^2);
        
    elseif param.repeat == 0
        mdot_in = param.mass_injected*(1/(param.sig_time*sqrt(2*pi)))*exp(-0.5*((time(end)-param.peak_time)./param.sig_time).^2);
        
    else
        disp('error: choose if recharge pulses are repeated')
    end
    
elseif strcmp(param.recharge, '2chamber')
    mdot_in = param.cond*(param.P_deep-P);
    
else
    disp('error: recharge conditions not defined')
end

if sw.plot_recharge == 1
plot(time(end)./(3600*24*365),mdot_in,'o')
xlabel('time')
ylabel('mass recharge rate (kg/s)')
hold on
x0=1e3;
y0=1e3;
width=550;
height=250;
set(gcf,'position',[x0,y0,width,height])
else
end

if ~isempty(storeTime) % if it's not the very first time step
    first_time_step = 0;
    if storeTime(end)==time % if the time step is repeated
        repeat = 1;
        storeTime(end) = time;
        storeTemp(end) =  T;
        storeP(end) = P;
        storeV(end) = V;
        dt = 0;
    else % if it's a new time step
        
        repeat = 0;
        storeTime = [storeTime; time];
        storeTemp = [storeTemp; T];
        storeP    = [storeP; P];
        storeV    = [storeV; V];
        dt = time-storeTime(end);
        
    end
else % if it is the first time step
    
    storeTime = [storeTime; time];
    storeTemp = [storeTemp; T];
    storeP    = [storeP; P];
    storeV    = [storeV; V];
    first_time_step = 1;
    dt = time;
end

cross=find(abs(diff(sign(diff(storeTime))))>0, 1); % number of time interutpions


if ~isempty(cross)
    cross_time=  storeTime(end);
    ind_cross = find(storeTime<cross_time);
    storeTime = [storeTime(ind_cross); cross_time];
    storeTemp = [storeTemp(ind_cross); storeTemp(end)];
    storeP    = [storeP(ind_cross); storeP(end)];
    storeV    = [storeV(ind_cross); storeV(end)];
end



dV_dP          = V/param.beta_r;
dV_dT          = -V*param.alpha_r;

rho_m          = y(5);
drho_m_dP      = rho_m/param.beta_m;
drho_m_dT      = -rho_m*param.alpha_m;

rho_x          = y(6);
drho_x_dP      = rho_x/param.beta_x;
drho_x_dT      = -rho_x*param.alpha_x;

[eps_x, deps_x_dT,deps_x_deps_g] = crystal_fraction(T,eps_g);

[rho_g,drho_g_dP,drho_g_dT]      = eos_g(P,T);

rho            = (1-eps_g-eps_x)*rho_m + eps_g*rho_g + eps_x*rho_x;
drho_dP        = (1-eps_g-eps_x)*drho_m_dP + eps_g*drho_g_dP + eps_x*drho_x_dP;
drho_dT        = (1-eps_g-eps_x)*drho_m_dT + eps_g*drho_g_dT + eps_x*drho_x_dT;
drho_deps_g    = -rho_m + rho_g;
drho_deps_x    = -rho_m + rho_x;

% exsolution
[m_eq,dm_eq_dP,dm_eq_dT] = exsolve(P,T);
c              = ((1-eps_g-eps_x)*rho_m*param.c_m + eps_g*rho_g*param.c_g + eps_x*rho_x*param.c_x)/rho;
dc_dP          = (1/rho)*((1-eps_g-eps_x)*param.c_m*drho_m_dP + eps_g*param.c_g*drho_g_dP + eps_x*param.c_x*drho_x_dP) - (c/rho)*drho_dP;
dc_dT          = (1/rho)*((1-eps_g-eps_x)*param.c_m*drho_m_dT + eps_g*param.c_g*drho_g_dT + eps_x*param.c_x*drho_x_dT) - (c/rho)*drho_dT;
dc_deps_g      = (1/rho)*(-rho_m*param.c_m  +  rho_g*param.c_g) - (c/rho)*drho_deps_g;
dc_deps_x      = (1/rho)*(-rho_m*param.c_m  +  rho_x*param.c_x) - (c/rho)*drho_deps_x;


% boundary conditions
[Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out, Hdot_in, Hdot_out, P_loss] ...
    = boundary_conditions(sw,param,P,T,eps_g,eps_x,V,rho_m,rho_x,m_eq,rho_g,rho,...
    c,m_eq_in, T_in, mdot_in,time);


% During repose
if sw.eruption == 0
    
    
    if eps_g > 1e-8
        phase = 3;
        
        mf_ex = eps_g*V*rho_g;
        mf_dis = m_eq*(1-eps_x-eps_g)*V*rho_m;
        mf = mf_ex + mf_dis;
        conc = m_eq;
        
    else
        
        if first_time_step == 1
            mf_last = param.mf_0;
        else
            mf_last = storeMF(end);
        end
        
        
        mf = mf_last + Mdot_in*m_eq_in*dt;
        conc = mf/((1-eps_x)*V*rho_m);
        
        if conc > m_eq
            eps_g = (conc-m_eq)*(1-eps_x)*rho_m/rho_g;
            phase = 3;
            
        else
            eps_g = 0;
            phase = 2;
            
        end
    end
end


% During an eruption
if sw.eruption == 1
    
    
    if eps_g > 1e-8
        
        phase = 3;
        
        mf_ex = eps_g*V*rho_g;
        mf_dis = m_eq*(1-eps_x-eps_g)*V*rho_m;
        mf = mf_ex + mf_dis;
        conc = m_eq;
    else
        
        if first_time_step == 1
            mf_last = param.mf_0;
            conc_last = param.conc_0;
        else
            mf_last = storeMF(end);
            conc_last = storeCONC(end);
        end
        
        % calculate total mass fluids at some time during eruption
        mf = mf_last + Mdot_in*m_eq_in*dt - Mdot_out*conc_last*dt;
        conc = mf/((1-eps_x)*V*rho_m);
        
        
        % if concentration is greater than m_eq, exsolve some gas and do 3
        % phase solution
        if conc > m_eq
            eps_g = (conc-m_eq)*(1-eps_x)*rho_m/rho_g;
            phase = 3;
            
        else % otherwise do 2 phase solution
            eps_g = 0;
            phase = 2;
            
        end
        
    end
end

% store mass of volatiles
if first_time_step == 0 % if it's not the very first time step
    if repeat == 1 % if the time step is repeated
        storeMF(end) = mf;
        storeCONC(end) = conc;
    else % if it's a new time step
        storeMF = [storeMF; mf];
        storeCONC = [storeCONC; conc];
    end
elseif first_time_step == 1  % if it is the first time step
    storeMF = [storeMF; mf];
    storeCONC = [storeCONC; conc];
end

if ~isempty(cross)
    storeMF    = [storeMF(ind_cross); storeMF(end)];
    storeCONC    = [storeCONC(ind_cross); storeCONC(end)];
end




% Build matrices

% coefficients in the system of unknowns Ax = B, here x= [dP/dt dT/dt deps_g/dt]
% note: P, T, and phi are y(1), y(2) and y(3) respectively
% values matrix A

% conservation of (total) mass
a11 = (1/rho)*drho_dP     + (1/V)*dV_dP;
a12 = (1/rho)*drho_dT     + (1/V)*dV_dT + (1/rho)*drho_deps_x*deps_x_dT;
a13 = (1/rho)*drho_deps_g               + (1/rho)*drho_deps_x*deps_x_deps_g;

% conservation of volatile mass
a21 = (1/rho_g)*drho_g_dP + (1/V)*dV_dP ...
    + (m_eq*rho_m*(1-eps_g-eps_x))/(rho_g*eps_g)*((1/m_eq)*dm_eq_dP + (1/rho_m)*drho_m_dP + (1/V)*dV_dP);
a22 = (1/rho_g)*drho_g_dT + (1/V)*dV_dT ...
    + (m_eq*rho_m*(1-eps_g-eps_x))/(rho_g*eps_g)*((1/m_eq)*dm_eq_dT + (1/rho_m)*drho_m_dT + (1/V)*dV_dT ...
    - deps_x_dT/(1-eps_g-eps_x));
a23 = 1/eps_g - (1+deps_x_deps_g)*m_eq*rho_m/(rho_g*eps_g);
% conservation of (total) enthalpy
a31 = (1/rho)*drho_dP      + (1/c)*dc_dP + (1/V)*dV_dP ...
    + (param.L_e*rho_g*eps_g)/(rho*c*T)*((1/rho_g)*drho_g_dP + (1/V)*dV_dP) ...
    - (param.L_m*rho_x*eps_x)/(rho*c*T)*((1/rho_x)*drho_x_dP + (1/V)*dV_dP);
a32 = (1/rho)*drho_dT      + (1/c)*dc_dT + (1/V)*dV_dT + 1/T ...
    + (param.L_e*rho_g*eps_g)/(rho*c*T)*((1/rho_g)*drho_g_dT  + (1/V)*dV_dT) ...
    - (param.L_m*rho_x*eps_x)/(rho*c*T)*((1/rho_x)*drho_x_dT + (1/V)*dV_dT)...
    + ((1/rho)*drho_deps_x + (1/c)*dc_deps_x - (param.L_m*rho_x)/(rho*c*T))*deps_x_dT;
a33 = (1/rho)*drho_deps_g  + (1/c)*dc_deps_g ...
    + (param.L_e*rho_g)/(rho*c*T) ...
    + ((1/rho)*drho_deps_x + (1/c)*dc_deps_x - (param.L_m*rho_x)/(rho*c*T))*deps_x_deps_g;

% values vector B
% conservation of (total) mass
b1  =  (Mdot_in - Mdot_out)/(rho*V) - P_loss;
% conservation of volatile mass
b2  =  (Mdot_v_in - Mdot_v_out)/(rho_g*eps_g*V) - P_loss*(1+(m_eq*rho_m*(1-eps_g-eps_x))/(rho_g*eps_g));
% conservation of (total) enthalpy
b3  =  (Hdot_in - Hdot_out)/(rho*c*T*V) - P_loss*(1-(param.L_m*rho_x*eps_x)/(rho*c*T)+(param.L_e*rho_g*eps_g)/(rho*c*T));


if phase == 2
    
    A          = [ a11 a12; a31 a32];
    A_P        = [ b1  a12;  b3  a32 ];
    A_T        = [ a11 b1 ; a31 b3 ];
    
    det_A          = det(A);
    dP_dt          = det(A_P)/det_A;
    dT_dt          = det(A_T)/det_A;
    deps_g_dt      = 0;
    
    dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss;
    drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt;
    drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt;
    
else
    
    A          = [ a11 a12 a13; a21 a22 a23; a31 a32 a33 ];
    A_P        = [ b1  a12 a13; b2  a22 a23; b3  a32 a33 ];
    A_T        = [ a11 b1  a13; a21 b2  a23; a31 b3  a33 ];
    A_eps_g    = [ a11 a12 b1 ; a21 a22 b2 ; a31 a32 b3  ];
    
    det_A          = det(A);
    dP_dt          = det(A_P)/det_A;
    dT_dt          = det(A_T)/det_A;
    deps_g_dt      = det(A_eps_g)/det_A;
    
    dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss;
    drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt;
    drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt;
    
end

dydz           = zeros(6,1);    % column vector
dydz(1)        = dP_dt;
dydz(2)        = dT_dt;
dydz(3)        = deps_g_dt;
dydz(4)        = dV_dt;
dydz(5)        = drho_m_dt;
dydz(6)        = drho_x_dt;

