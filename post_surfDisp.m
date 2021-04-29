
%run_n =1;
input      = dlmread(['input/input' num2str(run_n,'%0.5i') '.dat']);
output     = dlmread(['output/output' num2str(run_n,'%0.5i') '.dat']);


time_erupt = [];
dur_erupt  = [];
mass_erupt = [];
vol_erupt  = [];

q = ['output/stats' num2str(run_n,'%0.5i') '.dat'];

r = fopen(q);
tline = fgetl(r);
if tline == -1
    time_erupt = [];
    dur_erupt  = [];
    mass_erupt = [];
    vol_erupt  = [];
else
    stats      = dlmread(['output/stats' num2str(run_n,'%0.5i') '.dat']);
    % get stats
    time_erupt = stats(:,1)';
    dur_erupt  = stats(:,2)';
    mass_erupt = stats(:,3)';
    vol_erupt  = stats(:,4)';
end
fclose(r);

%eta         = dlmread(['output/eta' num2str(run_n,'%0.5i') '.dat']);

% get input
T_0               = input(2);
a0                = input(4);  % initial chamber radius (m)
mdot_in           = input(7);  % mass inflow rate (kg/s)
tau_inj           = input(10);
tau_cool          = input(11);
tau_visco_elastic = input(12);
P_lit             = input(13); % lithostatic pressure (Pa)
Tb                = input(14); % temperature at edge of shell
P_crit            = input(15); % critical overpressure (Pa)
rho_r             = input(16); % density of crust (kg/m^3)
mdot_out          = input(17); % mass outflow during eruptions (kg/s)
wt_g              = input(18); % water content in resident magma
m_eq_in           = input(19); % water content of inflowing magma
d                 = input(21); % depth (m)
beta_r            = 1e10;      % bulk modulus country rock (Pa)
kappa             = 1e-6; % thermal diffusivity
pr                = 0.25; % poisson's ratio
mu                = 3*beta_r*(1-2*pr)/(2*(1+pr)); % convert elastic moduli to get shear mod (Pa)


% get output
time     = output(:,1);
P        = output(:,2);
T        = output(:,3);
eps_g    = output(:,4);
V        = output(:,5);
 a        = (3*V./(4*pi)).^(1/3);
rho_m    = output(:,6);
rho_x    = output(:,7);
eps_x    = output(:,8);
m_eq     = output(:,9);
rho_g    = output(:,10);
rho      = output(:,11);
c        = output(:,12);


% First find an effective radius and viscosity for the shell based on initial
% temperature/viscosity profile

% viscous relaxation stuff
nn         = 1.9;   % power law exponent
AA         = 2e-4*(1e6)^-nn;  % material dependent constant (Pa^-n s^-1)
G          = 141e3; % activation energy for creep (J/mol)
M          = 8.314472;   % Molar gas constant (which gas?) (J/mol/K = kg m^2 s^-2 mol^-1 K^-1)
GLQ_n      = 64;
[quadpts,weights] = GLQ_points_weights_hard(GLQ_n);
cc         = 10*a0;
b          = a0 +cc;
quadpts_r  = (b-a0)/2.*quadpts + (a0+b)/2;
Trt        = (a0.*T_0.*(b-quadpts_r) + b.*Tb.*(quadpts_r-a0))./(quadpts_r.*(b-a0));
dev_stress = 20e6;%P_0;%200e6;
eta_rt     = (dev_stress.^(1-nn)./AA).*exp(G./M./Trt);
I          = (b-a0)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
eta_r0     = 3*a0^3*I;

%eta = eta_r0; % use this as the effective viscosity
% ind_cold = find(log10(eta_rt)>20); % find location beyond which too cold
%R2 = min(quadpts_r(ind_cold)); % radius of viscoelastic shell
% 
%   figure
%   plot(quadpts_r,log10(eta_rt),'o')

%% Solution from Segall 2010 Ch 7


mu = param.mu_r;
pr = param.pr_r;

if sw.visc_relax == 1
% Determine radius of shell and viscosity of shell
if param.shellOverride == 1
    R2  = param.R_shell;
    eta = param.eta_shell;
    
%     ind_shell=find(quadpts_r<R2);
%     eta_log_mean = mean(log10(eta_rt(ind_shell)));
% eta = 10^eta_log_mean
else
ind_shell = find(log10(eta_rt)<param.eta_cutoff);
eta_log_mean = mean(log10(eta_rt(ind_shell)));
eta = 10^eta_log_mean
ind_R2 = find(eta_rt==max(eta_rt(ind_shell)));
R2 = quadpts_r(ind_R2)
end

tau = (3*eta*(1-pr)*(R2/a(1))^3)/(mu*(1+pr));

% coordiantes at surface overhead chamber

uz_s=[]; 
for j = 2:length(time)
    t = time(j);
    R1 = a(j);      
    uz_s_t=0;
    for i = 1:j-1
        t_since = t-time(i);
        P0      = P(i+1)-P(i);                  
        uz_s_t = uz_s_t+ ((1-pr)*P0*(R1^3)/(mu*d^2))*(exp(-t_since/tau)+((R2/R1)^3)*(1-exp(-t_since/tau)));                
    end
    % cumulative displacements
    uz_s(j)=uz_s_t; 
end

% figure
% plot(time(1:length(uz_s))./(365*24*3600),uz_s.*100,'b','LineWidth', 2)
% title('Segall shell vertical surface displacements')
%  xlabel('time in years')
%  ylabel('cumulative displacements (cm)')
%  set(gca,'FontSize',16)



%% Plotting
% figure
% plot(time(1:length(uz_s))./(365*24*3600),uz_s.*100,'LineWidth', 2)
% hold on
% plot(time(1:length(uz_s))./(365*24*3600),uz_s.*100,'LineWidth', 2)
% 
% xlabel('time in years')
% ylabel('cumulative displacements (cm)')
% legend('vertical','radial')
% title('position at center')
% set(gca,'FontSize',16)
% 
% Compute elastic displacements due to Mogi source for comparison

elseif sw.visc_relax == 0

uz_s_e=zeros(1,length(time)-1); ur_s_e=zeros(1,length(time)-1);
R2_vec=[];    ur_c_t=0;uz_s_t=0;ur_s_t=0;
for j = 2:length(time)
    R1 = a(j-1);
    
    P0      = P(j)-P(j-1);
    
    uz = ((1-pr)*P0.*(R1.^3)./mu).*(d./((r^2 + d^2).^(3/2)));
    ur = ((1-pr)*P0.*(R1.^3)./mu).*(r./((r^2 + d^2).^(3/2)));
    
    uz_s_t = uz_s_t + uz;
    ur_s_t = ur_s_t + ur;
    
    % cumulative displacements
    uz_s_e(j)=uz_s_t;
    ur_s_e(j)=ur_s_t;
end

uz_s = uz_s_e;
ur_s = ur_s_e;

end


figure
plot(time(1:length(uz_s))./(365*24*3600),uz_s.*100,'LineWidth', 2)
hold on
%plot(time(1:end)./(365*24*3600),uz_s_e.*100,'LineWidth', 2)
%legend(['R2 = ' R2])
xlabel('time in years')
ylabel('cumulative displacements (cm)')
title('comparing models')
set(gca,'FontSize',16)
%axis([0 end_time/(3600*24*365) min(uz_s.*100) max(uz_s.*100)])
%ylim([0 18])
plot(t_d,Data001(:,2),'ro')

% % Scale the displacements like in Le Mevel et al 2015
% %Viscoelastic solution
% % no gas
% ind_t0 = 11;% 10, 1 spike
% ind_t0 = 7;% 25,10 middle spike
% ind_t0 = 8;% 50,25 middle spike
% % with gas
% ind_t0 = 13;% 10, 1 spike
% %ind_t0 = 1;% 25,10 middle spike
% %ind_t0 = 1;% 50,25 middle spike
% 
% t0 = time(ind_t0);
% du = uz_s-uz_s(ind_t0);
% du_max = du(46);
% %max(du);
% ind_max=46;
% %find(du==du_max);
% t_max = time(ind_max);
% % Elastic solution
% du_e = uz_s_e-uz_s_e(ind_t0);
% du_max_e = du_e(46);
% %max(du_e);
% ind_max_e=46;
% %find(du_e==du_max_e);
% t_max_e = time(ind_max_e);
% 
% figure
% plot((time-t0)./(t_max-t0),du./du_max,'LineWidth', 2)
% hold on
% plot((time-t0)./(t_max_e-t0),du_e./du_max_e,'LineWidth', 2)
% legend('Viscoelastic','Elastic')
% xlabel('(t-t0)/(t_max-t0)')
% ylabel('uz/uz_max')
% title('normalized displacements')
% set(gca,'FontSize',16)
% axis([0 1 0 1]);



% % Displacements at radial coordinate = 0.1*depth
% %
% % r = 0.1*d; % radial distance from center, to evaluate surface displacements
% %
% % uz_s=zeros(length(time)-1); ur_s=zeros(length(time)-1);
% % uz_s=[]; ur_s=[]; R2_vec=[];
% % for j = 2:length(time)
% % t = time(j);
% % R1 = a(j);
% %     ur_c_t=0;uz_s_t=0;ur_s_t=0;
% % for i = 1:j-1
% %     t_since = t-time(i);
% %     P0      = P(i+1)-P(i);
% %     ur_c_t = ur_c_t+disp_chamber_wall(P0,t_since,R1,R2,eta,pr,mu);
% %     [uz, ur] = disp_surf(P0,t_since,R1,R2,eta,pr,mu,d,r);
% %     uz_s_t = uz_s_t + uz;
% %     ur_s_t = ur_s_t + ur;
% %     radius_p = a0+ur_c_t;
% % end
% % % cumulative displacements
% % uz_s(j)=uz_s_t;
% % ur_s(j)=ur_s_t;
% % end
% 
% % figure
% % plot(time(1:length(uz_s))./(365*24*3600),uz_s.*100, 'LineWidth', 2)
% % hold on
% % plot(time(1:length(uz_s))./(365*24*3600),ur_s.*100,'LineWidth', 2)
% % xlabel('time in years')
% % ylabel('cumulative displacements (cm)')
% % legend('vertical','radial')
% % title('radial position = 10% chamber depth')
% 
% 
% 
% 
