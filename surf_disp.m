a   = (3*V./(4*pi)).^(1/3);
a0  = a(1);  % initial chamber radius (m)
mu = param.mu_r;
pr = param.pr_r;
d = param.depth;

% Initial temperature and viscosity profile around chamber
maxn = 10000;
GLQ_n =64;
[quadpts,weights] = GLQ_points_weights_hard(GLQ_n);
cc    = 10*a0; %outer shell radius (m)
b     = a0 +cc;
quadpts_r = (b-a0)/2.*quadpts + (a0+b)/2;
storeTime =0;
storeTemp = T_0;
Trt = heat_conduction_chamber_profile(maxn,a0,cc,quadpts_r,param.kappa,param.Tb);
nn   = 1.9;   % power law exponent
AA   = 2e-4*(1e6)^-nn;  % material dependent constant (Pa^-n s^-1)
G    = 141e3; % activation energy for creep (J/mol)
M    = 8.314472;   % Molar gas constant (which gas?) (J/mol/K = kg m^2 s^-2 mol^-1 K^-1)
dev_stress = 20e6;
eta_rt     = (dev_stress.^(1-nn)./AA).*exp(G./M./Trt);

eta_rt = sort(eta_rt);
quadpts_r = sort(quadpts_r);
Trt = sort(Trt, 'descend');


if sw.visc_relax == 1
    
    if strcmp(param.viscModel, 'shell')
        
        
        % Determine radius of shell and viscosity of shell
        if param.shellOverride == 1
            R2  = param.R_shell;
            eta = param.eta_shell;
            
        else
            t_R =0; i = 2;
            while t_R<param.tObs
                R2 = quadpts_r(i);
                T_S = Trt(i);
                eta_shell   = calc_viscosity(a(1),R2, T_0, T_S);
                t_R = (3*eta_shell*(1-param.pr_r)*(R2/a(1))^3)/(param.mu_r*(1+param.pr_r));
                i = i+1;
            end
            
            eta = eta_shell;
            
        end
        
        tau = (3*eta*(1-pr)*(R2/a(1))^3)/(mu*(1+pr));
        
        
        % compute vertical displacments at surface directly overhead chamber
        uz_s=[];
        for j = 2:length(time)
            t = time(j);
            R1 = a(j);
            uz_s_t=0;ur_c_t=0;
            for i = 1:j-1
                t_since = t-time(i);
                P0      = P(i+1)-P(i);
                uz_s_t = uz_s_t+ ((1-pr)*P0*(R1^3)/(mu*d^2))*(exp(-t_since/tau)+((R2/R1)^3)*(1-exp(-t_since/tau)));
                
                ur_c_t = ur_c_t + ((P0*R1^3)/(4*mu))*( exp(-t_since/tau)/(R1^2) +...
                    ( (3*(1-pr)/(R1^2))*(R2/R1)^3 - 2*(1-2*pr)*R1/(R1^3))*((1-exp(-t_since/tau))/(1+pr)));
                
            end
            % cumulative vertical displacements
            uz_s(j)=uz_s_t; ur_c(j)=ur_c_t;
        end
        
        % figure
        % plot(time(1:length(uz_s))./(365*24*3600),uz_s,'b','LineWidth', 2)
        % title('Segall shell vertical surface displacements')
        %  xlabel('time in years')
        %  ylabel('cumulative displacements (cm)')
        %  set(gca,'FontSize',16)
        %
        
        
        % Plotting
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
        
        
    elseif strcmp(param.viscModel, 'full')
        
        
        b = 11*a(1);
        eta   = calc_viscosity(a(1), b, T_0, param.Tb);
        
        % coordiantes at surface overhead chamber
        r = 0; % radial distance from center
        z = 0; % depth
        R=sqrt(r^2 + (z-d)^2);
        Rs=sqrt(r^2 + (z+d)^2);
        
        K = param.beta_r;
        tau = eta/mu;
        alp = (3*K+mu)/(3*K);
        tau_a = alp*tau;
        A=1+2*mu/K;
        B=(2*(mu^2))/(K*(3*K+mu));
        uz_s=[]; ur_s=[]; R2_vec=[];ur_c=[];uz_s2=[];
        for j = 2:length(time)
            t = time(j);
            R1 = a(j);
            % coordinates at chamber wall
            r_c = R1;
            z_c = d;
            R_c=sqrt(r_c^2 + (z_c-d)^2);
            Rs_c=sqrt(r_c^2 + (z_c+d)^2);
            ur_c_t=0;uz_s_t=0;ur_s_t=0;uz_s_t2=0;
            for i = 1:j-1
                t_since = t-time(i);
                P0      = P(i+1)-P(i);
                
                % Equation 14 in Bonafede and Ferrari (2009)
                % Note that uplift = -uz (see B+F 09)
                uz_s_t = uz_s_t+ ((P0*R1^3)/(4*mu))*( ((z-d)/(R^3) +(2*z/(Rs^3)) -...
                    (6*z*(z+d)^2)/(Rs^5))*(1+t_since/tau) - ((z+d)/(Rs^3))*(A-B*exp(-t_since/tau_a) + t_since/tau) );
                % uz_s_t2 = uz_s_t2+((P0*R1^3)/(4*mu*d^2))*(-1-2*t_since/tau -A+B*exp(-t_since/tau_a));
                ur_c_t = ur_c_t+ ((P0*R1^3)/(4*mu))*( ((r_c)/(R_c^3) - (6*r_c*z_c*(z_c+d))/(Rs_c^5))*(1+t_since/tau) +...
                    ((r_c)/(Rs_c^3))*(A-B*exp(-t_since/tau_a) + t_since/tau) );
            end
            % cumulative displacements
            uz_s(j)=-uz_s_t;
            uz_s2(j)=-uz_s_t2;
            ur_c(j)=ur_c_t;
            
        end
    end
    
elseif sw.visc_relax == 0
    r = 0;
    uz_s_e=zeros(1,length(time)-1); ur_s_e=zeros(1,length(time)-1);
    R2_vec=[];    ur_c_t=0;uz_s_t=0;ur_s_t=0;
    for j = 2:length(time)
        R1 = a(j-1);
        
        P0      = P(j)-P(j-1);
        
        uz = ((1-pr)*P0.*(R1.^3)./mu).*(d./((r^2 + d^2).^(3/2)));
        ur = ((1-pr)*P0.*(R1.^3)./mu).*(r./((r^2 + d^2).^(3/2)));
        
        uz_s_t = uz_s_t + uz;
        ur_s_t = ur_s_t + ur;
        
        % cumulative elastic displacements
        uz_s_e(j)=uz_s_t;
        ur_s_e(j)=ur_s_t;
    end
    
    uz_s = uz_s_e;
    ur_s = ur_s_e;
    
end



% figure
% plot(time(1:length(ur_c))./(365*24*3600),ur_c,'LineWidth', 2)
% xlabel('time in years')
% ylabel('chamber wall displacement (m)')
% set(gca,'FontSize',16)


