function [P_loss] = viscous_P_loss_shell(quadpts_r,Trt,param)



global storeTime storeSumk storeSumk_2 storeSumk_old storeSumk_2_old lengthTime maxTime
global storeSumk_oldold storeSumk_2_oldold switch_Tprofile storeP storeV

% time grid
time = storeTime;
time_index= length(time);

% if maxTime < time(end)
%     maxTime = time(end);
%     storeSumk_oldold=storeSumk_old;
%     storeSumk_2_oldold=storeSumk_2_old;
%     storeSumk_old=storeSumk;
%     storeSumk_2_old=storeSumk_2;
% end


r = ((3/(4*pi))*storeV).^(1/3); % chamber radii at each timestep


mu = param.mu_r;
pr = param.pr_r;

% Determine radius of shell and viscosity of shell
if param.shellOverride == 1
    
    R2  = param.R_shell;
    eta = param.eta_shell;
    tau = (3*eta*(1-pr)*(R2/r(1))^3)/(mu*(1+pr));

else
    
    t_R =0; i = 2;
    while t_R<param.tObs
        R2 = quadpts_r(i);
        T_S = Trt(i);
        eta_shell   = calc_viscosity(r(end),R2, Trt(1), T_S);
        t_R = (3*eta_shell*(1-param.pr_r)*(R2/r(end))^3)/(param.mu_r*(1+param.pr_r));
        i = i+1;
    end
    eta = eta_shell;
    tau = (3*eta*(1-pr)*(R2/r(end))^3)/(mu*(1+pr));

end



% WHAT I SHOULD DO HERE IS DO ALL THE STEP EXCEPT THE LAST 2 SO THAT IF
% ADAPTIVE TIME STEP I AM COVERED
if time_index == 2
    
    
    R1 = r(end);
    
    dP = storeP(2)-storeP(1);
    dt = time(2)-time(1);
    P_loss = ((dP*R1*exp(-dt/tau))/(4*mu*(1+pr)*tau))*(1-5*pr+3*(1-pr)*(R2/R1)^3);
    
    
elseif time_index > 2 && time_index>lengthTime
    P_loss = 0;
    for k = 1:time_index-1
        
        dP = storeP(k+1)-storeP(k);
        dt = time(time_index)-time(k);
        R1 = r(k);
        
        P_loss_k = ((dP*R1*exp(-dt/tau))/(4*mu*(1+pr)*tau))*(1-5*pr+3*(1-pr)*(R2/R1)^3);
        P_loss = P_loss + P_loss_k;
        
    end
    
elseif time_index > 2 && time_index==lengthTime
    
    P_loss = 0;
    for k = 1:time_index-1
        
        dP = storeP(k+1)-storeP(k);
        dt = time(time_index)-time(k);
        R1 = r(k);
        
        P_loss_k = ((dP*R1*exp(-dt/tau))/(4*mu*(1+pr)*tau))*(1-5*pr+3*(1-pr)*(R2/R1)^3);
        P_loss = P_loss + P_loss_k;
        
    end
    
    
elseif time_index > 2 && time_index<lengthTime
    
    P_loss = 0;
    for k = 1:time_index-2
        
        dP = storeP(k+1)-storeP(k);
        dt = time(time_index)-time(k);
        R1 = r(k);
        
        P_loss_k = ((dP*R1*exp(-dt/tau))/(4*mu*(1+pr)*tau))*(1-5*pr+3*(1-pr)*(R2/R1)^3);
        P_loss = P_loss + P_loss_k;
        
    end
    
elseif time_index <2
    P_loss = 0;
end

P_loss = (3/r(end))*P_loss;

lengthTime=time_index;
