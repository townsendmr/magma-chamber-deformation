function [P_loss] = viscous_P_loss_full(eta_r,param)

% Starting with Bonafede and Ferrari (2009) model for viscoelastic
% deformation due to spherical magma chamber in viscoelastic half space
% Here, I am only calculating the viscous component. The elastic component
% has been subtracted out and is accounted for elsewhere.

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


r = ((3/(4*pi))*storeV).^(1/3); % chamber radius 

K = param.beta_r;
mu = param.mu_r;
d = param.depth;
B = (2*mu^2)/(K*(3*K+mu));


% WHAT I SHOULD DO HERE IS DO ALL THE STEP EXCEPT THE LAST 2 SO THAT IF
% ADAPTIVE TIME STEP I AM COVERED
if time_index == 2
    
    dP = storeP(2)-storeP(1);
    dt = time(2)-time(1);
    a = r(1);
    alp = (3*K+mu)/(3*K);    
    tau = eta_r/mu;
    tau_a = tau*alp;
    P_loss = ((dP*a^3)/(4*mu))*((1/tau)*((a^-2)+(a*(a^2-8*d^2))/((a^2+4*d^2)^(5/2)))+...
        ((1/tau_a)*a*B/((a^2+4*d^2)^(3/2)))*exp(-dt/tau));
   
          
elseif time_index > 2 && time_index>lengthTime 
 P_loss = 0;
    for k = 1:time_index-1
               
    dP = storeP(k+1)-storeP(k);
    dt = time(time_index)-time(k);
    a = r(k);
    alp = (3*K+mu)/(3*K);   
    tau = eta_r/mu;
    tau_a = tau*alp;
    P_loss_k = ((dP*a^3)/(4*mu))*((1/tau)*((a^-2)+(a*(a^2-8*d^2))/((a^2+4*d^2)^(5/2)))+...
        ((1/tau_a)*a*B/((a^2+4*d^2)^(3/2)))*exp(-dt/tau));
    P_loss = P_loss + P_loss_k;

    end   
        
elseif time_index > 2 && time_index==lengthTime 

 P_loss = 0;
    for k = 1:time_index-1
               
    dP = storeP(k+1)-storeP(k);
    dt = time(time_index)-time(k);
    a = r(k);
    alp = (3*K+mu)/(3*K);   
    tau = eta_r/mu;
    tau_a = tau*alp;
    P_loss_k = ((dP*a^3)/(4*mu))*((1/tau)*((a^-2)+(a*(a^2-8*d^2))/((a^2+4*d^2)^(5/2)))+...
        ((1/tau_a)*a*B/((a^2+4*d^2)^(3/2)))*exp(-dt/tau));
    P_loss = P_loss + P_loss_k;

    end      
    

elseif time_index > 2 && time_index<lengthTime

   P_loss = 0;
    for k = 1:time_index-2
               
    dP = storeP(k+1)-storeP(k);
    dt = time(time_index)-time(k);
    a = r(k);
    alp = (3*K+mu)/(3*K);   
    tau = eta_r/mu;
    tau_a = tau*alp;
    P_loss_k = ((dP*a^3)/(4*mu))*((1/tau)*((a^-2)+(a*(a^2-8*d^2))/((a^2+4*d^2)^(5/2)))+...
        ((1/tau_a)*a*B/((a^2+4*d^2)^(3/2)))*exp(-dt/tau));
    P_loss = P_loss + P_loss_k;

    end      

elseif time_index <2
    P_loss = 0;
end

P_loss = (3/r(end))*P_loss;

lengthTime=time_index;
