function [value,isterminal,direction] = stopChamber(y,eruption, P_lit, P_crit)

P     = y(1);
T     = y(2);
eps_g = y(3);

[eps_x, ~,~]= crystal_fraction(T,eps_g);
value1a = 1;
value1b = eps_x;
value1c = eps_x./(1-eps_g)-0.8;%1-eps_x-eps_g-0.001;



if eruption ==0
    value2 = (P-P_lit)-P_crit;
elseif eruption==1
    value2 = P_lit-P;
else

    disp('eruption not specified')
end

Q_out          = 1;%heat_conduction_chamber;
if isnan(Q_out)
    value3=1;
%     disp('Q_out is NaN')
else
    value3=1;
end
% value2=1;


value4 = eps_x-0.5;


value      = [value1a; value1b; value1c; value2; value3;value4];
isterminal = [1; 1; 1; 1; 1; 1];   % Stop the integration
direction  = [0; 0; 0; 0; 1; 0];   % detect all zeros (default)