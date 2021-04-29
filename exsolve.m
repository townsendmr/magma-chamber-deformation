function [meq,dmeqdP,dmeqdT] = exsolve(P,T) 

% % Henry's law
% satConst   = 4.11e-6;   % saturation (or Henry's) constant (Pa^-0.5)
% 
% meq        = satConst*P^0.5;
% dmeqdP     = 0.5*satConst*P^-0.5;
% dmeqdT     = 0;

% meq        = 0.05;%satConst*P^0.5;
% dmeqdP     = 0;%0.5*satConst*P^-0.5;
% dmeqdT     = 0;

% parametrization of Zhang and Behrens 2000 from Dufek and Bergantz 2005
meq    =   (P*1e-6)^0.5.*(0.4874 - 608/T + 489530/T^2) ...
         + (P*1e-6)*(-0.06062 + 135.6./T - 69200./T.^2)   ...
         + (P*1e-6)^1.5*(0.00253 - 4.154/T + 1509/T.^2);
dmeqdP = 0.5*(P*1e-6)^-0.5*(0.4874 - 608/T + 489530/T.^2) ...
          + (-0.06062 + 135.6/T - 69200/T.^2)   ...
          + 1.5*(P*1e-6)^0.5*(0.00253 - 4.154/T + 1509/T^2);
dmeqdT =   (P*1e-6)^0.5*( 608/T^2 - 2*489530/T^3) ...
          + (P*1e-6)*(-135.6/T^2 + 2*69200/T^3)   ...
          + (P*1e-6)^1.5*(4.154/T^2 - 2*1509/T^3);

meq     = 1e-2.*meq;
dmeqdP  = 1e-8.*dmeqdP;
dmeqdT  = 1e-2.*dmeqdT;