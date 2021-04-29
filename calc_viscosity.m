function eta_r = calc_viscosity(a, b, T_0, Tb)


nn         = 1.9;   % power law exponent
AA         = 2e-4*(1e6)^-nn;  % material dependent constant (Pa^-n s^-1)
G          = 141e3; % activation energy for creep (J/mol)
M          = 8.314472;   % Molar gas constant (which gas?) (J/mol/K = kg m^2 s^-2 mol^-1 K^-1)
GLQ_n      = 64;
[quadpts,weights] = GLQ_points_weights_hard(GLQ_n);
quadpts_r  = (b-a)/2.*quadpts + (a+b)/2;
Trt        = (a.*T_0.*(b-quadpts_r) + b.*Tb.*(quadpts_r-a))./(quadpts_r.*(b-a));

A = 4.25e7; % material-dependent constant for viscosity law (Pa s)
B = 8.31; % molar gas constant (J/mol/K)
G = 141e3; % activation energy for creep (J/mol)

%eta_rt     = (dev_stress.^(1-nn)./AA).*exp(G./M./Trt);
eta_rt     = A.*exp(G./B./Trt);
I          = (b-a)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
eta_r     = 3*a^3*I;
