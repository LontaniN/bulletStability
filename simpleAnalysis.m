
P = (Ix/Iy)*((p*d)/v0);
M = ky_2 * CMa;
T = CLa + kx_2 * CMpa;
G = g*d*cos(phi0)/v0^2;
H = CLa - CD - ky_2*CMqCMadot;

%% COMPUTE TRAJECTORY OF MODES
Ma = Ma0 * ones(1,length(sV));  % define constant Mach array for Animation
[alpha,beta,betaR,lambdaF,lambdaS] = trajectory(P,M,T,G,H,xi0,xi0_prime,sV);

%% COMPUTE STABILITY FACTORS
Sg = P^2 / (4*M);
Sd = 2*T / H;
SgLimit = 1/(Sd*(2-Sd));

%% DATA PRINT
fprintf('Sg = %0.2f\n',Sg)
fprintf('Sd = %0.2f\n',Sd)
fprintf('Roll rate = %0.0f RPM\n',RPM);
fprintf('Barrel twist rate = 1 x %0.2f turn x inches\n', 1/twistRateInch);
