%% COMPUTE TRAJECTORY OF MODES
Ma = Ma0 * ones(1,length(sV));  % define constant Mach array for Animation
[alpha,beta,betaR,lambdaF,lambdaS] = trajectory(P,M,T,G,H,xi0,xi0_prime,sV);
alphaTot = sqrt(alpha.^2 + beta.^2);
DR = i.*1/ky_2 * (P*G0/2) * (CLa/CMa) * sV.^2 .* (1 + 2/3*(CD.*sV) + 1/3*(CD.*sV).^2);
%% DYNAMIC STABILITY
Sd = 2*T / H;
SgLimit = 1/(Sd*(2-Sd));
%% DATA PRINT
fprintf('Sg = %0.2f\n',Sg)
fprintf('Sd = %0.2f\n',Sd)
fprintf('Roll rate = %0.0f RPM\n',RPM);
fprintf('Barrel twist rate = 1 x %0.2f turn x inches\n', 1/twistRateInch);
