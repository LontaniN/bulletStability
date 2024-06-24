function [alpha,beta,betaR,lambdaF,lambdaS,fastFreq,slowFreq] = trajectory(P,M,T,G,H,xi0,xi0_prime,s,V,d)

i = sqrt(-1);
%% COMPUTE THE TRAJECTORY
% phase derivatives wrt downrange
phiF_prime = 1/2 * (P + sqrt(P^2 - 4*M));
phiS_prime = 1/2 * (P - sqrt(P^2 - 4*M));

% exponential of the magnitude of the modes
lambdaF = -1/2 * (H - (P*(2*T - H))/(sqrt(P^2 - 4*M)));
lambdaS = -1/2 * (H + (P*(2*T - H))/(sqrt(P^2 - 4*M)));

% initial magnitudes of the modes

KF0_exp = -(( i*xi0_prime + phiS_prime*xi0 )/( phiF_prime - phiS_prime )); %* 1/exp(i*phiF0);
KS0_exp = (( i*xi0_prime + phiF_prime*xi0 )/( phiF_prime - phiS_prime )); %* 1/exp(i*phiS0);

% yaw of repose
betaR = (P*G)/(M + i*P*T);

% magnitudes of modes
fastTerm = KF0_exp * exp((lambdaF + i*phiF_prime).*s);
slowTerm = KS0_exp * exp((lambdaS + i*phiS_prime).*s);

fastFreq = (phiF_prime*V)/(2*pi*d);
slowFreq = (phiS_prime*V)/(2*pi*d);
% Compute the points of the trajectory
xi = fastTerm + slowTerm + i*betaR;
alpha = real(xi) * 180/pi;
beta = imag(xi) * 180/pi;

betaR = betaR * 180/pi;
end