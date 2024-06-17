clc
clear
close all

i = sqrt(-1);

addpath("source\")
addpath("DATA\")

%% FLAGS
flags.RollDamp = true;       % activate the damping in roll, the roll rate will decrease as the bullet travels
flags.Animation = false;      % activate the 3D bullet motion animation
flags.LowRes = false;        % plot a lower resolution of the 3D bullet for the animation
flags.gifExport = false;      % export the 3D animation as a GIF

%% PARAMETERS
% Bullet / Projectile
load("DATA\9x19Para.mat")
d = geom.DCENTR;
r = d/2;

% Downrange 
sMax_metres = 50;                    % m  % maximum downrange
sMax = sMax_metres/d;                 % maximum downrange in calibers
sV = linspace(0,sMax,12000);          % vector of downrange
Ns = length(sV);

% BARREL
twistRateInch = 1/9.84;                           % turn/inches
twistRate = twistRateInch * 2*pi/0.0254;          % rad/m
twistRateCalTurn = 1/(twistRateInch * d/0.0254);  % cal/turn
n = twistRateCalTurn;

deltaMax = deg2rad(5); % rad  % MAXIMUM FIRST YAW, USUALLY MEASURED THROUGH EXPERIMENTS

% Animation parameters
amp = 1;                              % it multiplies the angles of the epycyclic trajectory, use it to better capture the movement if the angles are very small
startPoint = 0;                       % [%]  % percentage of the flight from which you want the animation to begin, example: 66% means that the animation will start at 2/3 of the computed flight
speed = 4;                            % speed up the animation, must be an integer!
margine = d/2;                        % adjust the borders of the 3D plot of the animation
animationFileName = "Animation.gif";  % name of the gif file that will be generate if flags.gifExport is set to true

%%%%%%
CD = coeffs.CD;
CLa = coeffs.CLa;
CMa = coeffs.CMa;
CMqCMadot = coeffs.CMqCMadot;
Clp = coeffs.Clp;
CMpa = coeffs.CMpa;
%%%%%%

%% INITIAL CONDITIONS
% environment
Temp = 298;              % K
rho = 1.225;             % kg/m^2
g = 9.81;                % m/s^2
c = sqrt(1.4*287*Temp);  % m/s% speed of sound

% v0 comes with the DATA, otherwise specify it yourself
v(1) = v0;
Ma0 = v0/c;         % initial Mach number

phi0 = 0 * pi/180;  % rad    % elevation from local horizon

initPhase = deg2rad(0); % rad  % initial phase of the yaw at muzzle, angle between the vertically up-ward plane and the plane containing both the in-bore yaw and the bore axis

%% SCALE COEFFS TO HAVE STAR COEFFS (From theory)
adim = (rho*S*d)/(2*m);
CLa = adim*CLa;
CD = adim*CD;
CMa = adim*CMa;
CMpa= adim*CMpa;
Clp = adim*Clp;
CMqCMadot = adim*CMqCMadot;

%% IMPORTANT QUANTITIES
p = twistRate * v0;   % rad/s% roll rate
RPM = p * 60/(2*pi);  % RPM rotation per minute

kx_2 = m*d^2 / Ix;
ky_2 = m*d^2 / Iy;

%% AERODYNAMIC COEFFICIENTS AT MUZZLE
if length(states.MACH) > 1
    CMa_muzzle = interp1(states.MACH, CMa, Ma0);
    CLa_muzzle = interp1(states.MACH, CLa, Ma0);
    CD_muzzle = interp1(states.MACH, CD, Ma0);
    CMqCMadot_muzzle = interp1(states.MACH, CMqCMadot, Ma0);
else
    CMa_muzzle = CMa;
    CLa_muzzle = CLa;
    CD_muzzle = CD;
    CMqCMadot_muzzle = CMqCMadot;
end

P = (Ix/Iy)*((p*d)/v0);
M = ky_2 * CMa_muzzle;
T = CLa_muzzle + kx_2 * CMpa;
G = g*d*cos(phi0)/v0^2;
H = CLa_muzzle - CD_muzzle - ky_2*CMqCMadot_muzzle;

%% COMPUTE STABILITY FACTORS AT MUZZLE
Sg = P^2 / (4*M);
Sd = 2*T / H;
SgLimit = 1/(Sd*(2-Sd));

%% YAW AT MUZZLE, INITIAL PERTURBATION
eps = sqrt(1 - 1/Sg) * sin(deltaMax)/(2*(Iy/Ix) - 1); % rad % IN-BORE YAW
xi0 = sin(eps) * exp((i*initPhase));   % rad  % initial complex yaw at muzzle  
xi0_prime = i * ((p*d)/v0) * sin(eps) * exp((i*initPhase)); % rad/m  % initial complex yaw-rate at muzzle

betaDot = abs(v0/d * xi0_prime * i);  % rad/s  % initial yaw-rate

%% TRAJECTORY COMPUTATION
if length(states.MACH) == 1
    simpleAnalysis
else
    pseudoSimulation
end

%% MILLER'S FORMULA
L_inch = 1/0.0254 * L; % inches
d_inch = 1/0.0254 * d; % inches
m_grain = 15432 * m; % grains
l = L/d; % cal

t = sqrt(30*m_grain/(2 * d_inch^3 * l * (1+l^2)));
T = t * d_inch; % twist rate in inches x turn

fprintf('Advised twist rate from Miller formula: 1 x %0.2f [turn x inches]\n',T)

%% FLAT-FIRE DEFLECTION
% in modern bullets external ballistics, the main two quantities that contribute 
% to the deflection of the bullet from the expected trajectory are
% the aerodynamic jump JA and the lateral throw-off TL
% which depends on the in-bore yaw eps
CLa0 = CLa_muzzle/adim;
CMa0 = CMa_muzzle/adim;
JA = -i * ((2*pi/n) * (1/ky_2 - 1/kx_2) * (CLa0/CMa0) * sin(eps)) * exp(i*initPhase);

LN = geom.LNOSE;
LCYL = geom.LCENTR;
TL = i * (2*pi/n * (LN + LCYL/2 - XCG) * sin(eps)) * exp(i*initPhase);

Deflection = norm(JA+TL) * sMax_metres * 1e2; %cm
fprintf('The deflection at %d is %0.2f cm\n',sMax_metres,Deflection)

nMin_inches = 7;
nMax_inches = 16;
nMin = 1/(1/nMin_inches * d/0.0254);
nMax = 1/(1/nMax_inches * d/0.0254);
nV = linspace(nMin,nMax,100);

JA_v = NaN(1,length(nV));
TL_v = NaN(1,length(nV));
Def_v = NaN(1,length(nV));

for j = 1:length(nV)
    JA_v(j) = i * ((2*pi/nV(j)) * (1/ky_2 - 1/kx_2) * (CLa0/CMa0) * sin(eps)) * exp(i*initPhase) * sMax_metres * 1e2;
    TL_v(j) = -i * (2*pi/nV(j) * (LN + LCYL/2 - XCG) * sin(eps)) * exp(i*initPhase) * sMax_metres * 1e2;
    Def_v(j) = norm(JA_v(j) - TL_v(j)); 

    TL_v(j) = sign( dot( [real(JA_v(j)) , imag(JA_v(j))] , [real(TL_v(j)) , imag(TL_v(j))] ) ) * norm(TL_v(j));
    JA_v(j) = norm(JA_v(j)); % always considered positive for standard projectiles

end
nV_inches = linspace(nMin_inches,nMax_inches,100);
%% PLOTS
PLOTS
