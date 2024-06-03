clc
clear
close all

i = sqrt(-1);

addpath("source\")
addpath("DATA\")

%% FLAGS
flags.RollDamp = true;       % activate the damping in roll, the roll rate will decrease as the bullet travels
flags.maxInitialYaw = true;  % compute the trajectory starting with the maximum yaw angle specified, otherwise use the minimum one computed through Kent's equation
flags.Animation = true;      % activate the 3D bullet motion animation
flags.LowRes = false;        % plot a lower resolution of the 3D bullet for the animation
flags.gifExport = false;      % export the 3D animation as a GIF

%% PARAMETERS
% Bullet / Projectile
load("DATA\9x19Para1.mat")
d = geom.DCENTR;
r = d/2;

% Downrange 
sMax_metres = 100;                    % m  % maximum downrange
sMax = sMax_metres/d;                 % maximum downrange in calibers
sV = linspace(0,sMax,12000);          % vector of downrange
Ns = length(sV);

% BARREL
twistRateInch = 1/9.84;                     % turn/inches
twistRate = twistRateInch * 2*pi/0.0254;    % rad/m

deltaMax = deg2rad(5); % rad  % MAXIMUM YAW AT MUZZLE, IT DETERMINES THE INITIAL COMPLEX YAW AND YAW RATE

% Animation parameters
amp = 1;                              % it multiplies the angles of the epycyclic trajectory, use it to better capture the movement if the angles are very small
startPoint = 0;                       % [%]  % percentage of the flight from which you want the animation to begin, example: 66% means that the animation will start at 2/3 of the computed flight
speed = 2;                            % speed up the animation, must be an integer!
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

initPhase = deg2rad(0); % rad  % initial phase of the yaw at muzzle

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
eps = sqrt(1 - 1/Sg) / (2*(Iy/Ix) - 1) * sin(deltaMax); % rad % IN-BORE YAW
xi0_prime = i * ((p*d)/v0) * sin(eps) * exp((i*initPhase)); % rad/m  % initial complex yaw-rate at muzzle

if flags.maxInitialYaw
    xi0 = sin(deltaMax) * exp((i*initPhase));   % rad  % initial complex yaw at muzzle
else
    xi0 = sin(eps) * exp((i*initPhase));   % rad  % initial complex yaw at muzzle  
end

betaDot = abs(v0/d * xi0_prime * i);

%% TRAJECTORY COMPUTATION
if length(states.MACH) == 1
    simpleAnalysis
else
    pseudoSimulation
end

%% PLOTS
PLOTS
