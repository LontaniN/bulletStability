clc
clear
close all

i = sqrt(-1);

% flags.Stabilize = false; 
flags.RollDamp = true;     % activate the damping in roll, the roll rate will decrease as the bullet travels
flags.Animation = true;    % activate the 3D bullet motion animation
flags.LowRes = false;      % plot a lower resolution of the 3D bullet for the animation
flags.gifExport = false;   % export the 3D animation as a GIF
%% PARAMETERS
sMax = 15000;                % maximum downrange in calibers
sV = linspace(0,sMax,7000);  % vector of downrange
Ns = length(sV);

safeMargin = 1.1;             % safety margin to staying in the stable region

% barrel
twistRateInch = 1/12;                     % turn/inches
twistRate = twistRateInch * 2*pi/0.0254;  % rad/m

% Bullet / Projectile
load("DATA\9x19Para1.mat")

%%%%%%
CD = coeffs.CD;
CLa = coeffs.CLa;
CMa = coeffs.CMa;
CMqCMadot = coeffs.CMqCMadot;
Clp = coeffs.Clp;
CMpa = coeffs.CMpa;
%%%%%%

d = geom.DCENTR;
r = d/2;

% environment
Temp = 298;              % K
rho = 1.225;             % kg/m^2
g = 9.81;                % m/s^2
c = sqrt(1.4*287*Temp);  % m/s% speed of sound

%% INITIAL CONDITION
% v0 comes with the DATA, otherwise specify it yourself
v(1) = v0;
Ma0 = v0/c;         % initial Mach number

phi0 = 0 * pi/180;  % rad% elevation from local horizon

xi0 = 0 * pi/180;   % rad    % initial yaw at muzzle
betaDot = 25;       % rad/s  % initial yaw rate
xi0_prime = -i*d/v0 * betaDot;

%% SCALE COEFFS TO HAVE STAR COEFFS
adim = (rho*S*d)/(2*m);
CLa = adim*CLa;
CD = adim*CD;
CMa = adim*CMa;
CMpa= adim*CMpa;
Clp = adim*Clp;
CMqCMadot = adim*CMqCMadot;

%% COMPUTATIONS
p = twistRate * v0;   % rad/s% roll rate

RPM = p * 60/(2*pi);  % RPM rotation per minute

kx_2 = m*d^2 / Ix;
ky_2 = m*d^2 / Iy;

if length(states.MACH) == 1
    simpleAnalysis
else
    pseudoSimulation
end

%% PLOTS
PLOTS
