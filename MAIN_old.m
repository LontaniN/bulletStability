clc
clear
close all

i = sqrt(-1);

flags.Stabilize = true;   % compute the roll rate needed to stabilize the bullet and compute the new motion
flags.Animation = false;  % activate the 3D bullet motion animation
flags.gifExport = false;  % export the 3D animation as a GIF

%% PARAMETERS
sV = linspace(0,20000,10000);             % downrange vector
Ns = length(sV);

safeMargin = 1.1;                         % safety margin to staying in the stable region

% barrel
twistRateInch = 1/12;                     % turn/inches
twistRate = twistRateInch * 2*pi/0.0254;  % rad/m

load("DATA\SierraDataEx.mat")


if exist("m",'var') == 0
    % geometry
    d = 7.82 * 1e-3;         % m
    r = d/2;                 % m
    L = 3.98 * d;            % m% reference length
    S = pi * r^2;            % m^2% reference area


    % mass
    m = 10.88 * 1e-3;        % kg
    Ix = 7.2131e-08;         % hardcoded for test
    Iy = 5.3743e-07;

    if exist("Ix")==0
        Ix = 1/2 * m * r^2;  % kg*m^2% approximation of solid cylinder
        Iy = 1/12 * m * (3*r^2 + L^2);
    end

    % aerodynamics
    CLa = 2.75;              % 1/rad
    CD = 0.331;
    CMa = 2.63;
    CMpa = -0.22;
    CMqCMadot = -8.2;
end

% environment
Temp = 298;              % K
rho = 1.225;             % kg/m^2
g = 9.81;                % m/s^2
c = sqrt(1.4*287*Temp);  % m/s% speed of sound

%% INITIAL CONDITION
v0 = 792.48;        % m/s % muzzle velocity
Ma = v0/c * ones(Ns,1);

phi0 = 0 * pi/180;  % rad % elevation from local horizon

xi0 = 0 * pi/180;   % rad % initial yaw at muzzle
betaDot = 25;       % rad/s % initial yaw rate
xi0_prime = -i*d/v0 * betaDot;
%% SCALE COEFFS TO HAVE STAR COEFFS
adim = (rho*S*d)/(2*m);
CLa_star = adim*CLa; CD_star = adim*CD; CMa_star = adim*CMa; CMpa_star = adim*CMpa; CMqCMadot_star = adim*CMqCMadot;

%% COMPUTATIONS
p = twistRate * v0;   % rad/s% roll rate
RPM = p * 60/(2*pi);  % RPM rotation per minute

kx_2 = m*d^2 / Ix;
ky_2 = m*d^2 / Iy;

P = (Ix/Iy)*((p*d)/v0);
M = ky_2 * CMa_star;
T = CLa_star + kx_2 * CMpa_star;
G = g*d*cos(phi0)/v0^2;
H = CLa_star - CD_star - ky_2*CMqCMadot_star;

%% COMPUTE TRAJECTORY OF MODES
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

%% PLOT
PLOTS
if flags.Stabilize
    %% CORRECTION FOR INSTABILITY
    if Sg<SgLimit
        SgLimit = SgLimit * safeMargin; % New Gyroscopic Stability Margin with a safety margin applied
        PLimit = sqrt(SgLimit*4*M) ;
        pLimit = PLimit * (Iy/Ix) * v0/d;
        RPMLimit = pLimit * 60/(2*pi);
        twistRateLimit = pLimit/v0; % rad/m
        twistRateLimit = twistRateLimit * 0.0254/(2*pi); % turn x inches

        fprintf('\n')
        fprintf('|| CORRECTION ||\n');
        fprintf('Sg needed = %0.2f\n',SgLimit)
        fprintf('Values needed for dynamic stability:\n');
        fprintf('Roll rate = %0.0f RPM\n',RPMLimit);
        fprintf('Barrel twist rate = 1 x %0.2f turn x inches\n', 1/twistRateLimit);

        %% COMPUTE TRAJECTORY OF MODES
        sVL = linspace(0,50000,10000); % downrange vector
        [alphaL,betaL,betaRL,lambdaFL,lambdaSL] = trajectory(PLimit,M,T,G,H,xi0,xi0_prime,sVL);

        %% PLOT
        figure('Name','Tip trajectory')
        plot3(sVL,betaL,alphaL)
        grid on
        xlabel('downrange [cal]')
        ylabel('yaw [deg]')
        zlabel('pitch [deg]')


        figure('Name','Tip trajectory')
        plot(betaL,alphaL)
        grid on
        xlabel('yaw [deg]')
        ylabel('pitch [deg]')
        axis equal


        SdVec = linspace(0,2,1000);
        SgVec = 1./(SdVec.*(2-SdVec));

        figure('Name','Stability factors')
        plot(SdVec,SgVec,'LineWidth',1.5)
        hold on
        plot(Sd,SgLimit,'Marker','.','MarkerSize',20);
        grid on
        ylim([0,6])
        xlabel('Sd, dynamic stability factor')
        ylabel('Sg, gyro stability factor')
        text(1,4,'Stable Region','HorizontalAlignment','center')
        text(0.4,0.5,'Slow mode unstable','HorizontalAlignment','center')
        text(1.6,0.5,'Fast mode unstable','HorizontalAlignment','center')
    end
end