clc
clear
close all

i = sqrt(-1);

flags.Stabilize = false;  % compute the roll rate needed to stabilize the bullet and compute the new motion
flags.RollDamp = true;    % activate the damping in roll, the roll rate will decrease as the bullet travels
flags.Animation = true;   % activate the 3D bullet motion animation
flags.gifExport = false;  % export the 3D animation as a GIF
%% PARAMETERS
sMax = 120000;                            % maximum downrange in calibers
sV = linspace(0,sMax,15000);              % vector of downrange
Ns = length(sV);

safeMargin = 1.1;                         % safety margin to staying in the stable region

% barrel
twistRateInch = 1/12;                     % turn/inches
twistRate = twistRateInch * 2*pi/0.0254;  % rad/m

load('SierraDataVar.mat')

if exist("m",'var') == 0
    % geometry
    d = 7.82 * 1e-3;                      % m
    r = d/2;                              % m
    L = 3.98 * d;                         % m% reference length
    S = pi * r^2;                         % m^2% reference area


    % mass
    m = 10.88 * 1e-3;                     % kg
    Ix = 7.2131e-08;                      % hardcoded for test
    Iy = 5.3743e-07;

    if exist("Ix")==0
        Ix = 1/2 * m * r^2;               % kg*m^2% approximation of solid cylinder
        Iy = 1/12 * m * (3*r^2 + L^2);
    end
end

% environment
Temp = 298;                               % K
rho = 1.225;                              % kg/m^2
g = 9.81;                                 % m/s^2
c = sqrt(1.4*287*Temp);                   % m/s% speed of sound

%% INITIAL CONDITION
v0 = 792;           % m/s% muzzle velocity
v(1) = v0;
Ma0 = v0/c;         % initial Mach number

phi0 = 0 * pi/180;  % rad% elevation from local horizon

xi0 = 0 * pi/180;   % rad% initial yaw at muzzle
betaDot = 25;       % rad/s% initial yaw rate
xi0_prime = -i*d/v0 * betaDot;
%% SCALE COEFFS TO HAVE STAR COEFFS

adim = (rho*S*d)/(2*m);
CLa(2,:) = adim*CLa(2,:);
CD(2,:) = adim*CD(2,:);
CMa(2,:) = adim*CMa(2,:);
CMpa_star = adim*CMpa;
Clp_star = adim*Clp;
CMqCMadot(2,:) = adim*CMqCMadot(2,:);

%% COMPUTATIONS
p = twistRate * v0;     % rad/s% roll rate
p = ones(Ns,1) * p;

RPM = p * 60/(2*pi);    % RPM rotation per minute

kx_2 = m*d^2 / Ix;
ky_2 = m*d^2 / Iy;

deltaS = 0;
GyroCriteria = inf;
for j = 1:Ns
    s = sV(j);
    if j>2
        if alpha(j-1)<25 && beta(j-1)<25  % Checks if the angles stays within reasonable limit for the model
            v(j) = v0 * exp(-trapz(sV(2:j),CD_star));
            deltaS = sV(j) - sV(j-1);
            if flags.RollDamp
                Kp = -(kx_2 * Clp_star + CD_star(j-1));
                pdv_0 = p(1)*d/v0;
                p(j) = v(j)/d * pdv_0 * exp(-Kp * sV(j));
            end
        else
            break;
        end

    else
        v(j) = v0;
    end

    Ma(j) = v(j)/c;
    if Ma(j) < 1.2 && Ma(j-1) > 1.2
        sTransonic = sV(j);               % saves the downrange for which the bullet enters the transonic regime
    end
    if Ma(j) < 0.85 && Ma(j-1) > 0.85
        sSubsonic = sV(j);                % saves the downrange for which the bullet enters the subsonic regime
    end

    % Interpolate coefficients
    CLa_star(j) = interp1(CLa(1,:), CLa(2,:), Ma(j));
    CD_star(j) = interp1(CD(1,:), CD(2,:), Ma(j));
    CMa_star = interp1(CMa(1,:), CMa(2,:), Ma(j));
    CMqCMadot_star(j) = interp1(CMqCMadot(1,:), CMqCMadot(2,:), Ma(j));

    P(j) = (Ix/Iy)*((p(j)*d)/v(j));
    M(j) = ky_2 * CMa_star;
    T(j) = CLa_star(j) + kx_2 * CMpa_star;
    G(j) = g * d * cos(phi0)/v(j)^2;
    H(j) = CLa_star(j) - CD_star(j) - ky_2*CMqCMadot_star(j);

    GyroCriteria(j) = P(j)^2 - 4*M(j);

    % compute stability factor
    Sg(j) = P(j)^2 / (4*M(j));
    Sd(j) = 2*T(j) / H(j);
    SgLimit(j) = 1/(Sd(j)*(2-Sd(j)));


    % compute trajectory of modes
    [alpha(j),beta(j),betaR(j),lambdaF(j),lambdaS(j)] = trajectoryPseudoSim(P(j),M(j),T(j),G(j),H(j),xi0,xi0_prime,s);

end

%% DATA PRINT
fprintf('Roll rate = %0.0f RPM\n',RPM(1));
fprintf('Barrel twist rate = 1 x %0.2f turn x inches\n', 1/twistRateInch);

%% PLOT
PLOTS

if flags.Stabilize
    %% CORRECTION FOR INSTABILITY
    if Sg(1)<SgLimit(1)  % the minimum Sg occurs at the muzzle
        PLimit = sqrt(SgLimit(1)*4*M(1));
        pLimit = PLimit * (Iy/Ix) * v0/d;
        RPMLimit = pLimit * 60/(2*pi);
        twistRateLimit = pLimit/v0; % rad/m
        twistRateLimit = twistRateLimit * 0.0254/(2*pi); % turn x inches

        fprintf('\n')
        fprintf('|| CORRECTION ||\n');
        fprintf('Sg needed = %0.2f at muzzle exit\n',SgLimit(1))
        fprintf('Values needed for dynamic stability:\n');
        fprintf('Roll rate = %0.0f RPM\n',RPMLimit);
        fprintf('Barrel twist rate = 1 x %0.2f turn x inches\n', 1/twistRateLimit);

        %% COMPUTE TRAJECTORY OF MODES
        deltaS = 0;
        for j = 1:Ns
            s = sV(j);
            if j>2
                v(j) = v0 * exp(-trapz(sV(2:j),CD_star_corr));
                deltaS = sV(j) - sV(j-1);
            else
                v(j) = v0;
            end
            Ma(j) = v(j)/c;

            % Interpolate coefficients
            CLa_star_corr = interp1(CLa(1,:), CLa(2,:), Ma(j));
            CD_star_corr(j) = interp1(CD(1,:), CD(2,:), Ma(j));
            CMa_star_corr = interp1(CMa(1,:), CMa(2,:), Ma(j));
            CMqCMadot_star_corr = interp1(CMqCMadot(1,:), CMqCMadot(2,:), Ma(j));

            P(j) = (Ix/Iy)*((safeMargin*pLimit*d)/v(j));
            M(j) = ky_2 * CMa_star_corr;
            T(j) = CLa_star_corr + kx_2 * CMpa_star;
            G(j) = g*d*cos(phi0)/v(j)^2;
            H(j) = CLa_star_corr - CD_star_corr(j) - ky_2*CMqCMadot_star_corr;

            % compute stability factor
            Sg(j) = P(j)^2 / (4*M(j));
            Sd(j) = 2*T(j) / H(j);
            SgLimit(j) = 1/(Sd(j)*(2-Sd(j)));

            % compute trajectory of modes
            [alpha(j),beta(j),betaR(j),lambdaF(j),lambdaS(j)] = trajectoryPseudoSim(P(j),M(j),T(j),G(j),H(j),xi0,xi0_prime,s);

        end
        %% PLOT
        PLOTS
    end
end