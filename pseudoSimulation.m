p = ones(Ns,1) * p;
deltaS = 0;
GyroCriteria = inf;
for j = 1:Ns
    s = sV(j);
    if j>2
        if alpha(j-1)<25 && beta(j-1)<25
            v(j) = v0 * exp(-trapz(sV(2:j),CD_V));
            deltaS = sV(j) - sV(j-1);
            if flags.RollDamp
                Kp = -(kx_2 * Clp + CD_V(j-1));
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
        sTransonic = sV(j);
    end
    if Ma(j) < 0.85 && Ma(j-1) > 0.85
        sSubsonic = sV(j);
    end

    % Interpolate coefficients
    CLa_V(j) = interp1(states.MACH, CLa, Ma(j));
    CD_V(j) = interp1(states.MACH, CD, Ma(j));
    CMa_V(j) = interp1(states.MACH, CMa, Ma(j));
    CMqCMadot_V(j) = interp1(states.MACH, CMqCMadot, Ma(j));

    P(j) = (Ix/Iy)*((p(j)*d)/v(j));
    M(j) = ky_2 * CMa_V(j);
    T(j) = CLa_V(j) + kx_2 * CMpa;
    G(j) = g * d * cos(phi0)/v(j)^2;
    H(j) = CLa_V(j) - CD_V(j) - ky_2*CMqCMadot_V(j);

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