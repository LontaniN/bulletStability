
%% EPYCYCLIC TRAJECTORY
figure('Name','Tip trajectory')
plot3(sV(1:length(alpha)),beta,alpha,"LineWidth",1.1)
grid on
hold on
xlabel('downrange [cal]')
ylabel('yaw [deg]')
zlabel('pitch [deg]')
if exist('sTransonic','var')
    % transonic plane
    [y,z] = meshgrid(-20:40:20); % Generate x and y data
    x = sTransonic*ones(size(y, 1)); % Generate z data
    plane = surf(x, y, z); % Plot the surface
    plane.FaceColor = "black";
    plane.FaceAlpha = 0.1;
end
if exist('sSubsonic','var')
    % subsonic plane
    [y,z] = meshgrid(-20:40:20); % Generate x and y data
    x = sSubsonic*ones(size(y, 1)); % Generate z data
    plane = surf(x, y, z); % Plot the surface
    plane.FaceColor = "blue";
    plane.FaceAlpha = 0.1;
end
ylim([min(beta)*1.2,max(beta)*1.2])
zlim([min(alpha)*1.2,max(alpha)*1.2])


figure('Name','Tip trajectory')
plot(beta,alpha)
grid on
xlabel('yaw [deg]')
ylabel('pitch [deg]')
axis equal

%% YAW OF REPOSE
if length(Sd) > 1
    figure('Name','Yaw of repose')
    plot3(sV(1:length(alpha)),real(betaR),real(betaR),"LineWidth",1.1)
    grid on
    xlabel('downrange [cal]')
    ylabel('yaw component [deg]')
    zlabel('pitch component[deg]')
else
    fprintf('Yaw of repose = %0.3g deg\n', abs(betaR));
end
%% STABILITY DIAGRAM
SdVec = linspace(0,2,1000);
SgVec = 1./(SdVec.*(2-SdVec));

figure('Name','Stability factors')
plot(SdVec,SgVec,'LineWidth',1.5)
hold on
scatter(Sd,Sg,'filled')
% plot(Sd,Sg,'Marker','.','MarkerSize',20);
grid on
ylim([0,6])
xlim([0,2])
xlabel('Sd, dynamic stability factor')
ylabel('Sg, gyro stability factor')
text(1,4,'Stable Region','HorizontalAlignment','center')
text(0.4,0.5,'Slow mode unstable','HorizontalAlignment','center')
text(1.6,0.5,'Fast mode unstable','HorizontalAlignment','center')


if length(Sd) == 1
    annotation('textbox',...
        [0.2 0.76 0.15 0.11],...
        'String',{['Sg =' num2str(Sg,'%.2f')],['Sd =' num2str(Sd,'%.2f')]},...
        'FontSize',10,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',1,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0 0 0]);
end

if length(Sd) > 1
    figure('Name','Speed over downrange')
    plot(sV(1:length(v)),v)
    grid on
    xlabel('Downrange [cal]')
    ylabel('Speed [m/s]')

    figure('Name','Coeffs and Mach over downrange')
    yyaxis left
    plot(sV(1:length(Ma)),Ma)
    ylabel('Mach [-]')
    yyaxis right
    plot(sV(1:length(CD_V)),CD_V/adim)
    hold on
    plot(sV(1:length(CLa_V)),CLa_V/adim)
    ylabel('Coeffs [-]')
    xlabel('Downrange [cal]')
    grid on
end

%% ANIMATION
if flags.Animation

    if startPoint == 0
        startIndex = 1;
    else
        startIndex = startPoint/100 * length(alpha);
    end

    maxPitchDispl = max(amp*deg2rad(alpha(startIndex:end)))*L + margine;
    maxYawDispl = max(amp*deg2rad(beta(startIndex:end)))*L + margine;

    if flags.LowRes
        % LOW RESOLUTION ANIMATION
        f = figure();
        sl = surf(Z,Y,X);
        sl.FaceAlpha = 0.55;
        sl.FaceColor = "black";
        sl.EdgeAlpha = 1;
        obj = sl;
        axis equal
        xlabel('x')
        view(-129,17)

        ylim([-maxYawDispl, maxYawDispl])
        zlim([-maxPitchDispl, maxPitchDispl])

        for j = startIndex:round(speed):length(alpha)

            if j == startIndex
                deltaAlpha = alpha(j);
                deltaBeta = beta(j);
            else
                deltaAlpha = alpha(j) - alpha(j-round(speed));
                deltaBeta = beta(j) - beta(j-round(speed));
            end

            if Ma(j) >= 1.2
                set(f.Children.Title,'String','Supersonic','Color','Red');
            elseif Ma(j) < 1.2 && Ma(j) > 0.85
                set(f.Children.Title,'String','Transonic','Color',"#D95319");
            else
                set(f.Children.Title,'String','Subsonic','Color','black');
            end

            rotate(obj,[0,1,0],deltaAlpha*amp);
            rotate(obj,[0,0,1],deltaBeta*amp);

            pause(0.0005)
            if flags.gifExport
                if exist(animationFileName,'file') && j == startIndex
                    error("The GIF file already exists, change the file name in MAIN.m")
                end
                if floor(i/3) == i/3
                    exportgraphics(f,animationFileName,'Append',true,'Resolution',120);
                end
            end
        end

    else
        % ANIMATION WITH LIGHTNING
        f = figure();
        s = [-45 170];
        k = [0.6 .5 .6 5];
        sl = surfl(Z,Y,X,s,k,'light');
        sl(2).Color = 'w';
        sl(1).EdgeColor = 'none';
        sl(1).FaceColor = '#af7045';
        obj = sl(1);
        axis equal
        xlabel('x')
        view(-129,17)
        % round(2/3*length(alpha)):2

        ylim([-maxYawDispl, maxYawDispl])
        zlim([-maxPitchDispl, maxPitchDispl])
        
        i = 0;
        for j = startIndex:round(speed):length(alpha)
            i = i + 1;
            if j == startIndex
                deltaAlpha = alpha(j);
                deltaBeta = beta(j);
            else
                deltaAlpha = alpha(j) - alpha(j-round(speed));
                deltaBeta = beta(j) - beta(j-round(speed));
            end

            if Ma(j) >= 1.2
                set(f.Children.Title,'String','Supersonic','Color','Red');
            elseif Ma(j) < 1.2 && Ma(j) > 0.85
                set(f.Children.Title,'String','Transonic','Color',"#D95319");
            else
                set(f.Children.Title,'String','Subsonic','Color','black');
            end

            rotate(obj,[0,1,0],deltaAlpha*amp);
            rotate(obj,[0,0,1],deltaBeta*amp);

            pause(0.0005)
            if flags.gifExport
                if exist(animationFileName,'file') && j == startIndex
                    error("The GIF file already exists, change the file name in MAIN.m")
                end
                if floor(i/3) == i/3
                    exportgraphics(f,animationFileName,'Append',true,'Resolution',120);
                end
            end
        end
    end
end

