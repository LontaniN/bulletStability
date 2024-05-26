
%% TRAJECTORY
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
    plot(sV(1:length(CD_star)),CD_star/adim)
    hold on
    plot(sV(1:length(CLa_star)),CLa_star/adim)
    ylabel('Coeffs [-]')
    xlabel('Downrange [cal]')
    grid on
end


L = 0.0284; % HARCODED FOR TEST
d = 0.0067056;
if flags.Animation

    amp = 1;
    margine = d/2;

    if flags.LowRes
        % LOW RESOLUTION ANIMATION
        figure()
        sl = surf(Z,Y,X);
        sl.FaceAlpha = 0.55;
        sl.FaceColor = "black";
        sl.EdgeAlpha = 1;
        obj = sl;
        axis equal
        xlabel('x')
        view(-129,17)

        startIndex = 1;
        speed = 2; % jump od indexes during for loop

        maxPitchDispl = max(deg2rad(alpha(startIndex:end)))*L + margine;
        maxYawDispl = max(deg2rad(beta(startIndex:end)))*L + margine;
        ylim([-maxYawDispl, maxYawDispl])
        zlim([-maxPitchDispl, maxPitchDispl])

        for j = startIndex:speed:length(alpha)


            if j == startIndex
                deltaAlpha = alpha(j);
                deltaBeta = beta(j);
            else
                deltaAlpha = alpha(j)-alpha(j-1);
                deltaBeta = beta(j)-beta(j-1);
            end

            rotate(obj,[0,1,0],-deltaAlpha*amp);
            rotate(obj,[0,0,1],-deltaBeta*amp);

            pause(0.0005)
            if flags.gifExport
                exportgraphics(gcf,'testAnimated.gif','Resolution',80,'Append',true);
            end
        end

    else
        % ANIMATION WITH LIGHTNING
        figure()
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

        startIndex = 1;
        speed = 2; % jump od indexes during for loop

        maxPitchDispl = max(deg2rad(alpha(startIndex:end)))*L + margine;
        maxYawDispl = max(deg2rad(beta(startIndex:end)))*L + margine;
        ylim([-maxYawDispl, maxYawDispl])
        zlim([-maxPitchDispl, maxPitchDispl])

        for j = startIndex:speed:length(alpha)


            if j == startIndex
                deltaAlpha = alpha(j);
                deltaBeta = beta(j);
            else
                deltaAlpha = alpha(j)-alpha(j-1);
                deltaBeta = beta(j)-beta(j-1);
            end

            rotate(obj,[0,1,0],-deltaAlpha*amp);
            rotate(obj,[0,0,1],-deltaBeta*amp);

            pause(0.0005)
            if flags.gifExport
                exportgraphics(gcf,'testAnimated.gif','Resolution',80,'Append',true);
            end
        end
    end
end

