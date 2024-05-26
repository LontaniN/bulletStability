clc
clear
close all

%% OGIVE
LNOSE = 0.014618;
DNOSE = 0.0067056;

% TNOSE = ;
nPOWER = 1/2;

xOgive = linspace(0,LNOSE,300);
rOgive = DNOSE/2 .* (xOgive./LNOSE).^nPOWER;

%% BODY
LCENTR = 0.009723;
DCENTR = DNOSE;

xBody = [0,LCENTR] + xOgive(end);
rBody = DCENTR/2 * ones(1,2);

%% BOATTAIL
LAFT = 0.00402;
DAFT = 0.005646;

xAft = [0,LAFT] + xBody(end);
rAft = [DCENTR/2,DAFT/2];


%% PLOT
x = [xOgive, xBody(2), xAft(2), xAft(end)+1e-5];
r = [rOgive, rBody(2), rAft(2), 0];

figure()
plot(x,r,'LineWidth',1.5,'Color','black');
grid on
hold on
plot(x,-r,'LineWidth',1.5,'Color','black');
axis equal

q = 60;
xq = linspace(0,x(end),q);
rq = interp1(x,r,xq);
[X,Y,Z] = cylinder(rq,q);
Z = Z*x(end);

figure()
s = [-45 170];
k = [0.6 .5 .6 5];
sl = surfl(Z,Y,X,s,k,'light');
sl(2).Color = 'w';
sl(1).EdgeColor = 'none';
sl(1).FaceColor = '#af7045';
%     'FaceAlpha',0.75,'EdgeAlpha',0,'FaceColor','black')
view(-135,25)
grid on
axis equal



