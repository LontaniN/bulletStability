load('DATA/9x19Para_FMJ_124gr.mat')
LNOSE = geom.LNOSE;
DNOSE = geom.DNOSE;
LCENTR = geom.LCENTR;
DCENTR = geom.DCENTR;
LAFT = geom.LAFT;
DAFT = geom.DAFT;

%% GEOMETRY BUILDUP
L = LNOSE + LCENTR + LAFT;
x = linspace(0,L,1000);

% boat tail geometry function (ONLY CYLINDRICAL BOAT)
xboat = [LNOSE + LCENTR, L     ];
yboat = [DCENTR/2      , DAFT/2];
a_boat = (yboat(2) - yboat(1))/(xboat(2) - xboat(1));
b_boat = yboat(2) - (a_boat*xboat(2));

% radius of the entire rocket with axial coordinate
if LAFT == 0
    bulletRFUN = @(x) rOgive(x).*(x<LNOSE) + (DNOSE/2).*( x>=LNOSE & x<L ) ;
else
    bulletRFUN = @(x) rOgive(x).*(x<LNOSE) + (DNOSE/2).*( x>=LNOSE & x<(LNOSE+LCENTR) ) + (a_boat*x + b_boat).*(x>=(LNOSE+LCENTR) & x<L);
end

%% PLOT
Projectile2D = figure();
plot(-x,bulletRFUN(x),'LineWidth',1.5,'Color','black');
grid on
hold on
plot(-x,-bulletRFUN(x),'LineWidth',1.5,'Color','black');
plot(-XCG,0,'.','MarkerSize',25,'Color','red')
axis equal
xlabel('[m]')
ylabel('[m]')

q = 60;
xq = linspace(0,x(end),q);
rq = interp1(x,bulletRFUN(x),xq);
[X,Y,Z] = cylinder(rq,q);
Z = Z*x(end);

Projectile3D = figure();
s = [-45 170];
k = [0.6 .5 .6 5];
sl = surfl(Z,Y,X,s,k,'light');
sl(2).Color = 'w';
sl(1).EdgeColor = 'none';
sl(1).FaceColor = '#af7045';
%     'FaceAlpha',0.75,'EdgeAlpha',0,'FaceColor','black')
view(-135,25)
grid on
xlabel('[m]')
ylabel('[m]')
zlabel('[m]')
axis equal

%% TEST

% Projectile3D = figure();
% s = [-45 170];
% k = [0.5 .3 .4 6];
% sl = surfl(X,Y,-Z+L,s,k,'light');
% sl(2).Color = 'w';
% sl(1).EdgeColor = 'none';
% sl(1).FaceColor = '#af7045';
% %     'FaceAlpha',0.75,'EdgeAlpha',0,'FaceColor','black')
% view(-135,25)
% grid on
% xlabel('[m]')
% ylabel('[m]')
% zlabel('[m]')
% axis equal
% 
% hold on
% 
% load('DATA/9x19Para_FMJ_158gr.mat')
% sl = surfl(X+0.015,Y,-Z+L,s,k);
% sl(1).EdgeColor = 'none';
% sl(1).FaceColor = '#af7045';
% 
% load('DATA/9x39mm_SP5.mat')
% sl = surfl(X+0.03,Y,-Z+L,s,k);
% sl(1).EdgeColor = 'none';
% sl(1).FaceColor = '#af7045';
% 
% load('DATA/7mm_twenty_nine_hunt.mat')
% sl = surfl(X+0.045,Y,-Z+L,s,k);
% sl(1).EdgeColor = 'none';
% sl(1).FaceColor = '#af7045';
% 
% exportStandardizedFigure(gcf,'All_proj',0.85,'overwriteFigure',true,'addMarkers',false,'changeColors',false,'changeLineStyle',false)