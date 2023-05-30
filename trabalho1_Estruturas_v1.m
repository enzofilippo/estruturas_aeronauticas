clc
clear all
close all

syms x real;

% - Propriedades atmosf?ricas - %
rho = 1.225;
g = 9.81;

% - Propriedades do avi?ozinho - %
m = 5500;
W = m*g;
V = 100;
c = 2;
b = 12;
S = b*c;
n = 7.5;
alpha = deg2rad(0)
lambda = 1;

% - Cargas externas - %
% Sustenta??o %
L = n*W;

% Polar de arrasto %
Cd0 = 0.2;
k = 0.3;
k1 = 0.4;
CL = 2*L/(rho*V^2*S);
CD = Cd0 + k1*CL + k*CL^2;
D = (1/2)*rho*V^2*CD*S;

% Momentos %
mx = 0;
my = 0;
mz = 0;

% - Distribui??o de sustenta??o e arrasto - %

Ldist = plotDistribuicaoSchrenk(L, b, lambda,'Sustenta??o [kN/m]');
Ddist = plotDistribuicaoSchrenk(D, b, lambda,'Arrasto [kN/m]');

% - C?lculo das cargas externas nos eixos da aeronave - %

fy = Ddist*cos(alpha)- Ldist*sin(alpha);
vpa(fy,4)
fz = Ddist*sin(alpha)+ Ldist*cos(alpha);
vpa(fz,4)
plotCargasExternas(b, fy, fz)

% - C?lculo dos esfor?os internos nos eixos da aeronave - %
Vy=distV(fy,b);
Vz=distV(fz,b);

My=distM(my,Vz,1,b);
Mz=distM(mz,Vy,-1,b);

plotEsforcosInternos(Vy,Vz,b,'Esfor?o cortante','V')
plotEsforcosInternos(-fy,-fz,b,'for?as','F')
plotEsforcosInternos(My,Mz,b,'Momento fletor','M')


% - Propriedades do material 1 - %
E1 = 70e9;
v1 = 0.3;
sigmae1 = 400e6;

% - Propriedades do material 2 - %
E2 = 110e9;
v2 = 0.3;
sigmae2 = 1100e6;

% - Propriedades da se??o transversal - %
h = 0.2;
t1 = 0.05;
t2 = 0.08;
t3 = 0.04;
coordinateMatrix = plotSecaoTransversal(c,h,t1,t2,t3);

% - Propriedades de ?rea ponderadas da se??o transversal e posi??o do centroide de ?rea ponderado - %
E0 = E1;
[A_pond,centroide_pond,I_pond] = PropriedadesPonderadas(coordinateMatrix,E0,E1,E2);

% - Eixo neutro - %
angEixoNeutro = EixoNeutro(My,Mz,I_pond); %Em ?

% - Distribui??o de Tens?o - %
sigmaX1 = distTensao(My,Mz,I_pond,E1,E0);
sigmaX2 = distTensao(My,Mz,I_pond,E2,E0);

% - Distribui??o de deflex?es - %
w = def(My,E0,I_pond(:,1),-1);
v = def(Mz,E0,I_pond(:,2),1);

plotDeflexao(v,b,'V','y')
plotDeflexao(w,b,'W','z')
plotDeformada(v,b,'Vista Superior','y')
plotDeformada(w,b,'Vista Lateral','z')

figure()
hold on
plotTensoes(sigmaX1,coordinateMatrix(1,:),"Horizontal",'b','r')
plotTensoes(sigmaX1,coordinateMatrix(2,:),"Horizontal",'b','r')
plotTensoes(sigmaX2,coordinateMatrix(3,:),"Vertical",'k','m')
plotTensoes(sigmaX2,coordinateMatrix(4,:),"Vertical",'k','m')
plotTensoes(sigmaX2,coordinateMatrix(5,:),"Vertical",'k','m')
plotTensoes(sigmaX2,coordinateMatrix(6,:),"Vertical",'k','m')

function [y0,z0] = getMidLine(coords,type)
    %coords = [y z w h]
    Npontos = 30
    if (type == "Horizontal")
        y0 = linspace(coords(1),coords(3),Npontos)
        z0 = coords(2) + coords(4)/2
    else
        y0 = coords(1) + coords(3)/2
        z0 = linspace(coords(2),coords(4),Npontos)
    end
end

function sigmaBar = evalSigma(sigma, y0, z0, type)
    syms x y z;
    Npontos = 30;
    sigmaBar = zeros(Npontos,1);
    sigmaBarRaiz = subs(sigma,x,0);
    if type == "Horizontal"
        sigmaBarRaizT = subs(sigmaBarRaiz,z,z0);
        for i = 1:length(y0)
            sigmaBar(i) = eval(subs(sigmaBarRaizT,y,y0(i)));
        end
    else
        sigmaBarRaizT = subs(sigmaBarRaiz,y,y0);
        for i = 1:length(z0)
           sigmaBar(i) = eval(subs(sigmaBarRaizT,z,z0(i)));
        end
    end
    
end

function plotTensoes(sigmaX,coords,type,colorBar,colorArrow)
    syms x y z;
    [y0, z0] = getMidLine(coords,type)
    sigmaBar = evalSigma(sigmaX, y0, z0, type)
    len = size(sigmaBar,1)
    zeroVec = zeros(len,1)
    micro = 1e-6
    hold on    
    if type == "Horizontal"
        plot3([0 0], [y0(1) y0(end)], [z0 z0], colorBar, 'LineWidth', 2); % Plot da barra
        hold on
        quiver3(zeroVec,y0,zeroVec+z0,sigmaBar*micro,zeroVec,zeroVec, colorArrow,'MaxHeadSize',0.05,'AutoScale', 'off','AutoScaleFactor',0.1,'LineWidth',0.5);
    else
        plot3([0 0], [y0 y0], [z0(1) z0(end)], colorBar, 'LineWidth', 2); % Plot da barra
        hold on
        quiver3(zeroVec,zeroVec+y0,z0',sigmaBar*micro,zeroVec,zeroVec, colorArrow,'MaxHeadSize',0.05,'AutoScale', 'off','AutoScaleFactor',0.1,'LineWidth',0.5);
    end
        % Set the axis labels and title
    xlabel('x [MPa]');
    ylabel('y [m]');
    zlabel('z [m]');
    grid on

    % Adjust the figure view for better visualization
    view(30, 30);
end

function plotCargasExternas(b,fy,fz)
Npontos = 50;
x = linspace(0, b/2, Npontos); % Adjust the number of points as needed
mili = 0.001
% Evaluate the function handle at each point
evaluatedY = subs(fy,x)
loadsY = eval(evaluatedY)*mili
evaluatedZ = subs(fz,x)
loadsZ = eval(evaluatedZ)*mili

% Create the figure and plot the beam as a line
figure;
plot3([0 b/2], [0 0], [0 0], 'k-', 'LineWidth', 2); % Beam line
hold on
len = size(loadsY,2)
zeroVec = zeros(1,len)
quiver3(x,zeroVec,zeroVec,zeroVec,loadsY,zeroVec, 'r','MaxHeadSize',0.05,'AutoScale', 'off','AutoScaleFactor',0.1,'LineWidth',0.5);
hold on
quiver3(x,zeroVec,zeroVec,zeroVec,zeroVec,loadsZ, 'b','MaxHeadSize',0.05,'AutoScale', 'off','AutoScaleFactor',0.1,'LineWidth',0.5);
grid on

% Set the axis labels and title
xlabel('x [m]');
ylabel('y [kN/m]');
zlabel('z [kN/m]');

% Adjust the figure view for better visualization
view(30, 30);
end

function plotEsforcosInternos(Ey,Ez,b,titulo,label)
syms x real;
mili = 0.001;
figure()
subplot(1,2,1)
hold on
fplot(Ey*mili,[0,b/2],'LineWidth',1.2)
hold off
grid
caption1 = sprintf('%s em y', titulo);
title(caption1);
xlabel('x [m]')
label1 = sprintf('%s_y [kN]', label);
ylabel(label1)

subplot(1,2,2)
hold on
fplot(Ez*mili,[0,b/2],'LineWidth',1.2)
hold off
grid
caption2= sprintf('%s em z', titulo);
title(caption2);
xlabel('x [m]')
label2 = sprintf('%s_z [kN]', label);
ylabel(label2)
end

function coordinateMatrix = plotSecaoTransversal(c,h,t1,t2,t3)
%[xPos yPos width height]
lowerBar = [0 0 c+t2 t1];
upperBar = [0 h c+t2 t1];
leftBar = [0 t1 t2 h-t1];
midBar1 = [t2/2+c/3-t3/2 t1 t3 h-t1];
midBar2 = [t2/2+2*c/3-t3/2 t1 t3 h-t1];
rightBar = [c t1 t2 h-t1];
coordinateMatrix = [lowerBar; upperBar; rightBar; leftBar; midBar1; midBar2];
colorVec = ['b'; 'b'; 'g'; 'g'; 'g'; 'g'];
figure()
hold on
for i = 1:size(coordinateMatrix,1)
    rectangle('Position',coordinateMatrix(i,:),'FaceColor', colorVec(i))
end

limSupX = (c+t2);
marginX = 0.1*limSupX;
xlim([-marginX,limSupX+marginX])
limSupY = (h+t1);
marginY = 0.1*limSupY;
ylim([-marginY,limSupY+marginY])
ylabel('Altura - z [m]')
xlabel('Profundidade - y [m]')
end

function [A,A_pond] = areaPond(EPondVec, coordinateMatrix)
    A_pond = 0;
    A = zeros(size(coordinateMatrix,1),1);
    for i = 1:size(coordinateMatrix,1)
        A(i) = coordinateMatrix(i,3)*coordinateMatrix(i,4);
        A_pond = A_pond + coordinateMatrix(i,3)*coordinateMatrix(i,4)*EPondVec(i);
    end
end

function [centroide_pond, Y0Z0] = calcCentroide(EPondVec, coordinateMatrix,A)
    Y0 = zeros(size(coordinateMatrix,1),1);
    Z0 = zeros(size(coordinateMatrix,1),1);
    Y0_pond = 0; Z0_pond = 0;
    for i = 1:size(coordinateMatrix,1)
        Y0(i) = coordinateMatrix(i,1) + coordinateMatrix(i,3)/2;
        Y0_pond = Y0_pond +  Y0(i)*EPondVec(i)*A(i);
        Z0(i) = coordinateMatrix(i,2) + coordinateMatrix(i,4)/2;
        Z0_pond = Z0_pond + Y0(i)*EPondVec(i)*A(i);        
    end
    centroide_pond = [Y0_pond Z0_pond];
    Y0Z0 = [Y0 Z0];
end

function IyyIzz = calcMomInerciaLocal(coordinateMatrix)
    Iyy = zeros(size(coordinateMatrix,1),1);
    Izz = zeros(size(coordinateMatrix,1),1);
    for i = 1:size(coordinateMatrix,1)
        Iyy(i) = (coordinateMatrix(i,3)*coordinateMatrix(i,4)^3)/12;
        Izz(i) = (coordinateMatrix(i,4)*coordinateMatrix(i,3)^3)/12;
    end
    IyyIzz = [Iyy Izz];
end

function [I_pond] = calcMomInerciaGlobal(EPondVec,coordinateMatrix,centroide_pond,IyyIzz,Y0Z0,A,A_pond)
    Iy0y0 = 0; Iz0z0 = 0;
    for i = 1:size(coordinateMatrix,1)
        Iy0y0 = Iy0y0 + EPondVec(i)*(IyyIzz(i,1)+((Y0Z0(i,2)^2)*A(i)));
        Iz0z0 = Iz0z0 + EPondVec(i)*(IyyIzz(i,2)+((Y0Z0(i,1)^2)*A(i)));
    end
    Iyy_pond = Iy0y0 - ((centroide_pond(:,2)^2)*A_pond);
    Izz_pond = Iz0z0 - ((centroide_pond(:,1)^2)*A_pond);
    I_pond = [Iyy_pond Izz_pond];
end

function [A_pond,centroide_pond,I_pond] = PropriedadesPonderadas(coordinateMatrix,E0,E1,E2)
EPondVec = [E1; E1; E2; E2; E2; E2]/E0;

%Calculo da ?rea A dos ret?ngulos e da ?rea ponderada total;
[A, A_pond] = areaPond(EPondVec, coordinateMatrix)

%C?lculo do centroide*
[centroide_pond, Y0Z0] = calcCentroide(EPondVec, coordinateMatrix,A)

%Momentos de in?rcia locais
IyyIzz = calcMomInerciaLocal(coordinateMatrix)

%Momento de in?rcia globais ponderados
I_pond = calcMomInerciaGlobal(EPondVec,coordinateMatrix,centroide_pond,IyyIzz,Y0Z0,A,A_pond)
end


function [dist] = plotDistribuicaoSchrenk(L, b, lambda, label)
syms x real;
mili = 0.001;
welip = ((4*L)/(pi*b))*sqrt(1-(2*x/b).^2);
wtrap = ((2*L)/(b*(1+lambda)))*(1-2*x*(1-lambda)/b);
weff = (welip+wtrap)/2;
figure()
hold on
fplot(welip*mili,[0,b/2],'--m')
fplot(wtrap*mili,[0,b/2],'-.r')
fplot(weff*mili,[0,b/2],'-b')
xlabel('x [m]')
ylabel(label)
title('Compara??o entre as distribui??es')
legend('Dist. el?ptica', 'Dist. trapezoidal', 'Dist. de Schrenk')
dist = weff;
end

function V = distV(f,b)
syms x real
V = int(-f);
const_V = subs(V,x,b/2);
V = V - const_V;
end

function M = distM(m, V,fator,b)
syms x real
M = int(-m +(fator*V));
const_M = subs(M,x,b/2);
M = M - const_M;
end

function angEixoNeutro = EixoNeutro(My,Mz,I_pond)
syms x real;
tangente = (Mz*I_pond(:,1))/(My*I_pond(:,2));
ang_rad = atan(tangente);
angEixoNeutro=rad2deg(ang_rad);
end

function sigmaX = distTensao(My,Mz,I_pond,E,E0)
syms x real;
syms y real;
syms z real;
sigmaX=(E/E0)*((-y*Mz/I_pond(:,2))+(z*My/I_pond(:,1)));
end

function  deflexao= def(M,E0,I,fator)
syms x real
deflexao=int(int(fator*M/(E0*I)));
C_def = subs(deflexao,x,0);
deflexao = deflexao - C_def;
end

function plotDeflexao(deflexao,b,titulo,label)
syms x real;
figure()
hold on
fplot(deflexao,[0,b/2])
hold off
grid
caption1=sprintf('%s (x)',titulo);
title(caption1)
xlabel('x (m)')
caption2= sprintf('%s (x) em %s', titulo,label);
ylabel(caption2)
end

function plotDeformada(deflexao,b,titulo,label)
syms x
envergadura = 0;
deform=envergadura+deflexao;

figure()
hold on
fplot(envergadura,[0,b/2],'-k','linewidth',2)
fplot(deform,[0,b/2],'--b','linewidth',2)
hold off
grid on
caption1=sprintf('%s',titulo);
title(caption1)
xlabel('x (m)')
caption2=sprintf('%s (m)',label);
ylabel(caption2)
legenda=sprintf('Deformada em %s',label);
legend('Indeformada',legenda)
end

