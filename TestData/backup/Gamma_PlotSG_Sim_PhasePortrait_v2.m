% This script plots the gamma and damping analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

RgbBlue = [0, 0.4470, 0.7410];
RgbRed = [0.8500, 0.3250, 0.0980];
RgbYellow = [0.9290, 0.6940, 0.1250];
RgbPurple = [0.4940, 0.1840, 0.5560];
RgbGreen = [0.4660, 0.6740, 0.1880];

%% Enables
Enable_SaveFigure = 0;

%% Load data
Gamma_SimData_SG = load('Gamma_SimData_SG_Unstable').Gamma_SimData_SG;

Time = Gamma_SimData_SG.time;
vdq_SG = Gamma_SimData_SG.signals(1).values;
w_SG_ = Gamma_SimData_SG.signals(5).values;
for i = 1:length(w_SG_)
    w(i) = w_SG_(:,:,i);
end
theta__ = Gamma_SimData_SG.signals(6).values;
for i = 1:length(theta__)
    theta(i) = theta__(:,:,i);
end

%% Calc
TimeRef = find(Time == 3);
TimeStart = find(Time == 3.35);
TimeEnd = find(Time == 3.5);

w_ref = w(TimeRef-100);
theta_SG_ref = theta(TimeRef-100);

% w = w(TimeStart+52000:TimeEnd) - w_ref;                     % The first few points with "noises" are discarded.
% theta = theta(TimeStart+52000:TimeEnd) - theta_SG_ref;
w = w(TimeStart:TimeEnd) - w_ref;                     % The first few points with "noises" are discarded.
theta = theta(TimeStart:TimeEnd) - theta_SG_ref;

% Normalize
if 0                        % We do not normalize the results by default
    W0 = 2*pi*60;
    J = 0.05*2/W0;
    K = 2.8253;            % Measured by runing MainSyn
    if 0
        w = w*W0*sqrt(J);
        theta = theta*sqrt(K);
    else
        w = w;
        theta = theta*sqrt(K)/sqrt(J)/W0;
    end
end

%% Plot settings

FigNum = 0;
LineWidth = 1.2;

FigRowMax = 1;
FigColumnMax = 2;

x_Limit = [-0.04,0.04];
x_Ticks = [-0.04,-0.02,0,0.02,0.04];
y_Limit = [-0.12,0.12];
y_Ticks = [-0.12,-0.06,0,0.06,0.12];

%% Plot

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.35 0.63]);

% Unstable phase portrait
plot(w,theta,'LineWidth',LineWidth); grid on; hold on;
xlabel('$\Delta \omega$ (pu)','interpreter','latex')
ylabel('$\Delta \theta$ (rad)','interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

% Unstable arrows
StepNum = 100;
StepSize = floor(length(w)/StepNum);
for i = 1:(StepNum-1)
    w1 = w(i*StepSize+1);
    theta1 = theta(i*StepSize+1);
    ScaleFactor(i) = abs(w1)^2*370;
    quiver(w1,theta1,w1*1.01,theta1*1.01,ScaleFactor(i),'Color','k','MaxHeadSize',0.8); grid on; hold on;
end

if Enable_SaveFigure
    print(gcf,'Gamma_Sim_SG.png','-dpng','-r600');
end





