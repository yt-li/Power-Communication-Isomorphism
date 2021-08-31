% This script plots the gamma and damping analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 0;

%% Load data
Gamma_SimData_SG = load('Gamma_SimData_SG').Gamma_SimData_SG;

Time = Gamma_SimData_SG.time;

vdq_SG = Gamma_SimData_SG.signals(1).values;
w_SG_ = Gamma_SimData_SG.signals(5).values;
for i = 1:length(w_SG_)
w_SG(i) = w_SG_(:,:,i);
end
pq_SG = Gamma_SimData_SG.signals(7).values;
s_SG = pq_SG(:,1) + 1i*pq_SG(:,2);
s_SG_abs = abs(s_SG);

%% Plot settings

FigNum = 0;
LineWidth = 1;

FigRowMax = 1;
FigColumnMax = 2;

Time = Time - 8;
t_Limit = [0,10];
t_Ticks = [0,2,4,6,8,10];
w_Limit = [0.996,1.004]; 
w_Ticks = [0.996,1,1.004];
s_Limit = [0,1];
s_Ticks = [0,0.5,1];

FigSize = [0.1 0.1 0.35 0.28];

%% Plot

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);
subplot(FigRowMax,FigColumnMax,1)
plot(Time,w_SG,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(w_Limit);
yticks(w_Ticks);

subplot(FigRowMax,FigColumnMax,2)
plot(Time,s_SG_abs,'LineWidth',LineWidth); grid on; hold on;
ylabel('$S$ (pu)','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(s_Limit);
yticks(s_Ticks);

if Enable_SaveFigure
    print(gcf,'Gamma_Sim_SG.png','-dpng','-r600');
end