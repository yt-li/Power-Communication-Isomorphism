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
pole_T1cl = load('Gamma_SG_pole_D0d1_T1cl').pole_T1cl;
pole_T12cl = load('Gamma_SG_pole_D0d1_T12cl').pole_T12cl;

%% Plot settings

FigNum = 0;
LineWidth = 1.5;

x_Limit = [-3,6]*1e-3;
x_Ticks = [-3,0,3,6]*1e-3;
y_Limit = [-3,3]; 
y_Ticks = [0.996,1,1.004];

FigSize = [0.1 0.1 0.35 0.45];

%% Plot

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);
scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',LineWidth); hold on; grid on;
scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',LineWidth); hold on; grid on;
legend({'Without $\Gamma$','With $\Gamma$'},'interpreter','latex')
xlabel('Real Part (Hz)','interpreter','latex')
ylabel('Imagary Part (Hz)','interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
% yticks(y_Ticks);

if Enable_SaveFigure
    print(gcf,'Gamma_Pole_SG.png','-dpng','-r600');
end