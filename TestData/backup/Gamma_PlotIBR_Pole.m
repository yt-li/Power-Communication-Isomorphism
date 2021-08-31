% This script plots the gamma and damping analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

%% Load data
pole_T1cl{1}    = load('Gamma_IBR_pole_PnQp_T1cl').pole_T1cl;
pole_T12cl{1}   = load('Gamma_IBR_pole_PnQp_T12cl').pole_T12cl;

pole_T1cl{2}    = load('Gamma_IBR_pole_P0Qp_T1cl').pole_T1cl;
pole_T12cl{2}   = load('Gamma_IBR_pole_P0Qp_T12cl').pole_T12cl;

pole_T1cl{3}    = load('Gamma_IBR_pole_PpQp_T1cl').pole_T1cl;
pole_T12cl{3}   = load('Gamma_IBR_pole_PpQp_T12cl').pole_T12cl;

pole_T1cl{4}    = load('Gamma_IBR_pole_PnQ0_T1cl').pole_T1cl;
pole_T12cl{4}   = load('Gamma_IBR_pole_PnQ0_T12cl').pole_T12cl;

pole_T1cl{5}    = load('Gamma_IBR_pole_P0Q0_T1cl').pole_T1cl;
pole_T12cl{5}   = load('Gamma_IBR_pole_P0Q0_T12cl').pole_T12cl;

pole_T1cl{6}    = load('Gamma_IBR_pole_PpQ0_T1cl').pole_T1cl;
pole_T12cl{6}   = load('Gamma_IBR_pole_PpQ0_T12cl').pole_T12cl;

pole_T1cl{7}    = load('Gamma_IBR_pole_PnQn_T1cl').pole_T1cl;
pole_T12cl{7}   = load('Gamma_IBR_pole_PnQn_T12cl').pole_T12cl;

pole_T1cl{8}    = load('Gamma_IBR_pole_P0Qn_T1cl').pole_T1cl;
pole_T12cl{8}   = load('Gamma_IBR_pole_P0Qn_T12cl').pole_T12cl;

pole_T1cl{9}    = load('Gamma_IBR_pole_PpQn_T1cl').pole_T1cl;
pole_T12cl{9}   = load('Gamma_IBR_pole_PpQn_T12cl').pole_T12cl;

% Notes:
% p means +0.5, n means -0.5, 0 means 0

%% Plot settings

FigNum = 0;
LineWidth = 1.5;

FigRowMax1 = 2;
FigColumnMax1 = 2;

FigRowMax2 = 3;
FigColumnMax2 = 3;

x_Limit = [-6,0];
x_Ticks = [-6,-4,-2,0];
y_Limit = [-1,1]; 
y_Ticks = [-1,-0.5,0,0.5,1];

FigSize = [0.1 0.1 0.35 0.45];
FigSize2 = [0.1 0.1 0.6 0.8];

FigTitle{1} = '$P=-0.5$ pu, $Q=0.5$ pu';
FigTitle{2} = '$P=0$ pu, $Q=0.5$ pu';
FigTitle{3} = '$P=0.5$ pu, $Q=0.5$ pu';
FigTitle{4} = '$P=-0.5$ pu, $Q=0$ pu';
FigTitle{5} = '$P=0$ pu, $Q=0$ pu';
FigTitle{6} = '$P=0.5$ pu, $Q=0$ pu';
FigTitle{7} = '$P=-0.5$ pu, $Q=-0.5$ pu';
FigTitle{8} = '$P=0$ pu, $Q=-0.5$ pu';
FigTitle{9} = '$P=0.5$ pu, $Q=-0.5$ pu';

%% Plot

FigNum = FigNum + 1;
figure(FigNum)
% set(gcf,'units','normalized','outerposition',FigSize);

subplot(FigRowMax1,FigColumnMax1,1)
scatter(real(pole_T1cl{1}),imag(pole_T1cl{1}),'x','LineWidth',LineWidth); hold on; grid on;
scatter(real(pole_T12cl{1}),imag(pole_T12cl{1}),'x','LineWidth',LineWidth); hold on; grid on;
legend({'Without $\Gamma$','With $\Gamma$'},'interpreter','latex')
xlabel('Real Part (Hz)','interpreter','latex')
ylabel('Imagary Part (Hz)','interpreter','latex')
title(FigTitle{1},'interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

subplot(FigRowMax1,FigColumnMax1,2)
scatter(real(pole_T1cl{3}),imag(pole_T1cl{3}),'x','LineWidth',LineWidth); hold on; grid on;
scatter(real(pole_T12cl{3}),imag(pole_T12cl{3}),'x','LineWidth',LineWidth); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex')
ylabel('Imagary Part (Hz)','interpreter','latex')
title(FigTitle{3},'interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

subplot(FigRowMax1,FigColumnMax1,3)
scatter(real(pole_T1cl{7}),imag(pole_T1cl{7}),'x','LineWidth',LineWidth); hold on; grid on;
scatter(real(pole_T12cl{7}),imag(pole_T12cl{7}),'x','LineWidth',LineWidth); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex')
ylabel('Imagary Part (Hz)','interpreter','latex')
title(FigTitle{7},'interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

subplot(FigRowMax1,FigColumnMax1,4)
scatter(real(pole_T1cl{9}),imag(pole_T1cl{9}),'x','LineWidth',LineWidth); hold on; grid on;
scatter(real(pole_T12cl{9}),imag(pole_T12cl{9}),'x','LineWidth',LineWidth); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex')
ylabel('Imagary Part (Hz)','interpreter','latex')
title(FigTitle{9},'interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

if Enable_SaveFigure
    print(gcf,'Gamma_Pole_IBR_2t2.png','-dpng','-r600');
end

% 3*3 plot
FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize2);
k = 0;
for r = 1:3
    for c = 1:3
        k = k + 1;
        subplot(FigRowMax2,FigColumnMax2,k)
        scatter(real(pole_T1cl{k}),imag(pole_T1cl{k}),'x','LineWidth',LineWidth); hold on; grid on;
        scatter(real(pole_T12cl{k}),imag(pole_T12cl{k}),'x','LineWidth',LineWidth); hold on; grid on;
        if r == 3
        xlabel('Real Part (Hz)','interpreter','latex')
        end
        if c == 1
        ylabel('Imagary Part (Hz)','interpreter','latex')
        end
        title(FigTitle{k},'interpreter','latex')
        xlim(x_Limit);
        xticks(x_Ticks)
        ylim(y_Limit);
        yticks(y_Ticks);
    end
end
legend({'Without $\Gamma$','With $\Gamma$'},'interpreter','latex')

if Enable_SaveFigure
    print(gcf,'Gamma_Pole_IBR_3t3.png','-dpng','-r600');
end

