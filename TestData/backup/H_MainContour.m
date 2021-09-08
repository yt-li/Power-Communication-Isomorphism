% This function is used to plot the wb-lambda contour: i.e., the
% trace of lambda when fixing wb to different values.

% Author(s): Yitong Li, Yunjie Gu

%%
clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

ColorRGB();

%% Enables
Enable_SaveFigure = 1;

%%
ColorStep = 4;
% ColorLower = [1,1,0]; ColorUpper = [1,0,0];       % yellow to red
ColorLower = [0,1,1]; ColorUpper = [0,0,1];       % light blue to dark blue
% GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStep)';
% GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStep)';
% GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStep)';
% ColorCell = {};
% for k = 1:ColorStep
%     ColorCell{k} = [GradRed(k),GradGreen(k),GradBlue(k)];
% end

Fbase = 60;
FreqLower = (sqrt(2)-1)*Fbase;
FreqUpper = 2*Fbase;

%% Plot
figure(1)

LineWidth = 1.5;
NumPoint = 200;
FigSize = [0.1 0.1 0.3 0.6];

set(gcf,'units','normalized','outerposition',FigSize);

% Plot zero point
plot(0, 0, 'k.', 'markersize', 20); grid on; hold on;

% Analysis region
AxisScale = 150;
x_vec = linspace(-AxisScale,0,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);

% When lambda = 0, wb = (sqrt(2)-1)*60Hz = 24.8528Hz, which is the lower
% limit for wb.
fb = (sqrt(2)-1)*Fbase;
FreqFactor = (fb-FreqLower)/(FreqUpper-FreqLower);
LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineWidth',LineWidth,'LineColor',LineColor); grid on; hold on;

fb = 2/3*Fbase;
FreqFactor = (fb-FreqLower)/(FreqUpper-FreqLower);
LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineWidth',LineWidth,'LineColor',LineColor); grid on; hold on;

fb = Fbase;
FreqFactor = (fb-FreqLower)/(FreqUpper-FreqLower);
LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineWidth',LineWidth,'LineColor',LineColor); grid on; hold on;

AxisScale = 250;
x_vec = linspace(-AxisScale,0,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);

fb = 1.5*Fbase;
FreqFactor = (fb-FreqLower)/(FreqUpper-FreqLower);
LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineWidth',LineWidth,'LineColor',LineColor); grid on; hold on;

fb = 2*60;
FreqFactor = (fb-FreqLower)/(FreqUpper-FreqLower);
LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineColor','k'); grid on; hold on;

xlabel('Real Part of $\lambda$ (Hz)','interpreter','latex')
ylabel('Imaginary Part of $\lambda$ (Hz)','interpreter','latex')
%xlim([-100,20]);
% xticks(t_Ticks)
%ylim([-100,150]);
% yticks(w_Ticks);

if Enable_SaveFigure
    print(gcf,'H_Contour.png','-dpng','-r600');
end

%% Plot pole
pole_G = load('pole_G').pole_G;
for k = 1:length(pole_G)
    scatter(real(pole_G{k}),imag(pole_G{k}),'x','LineWidth',1.5); hold on; grid on;
end
