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

%% Enables
Enable_SaveFigure = 1;

ColorRGB();

%% Plot
figure(1)

FigSize = [0.1 0.1 0.3 0.6];
set(gcf,'units','normalized','outerposition',FigSize);

% Plot zero point
plot(0, 0, 'k.', 'markersize', 20); grid on; hold on;

% Notes: 
% When lambda = 0, wb = (sqrt(2)-1)*60Hz = 24.8528Hz, which is the lower
% limit for wb.

% Analysis region
NumPoint = 1000;
AxisScale = 150;
x_vec = linspace(-AxisScale,0,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);

fb = 24.8528;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineColor','k'); grid on; hold on;

fb = 40;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineColor','k'); grid on; hold on;

fb = 60;
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1],'LineColor','k'); grid on; hold on;

% %
% fb = 120;
% NumPoint = 1000;
% AxisScale = 250;
% x_vec = linspace(-AxisScale,0,NumPoint);
% y_vec = linspace(-AxisScale,AxisScale,NumPoint);
% y_vec = flip(y_vec);
% ConMat = GetContourMatrix(fb,x_vec,y_vec);
% contour(x_vec,y_vec,ConMat,[1,1]); grid on; hold on;

% t1 = quiver(-100,0,120,0,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.1,'AlignVertexCenters','on'); grid on; hold on;
% t2 = quiver(0,-100,0,250,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.1,'AlignVertexCenters','on'); grid on; hold on;

xlabel('Real Part of $\lambda$ (Hz)','interpreter','latex')
ylabel('Imaginary Part of $\lambda$ (Hz)','interpreter','latex')
xlim([-100,20]);
% xticks(t_Ticks)
ylim([-100,150]);
% yticks(w_Ticks);

if Enable_SaveFigure
    print(gcf,'H_Contour.png','-dpng','-r600');
end

