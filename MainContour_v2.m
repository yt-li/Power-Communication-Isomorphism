% This function is used to plot the wb-lambda contour

% Author(s): Yitong Li, Yunjie Gu

%%
clear all
clc
close all

figure(1)

%
fb = 30;
NumPoint = 1000;
AxisScale = 100;
x_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1]); grid on; hold on;

%
fb = 50;
NumPoint = 1000;
AxisScale = 150;
x_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1]); grid on; hold on;
% contour(x_vec,y_vec,ConMat,[1,1],'Color','r'); grid on; hold on;

%
fb = 100;
NumPoint = 1000;
AxisScale = 200;
x_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = linspace(-AxisScale,AxisScale,NumPoint);
y_vec = flip(y_vec);
ConMat = GetContourMatrix(fb,x_vec,y_vec);
contour(x_vec,y_vec,ConMat,[1,1]); grid on; hold on;

stop

%% Output
Mag_F_Value;
Mag_F_Sign

Ang_F_Value1;
Ang_F_Value2;
Ang_F_Sign


F_Sign
F_Plot
