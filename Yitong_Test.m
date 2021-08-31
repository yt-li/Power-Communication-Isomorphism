clear all
clc
close all

% figure(1)
% plot([0,10],[0,10]); grid on; hold on;
% x = [0 1];
% y = [0 1];
% annotation('textarrow',x,y)

figure(1)
plot([1:10]); grid on; hold on;
theta = pi/2/2;
t1 = [2,3,10,10];
a = quiver(2,3,cos(theta),sin(theta),1); grid on; hold on;
