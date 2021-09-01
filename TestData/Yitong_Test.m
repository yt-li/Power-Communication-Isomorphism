clear all
clc
close all

% figure(1)
% plot([0,1],[0,1],'Color','r'); grid on; hold on;
% plot([0,1],[0,2],'Color','k'); grid on; hold on;
% h1 = colorbar
% caxis([20 50])
% 
% figure(2)
% 
% colormap(jet)


r = (0:.1:.9)';
g = r.^1.8;
b = r.^2.1;
mymap = [r g b]; 
h1 = colormap(mymap)
colorbar()
clear h1