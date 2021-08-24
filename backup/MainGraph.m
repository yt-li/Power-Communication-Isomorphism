clear all
clc
close all

BLUW = [0, 0.4470, 0.7410];
RED = [0.8500, 0.3250, 0.0980];


figure(1)
A = [1,2,10;
     2,0,0;
     10,0,0];
 
t1 = size(A,1)
NodeName = {'1','3','5'};
GraphData = graph(A,NodeName,'upper');
GraphPlot = plot(GraphData); grid on;

highlight(GraphPlot,[1 2],'EdgeColor','k');
highlight(GraphPlot,[1 2],'NodeColor',RED);