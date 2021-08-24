clear all
clc
close all

BLUW = [0, 0.4470, 0.7410];
RED = [0.8500, 0.3250, 0.0980];

% UserData = '3MachineModel_test_v3';
UserData = 'Nature_NETS_NYPS_68Bus_original';
LoadData();

Nmax1 = max(ListLine(:,1));
Nmax2 = max(ListLine(:,2));
Nmax = max(Nmax1,Nmax2);
N_Line = size(ListLine,1);

GraphMatrix = zeros(Nmax,Nmax);

for i = 1:N_Line
    m = ListLine(i,1);
    n = ListLine(i,2);
    if m ~= n
        GraphMatrix(m,n) = 1;
        GraphMatrix(n,m) = 1;
    end
end

figure(1)
GraphData = graph(GraphMatrix,'upper');
GraphPlot = plot(GraphData); grid on;

% highlight(GraphPlot,[1 2],'EdgeColor','k');
% highlight(GraphPlot,[1 2],'NodeColor',RED);