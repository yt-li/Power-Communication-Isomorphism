% This function plots the graph for a power system

% Set the RGB of the default color in matlab
BLUE = [0, 0.4470, 0.7410];
RED = [0.8500, 0.3250, 0.0980];

% Plot the graph
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;

% Change all edges and nodes to black first
highlight(GraphFigure,GraphData,'EdgeColor','k','LineWidth',1);
highlight(GraphFigure,GraphData,'NodeColor','k');

% Highlight the node types by colors
IndexVoltageNode = find(DeviceSourceType == 1);
IndexCurrentNode = find(DeviceSourceType == 2);
IndexFloatingNode = find(DeviceSourceType == 3);

highlight(GraphFigure,IndexVoltageNode,'NodeColor',BLUE);
highlight(GraphFigure,IndexCurrentNode,'NodeColor',RED);
% highlight(GraphPlot,IndexFloatingNode,'NodeColor','k');