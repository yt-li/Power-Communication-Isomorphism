% This function analyzes the xi

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

%% Load data
DataName = 'K_68Bus_SG_IBR_Load_Data';
% DataName = 'K_68Bus_SG_IBR_Data';
% DataName = 'K_68Bus_SG_IBR_17_Data';

Data = load(DataName).SaveData;

%%
KH              = Data.KH;
YbusVI          = Data.YbusVI;
YbusVIF         = Data.YbusVIF;
GbusVI          = Data.GbusVI;
GbusVIF         = Data.GbusVIF;
YbusOrigin      = Data.YbusOrigin;
Index_Vbus      = Data.Index_Vbus;
Index_Ibus      = Data.Index_Ibus;
Index_Fbus      = Data.Index_Fbus;
Index_Ebus      = Data.Index_Ebus;
Order_Old2New   = Data.Order_Old2New;
Order_New2Old   = Data.Order_New2Old;

NodeYZ = abs(diag(GbusVIF));
GraphFigure.NodeFontSize = 5.5;

for k = 1:length(Index_Vbus)
    NodeYZ(k) = 1/NodeYZ(k);    % Convert admittance to impedance at voltage node
end
NodeYZ = NodeYZ(Order_New2Old); % Convert the order back to its origin

%%
FigNum = 0;
ColorRGB();
FigSize = [0.1 0.1 0.5 0.75];

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);

%% Plot graph
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor',[0,0,0],'LineWidth',1.1);     % Change all edges to black by default
highlight(GraphFigure,GraphData,'NodeColor',[0,0,0]);                   % Change all nodes to black by default
highlight(GraphFigure,GraphData,'MarkerSize',5);

%% Set voltage node
highlight(GraphFigure,Index_Vbus,'NodeColor',[0,0,0]);

%% Set floating node
highlight(GraphFigure,Index_Fbus,'NodeColor',[0.7,0.7,0.7]);

%% Reduce pure empty node
highlight(GraphFigure,Index_Ebus,'MarkerSize',1);
% highlight(GraphFigure,Index_Fbus,'Marker','o');

%% Calculation
[Phi,Xi,Psi] = eig(KH);
PhiInv = inv(Phi);
Xi = diag(Xi);
[~,Index_XiMin] = min(real(Xi));
XiMin = Xi(Index_XiMin);
if abs(XiMin)<=1e-4                    % Check if xi_min is zero
    Xi_ = Xi;
    Xi_(Index_XiMin) = inf;
    [~,Index_XiMin] = min(real(Xi_));
    XiMin = Xi(Index_XiMin);
end

PhiRightMin = Phi(:,Index_XiMin);
PhiLeftMin = transpose(PhiInv(Index_XiMin,:));
FiedlerVec = PhiRightMin.*PhiLeftMin;
FiedlerAbs = abs(FiedlerVec);
% FiedlerAbsMax = max(FiedlerAbs);
FiedlerAbsMax = 0.7;
FiedlerAbsNorm = FiedlerAbs/FiedlerAbsMax;

%% Color bar
ColorStepSize = 100;
% ColorLower = [1,0.95,0];
% ColorUpper = [1,0,0];
ColorLower = [0,1,1];
ColorUpper = [0,0,1];
GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStepSize)';
GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStepSize)';
GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStepSize)';
% colormap([GradRed GradGreen GradBlue]);
% colorbar();
% caxis([0 FiedlerAbsMax]);

%% Set current node
highlight(GraphFigure,Index_Ibus,'NodeColor',[0,1,0]);
ColorFactor = FiedlerAbsNorm((length(Index_Vbus)+1):end);
for k = 1:length(Index_Ibus)
    NodeColor{k} = (ColorUpper - ColorLower) * ColorFactor(k) + ColorLower;
    highlight(GraphFigure,Index_Ibus(k),'NodeColor',NodeColor{k});
end

%% Set node label
for k = 1:length(Order_Old2New)
    GraphFigure.NodeLabel{k} = num2str(NodeYZ(k),2);
end

%%
% x = GraphFigure.XData';
% y = GraphFigure.YData';
% v = NodeYZ;
% 
% PlotHeatMap(x,y,v,1,[0,0.6]);

%% Save
if Enable_SaveFigure
    print(gcf,['Graph_' DataName '.png'],'-dpng','-r600');
end