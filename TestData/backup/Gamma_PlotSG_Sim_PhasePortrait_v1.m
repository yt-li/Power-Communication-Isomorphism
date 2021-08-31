% This script plots the gamma and damping analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

RgbBlue = [0, 0.4470, 0.7410];
RgbRed = [0.8500, 0.3250, 0.0980];
RgbYellow = [0.9290, 0.6940, 0.1250];
RgbPurple = [0.4940, 0.1840, 0.5560];
RgbGreen = [0.4660, 0.6740, 0.1880];

%% Enables
Enable_SaveFigure = 0;

%% Load data
Gamma_SimData_SG = load('Gamma_SimData_SG_Unstable').Gamma_SimData_SG;

Time = Gamma_SimData_SG.time;
vdq_SG = Gamma_SimData_SG.signals(1).values;
w_SG_ = Gamma_SimData_SG.signals(5).values;
for i = 1:length(w_SG_)
    w(i) = w_SG_(:,:,i);
end
theta__ = Gamma_SimData_SG.signals(6).values;
for i = 1:length(theta__)
    theta(i) = theta__(:,:,i);
end

%% Calc
TimeStart = find(Time == 3);
TimeEnd = find(Time == 3.5);

w_ref = w(TimeStart-100);
theta_SG_ref = theta(TimeStart-100);

w = w(TimeStart+52000:TimeEnd) - w_ref;                     % The first few points with "noises" are discarded.
theta = theta(TimeStart+52000:TimeEnd) - theta_SG_ref;

% Normalize
W0 = 2*pi*60;
J = 0.05*2/W0;
K = 2.8253;            % Measured by runing MainSyn

if 0
    w = w*W0*sqrt(J);
    theta = theta*sqrt(K);
else
    w = w;
    theta = theta*sqrt(K)/sqrt(J)/W0;
end

%% Find critical point

toler = 5e-4;
SaveIndex = [];
Counter = 0;
for k = 1:length(w)
    if abs(abs(w(k)) - abs(theta(k)))<=toler
        Counter = Counter + 1;
        SaveIndex(Counter) = k;
    end
end
Index(1) = SaveIndex(1);
Counter = 1;
for k = 2:length(SaveIndex)
    if w(SaveIndex(k))*w(SaveIndex(k-1))<0
        Counter = Counter+1;
        Index(Counter) = SaveIndex(k);
    end
end
CircleTheta{1} = pi/4 : pi/100 : pi*3/4;
CircleRadius{1} = abs(w(Index(1)))*sqrt(2)-5e-4;

CircleX{1} = CircleRadius{1} * cos(CircleTheta{1});
CircleY{1} = CircleRadius(1) * sin(CircleTheta{1});
CircleRadius2 = abs(w(Index(2)))*sqrt(2)+3e-4;
CircleTheta2 = pi*3/4 : pi/100 : pi*5/4;
CircleX2 = CircleRadius2 * cos(CircleTheta2) + 0;
CircleY2 = CircleRadius2 * sin(CircleTheta2) + 0;
CircleRadius3 = abs(w(Index(3)))*sqrt(2)-5e-4;
CircleTheta3 = pi/4 : pi/100 : pi*3/4;
CircleX3 = CircleRadius3 * cos(CircleTheta3) + 0;
CircleY3 = CircleRadius3 * sin(CircleTheta3) + 0;

% toler = 1e-4;
% for k = 1:length(w)
%     if abs(w(k)) <= toler
%         Index(1) = k;
%         break;
%     end
% end
% for k = 1:length(theta)
%     if abs(theta(k)) <= toler
%         Index(2) = k;
%         break;
%     end
% end
% CircleRadius1 = abs(theta(Index(1)))-1e-4;
% CircleTheta1 = pi/2 : pi/100 : pi;
% CircleX1 = CircleRadius1 * cos(CircleTheta1);
% CircleY1 = CircleRadius1 * sin(CircleTheta1);
% CircleRadius2 = abs(w(Index(2)))+1e-5;
% CircleTheta2 = pi : pi/100 : pi*3/2;
% CircleX2 = CircleRadius2 * cos(CircleTheta2);
% CircleY2 = CircleRadius2 * sin(CircleTheta2);

%% Plot settings

FigNum = 0;
LineWidth = 1.2;

FigRowMax = 1;
FigColumnMax = 2;

Time = Time - 8;
x_Limit = [-0.03,0.03];
x_Ticks = [0,2,4,6,8,10];
y_Limit = x_Limit;
y_Ticks = x_Ticks;

FigSize1 = [0.1 0.1 0.4 0.4];
FigSize2 = [0.1 0.1 0.4 0.7];

%% Plot

if 1
FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.35 0.65]);

plot(w,theta,'LineWidth',LineWidth); grid on; hold on;
xlabel('Normalized $\Delta \omega$','interpreter','latex')
ylabel('Normalized $\Delta \theta$','interpreter','latex')
xlim(x_Limit);
% xticks(t_Ticks)
ylim(y_Limit);
% yticks(w_Ticks);

plot(w(Index(1)),theta(Index(1)),'.','MarkerSize',15,'Color',RgbRed); grid on; hold on;
plot([0,w(Index(1))],[0,theta(Index(1))],'Color',RgbRed); grid on; hold on;
plot(CircleX,CircleY,'.','Color',RgbRed,'LineWidth',1); grid on; hold on;

plot(w(Index(2)),theta(Index(2)),'.','MarkerSize',15,'Color',RgbYellow); grid on; hold on;
plot([0,w(Index(2))],[0,theta(Index(2))],'Color',RgbYellow); grid on; hold on;
plot(CircleX2,CircleY2,'.','Color',RgbYellow,'LineWidth',1); grid on; hold on;

if Enable_SaveFigure
    print(gcf,'Gamma_Sim_SG.png','-dpng','-r600');
end

end




% if 0
% FigNum = FigNum + 1;
% figure(FigNum)
% set(gcf,'units','normalized','outerposition',FigSize2);
% 
% subplot(FigRowMax,FigColumnMax,1)
% plot(w,theta,'LineWidth',LineWidth); grid on; hold on;
% xlabel('Normalized $\Delta \omega$','interpreter','latex')
% ylabel('Normalized $\Delta \theta$','interpreter','latex')
% % plot([-0.03,0.03],[-0.03,0.03]); grid on; hold on;
% % plot([-0.03,0.03],[0.03,-0.03]); grid on; hold on;
% plot(CircleX,CircleY,'-.','Color','r'); grid on; hold on;
% xlim(x_Limit);
% % xticks(t_Ticks)
% ylim(y_Limit);
% % yticks(w_Ticks);
% 
% subplot(FigRowMax,FigColumnMax,2)
% plot(w,theta,'LineWidth',LineWidth); grid on; hold on;
% xlabel('Normalized  $\Delta \omega$','interpreter','latex')
% ylabel('Normalized $\Delta \theta$','interpreter','latex')
% % plot([-0.03,0.03],[-0.03,0.03]); grid on; hold on;
% % plot([-0.03,0.03],[0.03,-0.03]); grid on; hold on;
% plot(CircleX2,CircleY2,'Color','r'); grid on; hold on;
% xlim(x_Limit);
% % xticks(t_Ticks)
% ylim(y_Limit);
% % yticks(w_Ticks);
% 
% end

