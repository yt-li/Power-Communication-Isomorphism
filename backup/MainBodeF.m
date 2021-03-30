% This function plots the bode diagram of F

%%
clear all
clc
close all

%%
s = sym('s');
w0 = 2*pi*50;   % rad/s

lambda = [5,50,500]*2*pi;

for i = 1:length(lambda)
    F(i) = (1i*w0 - lambda(i))/(s + 1i*w0 - lambda(i));
end

%%
% Figure index
fn = 0;

% Set frequency
f_p = logspace(-1,3,1000);
f_pn = [-flip(f_p),f_p];
w_pn = f_pn*2*pi;
s_pn = 1i*w_pn;

%% Log-scaled axis
fn = fn+1;
figure(fn)

for i = 1:length(lambda)
    SimplexPS.bode_c(F(i),1j*w_pn,2*pi,'PhaseOn',1,'LineWidth',1.5);
end
legend('$\lambda$=5 Hz','$\lambda$=50 Hz','$\lambda$=500 Hz','interpreter','latex','location','southeast');

% Get subplot index
f11 = subplot(2,2,1);
f12 = subplot(2,2,2);
f21 = subplot(2,2,3);
f22 = subplot(2,2,4);

% Set axis
XLim_f_p = [1e-1,1e3];
XTick_f_p = logspace(-1,3,5);
set(f11,'XLim',-flip(XLim_f_p));
set(f21,'XLim',-flip(XLim_f_p));
set(f11,'XTick',-flip(XTick_f_p));
set(f21,'XTick',-flip(XTick_f_p));

set(f12,'XLim',XLim_f_p);
set(f22,'XLim',XLim_f_p);
set(f12,'XTick',XTick_f_p);
set(f22,'XTick',XTick_f_p);

YLim_mag = [1e-1,1e1];
set(f11,'YLim',YLim_mag);
set(f12,'YLim',YLim_mag);
set(f12,'YTickLabel',[]);

YLim_ang = [-180,180];
YTick_ang = [-180,-90,0,90,180];
set(f21,'YLim',YLim_ang);
set(f22,'YLim',YLim_ang);
set(f21,'YTick',YTick_ang);
set(f22,'YTick',YTick_ang);
set(f22,'YTickLabel',[]);

% Change position
x_shift = 0.03;
f11.Position(1) = f11.Position(1) + x_shift;
f21.Position(1) = f21.Position(1) + x_shift;

f12.Position(1) = f12.Position(1) - x_shift;
f22.Position(1) = f22.Position(1) - x_shift;

% Set label
f11.YLabel.String = 'Magnitude';
f21.YLabel.String = 'Phase (Degree)';
f21.XLabel.String = 'Negative Frequency (Hz)';
f22.XLabel.String = 'Positive Frequency (Hz)';


%% Linear axis
fn = fn+1;
figure(fn)

linewidth = 1;

% for i = 1:1
for i = 1:length(lambda)

    X = F(i);

    funcX = matlabFunction(X);
 
 	Xw = zeros(1,length(s_pn));
    for n = 1:length(s_pn)
        Xw(n) = funcX(s_pn(n));
    end
        
    ang_Xw = angle(Xw)/pi*180;    % Degree
    mag_Xw = abs(Xw);
    mag_Xw_log = 20*log10(mag_Xw);
    
    subplot(2,1,1)
    plot(f_pn,mag_Xw,'linewidth',1); grid on; hold on;
    % plot(f_pn,mag_Xw_log,'linewidth',1); grid on; hold on;
    ylabel('Magnitude')
    subplot(2,1,2)
    plot(f_pn,ang_Xw,'linewidth',1); grid on; hold on;
    ylabel('Phase (Degree)')
    xlabel('Frequency (Hz)')
    
end
legend('$\lambda$=5 Hz','$\lambda$=50 Hz','$\lambda$=500 Hz','interpreter','latex','location','southeast')