clear all
clc
close all

s = tf('s');

Fbase = 60;
Wbase = Fbase*2*pi;

% L filter impedance
Lf = 0.05/Wbase;
Rf = 0.01;
Zf = s*Lf + Rf;

% C filter impedance
Cf = 0.0002/Wbase;
Yc = s*Cf;
Zc = 1/Yc;

% Line impedance
Ll = 0.3/Wbase;
Rl = 0.03;

Zl = s*Ll + Rl;

%% Calculate bandwidth for Voltage node
pole_w = pole(1/(Zl+Lf+Rf));
wb_v = H_CalcBandwidth(pole_w);
fb_v = wb_v/2/pi;

%%
% PI control impedance
w_i = 1000*2*pi;
kp_i = Lf*w_i;
ki_i = Lf*w_i^2/4;
Z_PIi = kp_i + ki_i/(s-1i*Wbase);

% Current node impedance
Zi = Z_PIi+Zf;
% Zi = Zi*Zc/(Zi+Zc);       % Ignore the parallel capacitance

% Channel gain
G = Zi / (Zi + Zl);
G = minreal(G);
pole_w = pole(G);
pole_f = pole_w/2/pi;
figure(1)
scatter(real(pole_f),imag(pole_f),'x','LineWidth',1.5); hold on; grid on;

%% Calculate bandwidth for current node
for k = 1:length(pole_w)
    wb_i_(k) = H_CalcBandwidth(pole_w(k));
end
wb_i = min(wb_i_);
fb_i = wb_i/2/pi;

%%
if 0
w_i = linspace(100,1000,5)*2*pi;

pole_f = {};
for k = 1:length(w_i)
    kp_i = Lf*w_i(k);
    ki_i = Lf*w_i(k)^2/4;
    Z_PIi = kp_i + ki_i/(s-1i*Wbase);
    G = (Z_PIi+Zf) / (Zf + Z_PIi + Zl);
    G = minreal(G);
    pole_f{k} = pole(G)/2/pi;
end
save('pole_G','pole_G');

figure(2)
for k = 1:length(w_i)
    scatter(real(pole_f{k}),imag(pole_f{k}),'x','LineWidth',1.5); hold on; grid on;
end

% Notes:
% By sweeping w_i, we could find that the larger of w_i, the larger of w_b.

end
