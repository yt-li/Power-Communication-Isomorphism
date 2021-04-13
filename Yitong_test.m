clear all
clc
close all

% s = tf('s');
s = sym('s','real');

W0 = 2*pi*50;

% AC filter
X = 0.05;
L = X/W0;
R = X/5;

% AC load
Yload = 0.5;
Zload = 1/Yload;

% Current loop
wi = 2*pi*250;
kpi = wi*L;
kii = wi^2*L/4;

PIi = kpi + kii/s;

% PLL controller
w_pll = 2*pi*10;
kp_pll = w_pll;
ki_pll = w_pll^2/4;

PIpll = kp_pll + ki_pll/s;

%% Impedance method
Zin = PIi + (s+1i*W0)*L + R;
Zload = 1/0.5;

Gcl = 1/(Zin+Zload);

if ~strcmpi(class(s),'sym')
pole_sys = pole(Gcl)/2/pi;
pole_sys(1)
pole_sys(2)
end

%% Yunjie's method
w = sym('w');
Zin_ab = kpi + kii/(s - 1i*w) + s*L + R;
Gbus = - 1/((1/Zin_ab) + Yload);
dGbus_dw = diff(Gbus,w)

dGbus_dw = subs(dGbus_dw,W0);
S = 0.5*0.5;
Gamma = real(-S*dGbus_dw*exp(-pi/2));
Gamma = simplify(Gamma)

