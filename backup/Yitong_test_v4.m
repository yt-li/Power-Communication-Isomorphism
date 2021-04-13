clear all
clc
close all

s = sym('s');

W0 = 2*pi*50;

wi = 2*pi*250;
L = 0.05/W0;
R = 0.01;
kp = wi*L;
ki = wi^2*L/4;

Zin = kp + ki/(s-1i*W0) + s*L + R;
% Zin = (s-1i*W0)/((kp + s*L + R)*(s-1i*W0) + ki);
Yload = 0.5;
Ybus = 1/Zin + Yload;

% Method 1
% Gbus = -1/Ybus;
% Gbus_prime = diff(Gbus);
% Gbus = subs(Gbus,'s',1i*W0*(1+1e-10));
% Gbus_prime = subs(Gbus_prime,'s',1i*W0*(1+1e-10));

% Method 2
Gbus = -((kp+R+s*L)*(s-1i*W0)+ki)/( (s-1i*W0) + 0.5*((kp+R+s*L)*(s-1i*W0)+ki) );
Gbus_prime = diff(Gbus);
Gbus = subs(Gbus,'s',1i*W0);
Gbus_prime = subs(Gbus_prime,'s',1i*W0);

% Method 3
% Gbus = -1/Ybus;
% dW = W0*1e-5;
% Gbus_ = subs(Gbus,'s',1i*(W0+dW));
% Gbus = subs(Gbus,'s',1i*W0*(1+1e-10));
% Gbus_prime = (Gbus_ - Gbus)/(1i*dW);

% Gbus_prime = subs(Gbus_prime,'s',1i*W0*(1+1e-10));
Gbus = double(Gbus)
Gbus_prime = double(Gbus_prime)