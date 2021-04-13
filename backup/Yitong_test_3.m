clear all
clc
close all

w_pll = 2*pi*10;
kp_pll = w_pll;
ki_pll = w_pll^2/4;

Ai = [0,        0;
      ki_pll,   0];
Bi = [1;
      kp_pll];
Ci = [ki_pll,   0;
      0,        1];
Di = [kp_pll;
      0];
T_I_ss = ss(Ai,Bi,Ci,Di);

T = append(T_I_ss,T_I_ss);

gamma_p = 1;
gamma_n = -1;

gamma = [gamma_p,       gamma_n;
          conj(gamma_n), conj(gamma_p)];
      
feedin = [1,2];
feedout = [1,3];
Gcl = feedback(T,gamma,feedin,feedout);

pole_sys = pole(Gcl)/2/pi;
      
