clear all
clc
close all

pole_sys_toolbox = load('pole_sys').pole_sys;
pole_sys_toolbox_ = load('pole_sys_').pole_sys;

index = find(abs(imag(pole_sys_toolbox))<45);
pole_sys_toolbox = pole_sys_toolbox(index);
index = find(abs(real(pole_sys_toolbox))<20);
pole_sys_toolbox = pole_sys_toolbox(index);

index = find(abs(imag(pole_sys_toolbox_))<45);
pole_sys_toolbox_ = pole_sys_toolbox_(index);
index = find(abs(real(pole_sys_toolbox_))<20);
pole_sys_toolbox_ = pole_sys_toolbox_(index);

figure(999)
scatter(real(pole_sys_toolbox),imag(pole_sys_toolbox),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_sys_toolbox_),imag(pole_sys_toolbox_),'x','LineWidth',1.5); hold on; grid on;
legend('toolbox,w','toolbox,w0')
