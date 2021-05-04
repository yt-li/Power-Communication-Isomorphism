% This function is a backup for the main function

%% Calculate transfer function to compare the poles with toolbox

% Notes:
% The two machine system coincides with the toolbox well. But when the
% number of machines is larger than two, the poles calculated here are
% wrong. This might be caused by the numerical calculation of matlab.

s = tf('s');

TH = [];
for i = 1:N_Bus
    if DeviceType{i} == 0               % Synchronous generator
        THinv_i(i) = 1/(D + J*s);
    elseif DeviceType{i} == 11          % Grid-following inverter
        PI_pll = kp_pll + ki_pll/s;
        THinv_i(i) = PI_pll;
    else
        error(['Error']);
    end
    TH = blkdiag(TH,THinv_i(i));
end

K = double(K);
Gamma = double(Gamma);

% T1cl = inv(eye(N_Bus) + TH*K/s)*TH;
% T1cl = minreal(T1cl,1e-5);
% T2cl = inv(eye(N_Bus) + T1cl*Gamma);
% T2cl = minreal(T2cl,1e-5);
% 
% pole_sys = pole(T2cl)/2/pi;
% 
% % Check stability
% fprintf('Checking if the system is stable by the poles:\n')
% if isempty(find(pole_sys>1e-8, 1))
%     fprintf('Stable!\n')
% else
%     fprintf('Unstable!\n')
% end
% 
% % Plot pole map
% fprintf('Plotting pole map...\n')
% figure_n = figure_n+1;
% 
% figure(figure_n);
% scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
% xlabel('Real Part (Hz)');
% ylabel('Imaginary Part (Hz)');
% title('Global pole map');


%%
% ===========================================
% For test
% ===========================================

if 1

T1 = 1/(D+J*s);
T2 = 1/D;
T3 = 1/(J*s);

F1 = -s/T1;
F2 = -s/T2;
F3 = -s/T3;
    
figure_n = figure_n+1;
figure(figure_n)
[p1,~] = SimplexPS.nyquist_c(F1,s_pn,'Color','b');
[p2,~] = SimplexPS.nyquist_c(F2,s_pn,'Color','r');
[p3,~] = SimplexPS.nyquist_c(F3,s_pn,'Color','y');
scatter(real(xi_diag),imag(xi_diag),'x','LineWidth',1.5); hold on; grid on;
legend([p1,p2,p3],{'1','2','3'})
    
end