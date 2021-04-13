% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The device can not be floating bus, i.e., each bus should have a active
% device.
%
% PLi and QLi should be combined into the nodal admittance matrix.
%
% Devices should be listed in the same order as the bus.
%
% The flux inductor of SGs should be combined into the line impedances.
%
% Put voltage source and current source together.
%
% The reference direction should be doubled checked: Ybus, Gbus, T, 
%
% Add star-Y conversion

%% Clear
clear all
clc
close all

%% Load data
fprintf('Loading data...\n')

% For random test
% UserData = 'ExamplePowerSystem_v1_1.xlsx';    % large X/R ratio
% UserData = 'ExamplePowerSystem_v1_2.xlsx';    % small X/R ratio
% UserData = 'ExamplePowerSystem_v1_3.xlsx';      % With current node

% Pure SG system
% UserData = 'ExamplePowerSystem_v2.xlsx';      % pure SG system
UserData = 'ExamplePowerSystem_v2_1.xlsx';      % pure SG system, two machine
% UserData = 'ExamplePowerSystem_v2_2.xlsx';      % pure SG system, three machine

% Hybrid SG-IBR system
% UserData = 'ExamplePowerSystem_v3.xlsx';      % Hybrid SG IBR system

% ### Re-arrange basic settings
ListBasic = xlsread(UserData,'Basic');

Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;    % (A), base current
Zbase = Sbase/Ibase;    % (Ohm), base impedance
Ybase = 1/Zbase;        % (S), base admittance
Wbase = Fbase*2*pi;     % (rad/s), base angular frequency

% ### Re-arrange advanced settings
ListAdvance = xlsread(UserData,'Advance');
Flag_PowerFlowAlgorithm   	= ListAdvance(5);
Enable_CreateSimulinkModel	= ListAdvance(6);
Enable_PlotPole           	= ListAdvance(7);
Enable_PlotAdmittance     	= ListAdvance(8);
Enable_PrintOutput       	= ListAdvance(9);
Enable_Participation        = ListAdvance(10);

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplexPS.Toolbox.RearrangeListBus(UserData);

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplexPS.Toolbox.RearrangeListLine(UserData,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the device netlist
[DeviceBus,DeviceType,Para,N_Device] = SimplexPS.Toolbox.RearrangeListDevice(UserData,Wbase,ListBus);

% Get the device source type: 1-voltage source; 2-current source
for i = 1:N_Bus
    if 0<=DeviceType{i} && DeviceType{i}<=9
      	DeviceSourceType{i} = 1;
    elseif 10<=DeviceType{i} && DeviceType{i}<=19
    	DeviceSourceType{i} = 2;
    elseif 20<=DeviceType{i} && DeviceType{i}<=29
        DeviceSourceType{i} = 1;
    elseif DeviceType==90
        DeviceSourceType{i} = 1;  % An inf bus is regarded as voltage source.
    elseif DeviceType==100
        DeviceSourceType{i} = 2;  % Floating bus, i.e., no device, is regarded as a current source with i=0.
    else
     	error(['Error']);
    end
end

% Doing the delta-star conversion

%% Power Flow
fprintf('Doing power flow anlaysis...\n')

% ### Power flow analysis
[PowerFlow,~,V,I,Ang,P,Q,Vm] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
% Form of PowerFlow{i}: P, Q, V, xi, w

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);
ListPowerFlow_ = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

% Notes:
% Because V and I rather than PowerFlow are used later, the Load2SelfBranch
% will not work here. We have to make sure PLi and QLi = 0 in the netlist.

%% Network Matrix
fprintf('Calculating network matrix...\n')

% Calculate nodal admittance matrix
Ybus = YbusCalc_s(ListLine,Wbase,'albe');

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. The good news is that, the
% passive component, and the inner loops of inverters, are indeed symmtric
% in complex dq and alpha/beta frame.

% Seperate Ybus into four blocks according to source types
Ibus_n = N_Bus+1;   % Default: no ibus
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        Ibus_n = i;
        break;
    end
end
[Y11,Y12,Y21,Y22] = BlockMatrix(Ybus,Ibus_n-1,Ibus_n-1);

% Notes:
% It should be ensured that the buses are listed in the form like this:
% [Vbus1, Vbus2, Vbus3, Ibus4, Ibus5]

% Get Gbus
G11 = Y11 - Y12*inv(Y22)*Y21;
G12 = -Y12*inv(Y22);
G21 = inv(Y22)*Y21;
G22 = inv(Y22);

Gbus = [G11,G12;
        G21,G22];
    
Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.
    
% Get Gbus_prime
Gbus_prime = diff(Gbus,'s');

% Convert sym to num
Gbus = subs(Gbus,'s',1i*Wbase);
Gbus = double(Gbus);

Gbus_prime = subs(Gbus_prime,'s',1i*Wbase);
Gbus_prime = double(Gbus_prime);
    
% Update input and output
In = [V(1:Ibus_n-1);
      I(Ibus_n:end)];
Out = [I(1:Ibus_n-1);
       V(Ibus_n:end)];
   
% Get S matrix
S = conj(In)*transpose(In);

% Get epsilon
for i = 1:N_Bus
    if DeviceSourceType{i} == 1
        epsilon(i) = 0;
    elseif DeviceSourceType{i} == 2
        epsilon(i) = pi/2;
    else
        error(['Error']);
    end
end

% Get K matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        if n ~= m
            ang_K(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus(m,n));
            K(m,n) = abs(S(m,n))*abs(Gbus(m,n))*sin(ang_K(m,n));
        end
    end
end
for m = 1:N_Bus
    K_temp = 0;
    for n = 1:N_Bus
        if n~=m
            K_temp = K_temp - K(m,n);
        end
    end
    K(m,m) = K_temp;
end
K = vpa(K);
K = -K;             % For negative feedback
[K11,K12,K21,K22] = BlockMatrix(K,Ibus_n-1,Ibus_n-1);

% Get Gamma matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_Gamma(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        Gamma(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_Gamma(m,n));
    end
end
Gamma = vpa(Gamma);
Gamma = -Gamma;    % For negative feedback

%% Apparatus Matrix
fprintf('Calculating apparatus matrix...\n')

% Symbolic
s = sym('s');

% Synchronous generator
D = 5/Wbase^2;
J = 3.5*2/Wbase^2;

D = D*Wbase;
J = J*Wbase;
% Notes: Adding '*Wbase' is because the power swing equation rather than
% the torque swing equation is used, and P=T*w0 if w is not in per unit
% system.

% Droop
m = 1/D;
Tf = J*m;

% PI
w_pll = 190*2*pi;
kp_pll = w_pll;
ki_pll = w_pll^2/4;
PI_pll = kp_pll + ki_pll/s;
    
T_V = 1/(D + J*s);
F_V = -s/T_V;

T_I = PI_pll;
F_I = -s/T_I;

Hinv = eye(N_Bus);      % Let Hinv be identity matrix

%% Stability criterion
fprintf('Calculating stability criterion...\n')
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = BlockMatrix(KH,Ibus_n-1,Ibus_n-1);
[phi,xi] = eig(KH);
Gamma_Hphi = inv(phi)*Hinv*Gamma*phi;
[GH11,GH12,GH21,GH22] = BlockMatrix(Gamma_Hphi,Ibus_n-1,Ibus_n-1);
[~,sigma,~] = svd(Gamma_Hphi);
sigma_max = max(max(sigma));
for i = 1:length(xi)
    if xi(i,i)<-1e-9
        error(['Error']);
    end
end

% Current node
KH_I = K22;
GH_I = GH22;
[phi_I,xi_I] = eig(KH_I);
[~,sigma_I,~] = svd(GH_I);
sigma_I_max = max(max(sigma_I));
[zeta_m_I,w_min_I] = CalcZeta(T_I,diag(xi_I));

% Voltage node
if isempty(KH22)
    KH_V = KH;
    GH_V = Gamma_Hphi;
else
    KH_V = KH11+KH12*inv(KH22)*KH21;
    GH_V = GH11+KH12*inv(KH22)*GH21;
end
[phi_V,xi_V] = eig(KH_V);
[~,sigma_V,~] = svd(GH_V);
sigma_V_max = max(max(sigma_V));
[zeta_m_V,w_min_V] = CalcZeta(T_V,diag(xi_V));

%% Plot

figure_n = 0;

w_p = logspace(-1,1,500)*2*pi;
w_pn = [-flip(w_p),w_p];
s_pn = 1i*w_pn;

figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_V,s_pn);
scatter(real(diag(xi_V)),imag(diag(xi_V)),'x','LineWidth',1.5); hold on; grid on;

figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_I,s_pn);
scatter(real(diag(xi_I)),imag(diag(xi_I)),'x','LineWidth',1.5); hold on; grid on;

%% Output

xi
sigma

xi_I
sigma_I
zeta_m_I

xi_V
sigma_V
zeta_m_V

fprintf('Check if the system is stable by the proposed criterion:\n')
if isempty(zeta_m_I)
    FlagStable = min(zeta_m_V)>sigma_V_max;
else
    FlagStable = (min(zeta_m_I)>sigma_I_max) && (min(zeta_m_V)>sigma_V_max);
end
if FlagStable
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end


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

T1cl = inv(eye(N_Bus) + TH*K/s)*TH;
T1cl = minreal(T1cl,1e-5);
T2cl = inv(eye(N_Bus) + T1cl*Gamma);
T2cl = minreal(T2cl,1e-5);

pole_sys = pole(T2cl)/2/pi;

% Check stability
fprintf('Checking if the system is stable by the poles:\n')
if isempty(find(pole_sys>1e-8, 1))
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end

% Plot pole map
fprintf('Plotting pole map...\n')
figure_n = figure_n+1;

figure(figure_n);
scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Global pole map');
   
stop

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