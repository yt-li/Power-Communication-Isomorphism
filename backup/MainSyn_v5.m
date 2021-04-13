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
% The current node has to be active load for ensuring K > 0.

%% Clear
clear all
clc
close all

%% Select data
% For random test
% UserData = 'ExamplePowerSystem_v1_1.xlsx';    % large X/R ratio
% UserData = 'ExamplePowerSystem_v1_2.xlsx';    % small X/R ratio
% UserData = 'ExamplePowerSystem_v1_3.xlsx';      % With current node

% Pure SG system
UserData = 'ExamplePowerSystem_v2.xlsx';      % pure SG system, 5 machines
% UserData = 'ExamplePowerSystem_v2_1.xlsx';      % pure SG system, 2 machine
% UserData = 'ExamplePowerSystem_v2_2.xlsx';      % pure SG system, 3 machine

% Hybrid SG-IBR system
% UserData = 'ExamplePowerSystem_v3.xlsx';    	% SG + IBR
% UserData = 'ExamplePowerSystem_v3_1.xlsx';    % SG + IBR + floating bus
% UserData = 'ExamplePowerSystem_v3_2.xlsx';  	% SG + floating bus
% UserData = 'ExamplePowerSystem_v3_3.xlsx';      % One SG + One IBR + No load
% UserData = 'ExamplePowerSystem_v3_4.xlsx';      % One SG + One IBR

% Nature example
%UserData = 'NETS_NYPS_68Bus.xlsx';  % Original
%UserData = 'NatureExample_NETS_NYPS_68Bus_v1.xlsx'; % Pure SG
% UserData = 'NatureExample_NETS_NYPS_68Bus_v2.xlsx'; % Hybrid

%% Load data
fprintf('Loading data...\n')

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

% ### Re-arrange the device netlist
[DeviceBus,DeviceType,DevicePara,N_Device] = SimplexPS.Toolbox.RearrangeListDevice(UserData,Wbase,ListBus);

% Notes:
% The codes in this part are borrowed from the SimplexPS toolbox.

%% Power Flow
fprintf('Doing power flow anlaysis...\n')

% ### Power flow analysis
[PowerFlow,~,V_,I_,Ang,P,Q,Vm] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
% Form of PowerFlow{i}: P, Q, V, xi, w

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% % For printting later
% ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);
% ListPowerFlow_ = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

[V,I] = UpdateVI(PowerFlowNew);

% Notes:
% The codes in this part are borrowed from the SimplexPS toolbox. The V and
% I are updated based on the new power flow.

%% Apparatus data and nodal admittance matrix

% Symbolic
s = sym('s');

% Set
W0 = Wbase;

% Calculate nodal admittance matrix
fprintf('Calculating nodal admittance matrix...\n')
Ybus = YbusCalc_s_sym(ListLine,W0,'albe');

% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Convert Ybus to double
dW = 1e-5*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(Wbase+dW)));
Ybus = double(subs(Ybus,'s',1i*Wbase));

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. The good news is that, the
% passive component, and the inner loops of inverters, are indeed symmtric
% in complex dq and alpha/beta frame.

% Get the device source type:
% Notes:
% In the netlist, the device source type has to be listed like this:
% all voltage nodes, all current nodes, all floating bus nodes.
for i = 1:N_Bus
    if DeviceType{i}==0
      	DeviceSourceType{i} = 1;    % Voltage node
    elseif DeviceType{i} == 11
    	DeviceSourceType{i} = 2;    % Current node
    elseif DeviceType{i}==100       
        DeviceSourceType{i} = 3;    % Floating node     
        % DeviceSourceType{i} = 2;                                % ?????
    else
     	error(['Error']);
    end
end

% Find the index of the first floating bus/node
Fbus_1st = N_Bus+1;  % Default: no floating bus
for i = 1:N_Bus
    if DeviceSourceType{i} == 3
        Fbus_1st = i;
        break;
    end
end

% Check if the following buses are all floating buses
for i = Fbus_1st:N_Bus
    if DeviceSourceType{i} ~= 3
        error(['Error']);
    end
end

% Find the index of the first current bus/node
Ibus_1st = N_Bus+1;   % Default: no current bus
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        Ibus_1st = i;
        break;
    end
end

% Check if the following buses are all current buses
for i = Ibus_1st:(Fbus_1st-1)
    if DeviceSourceType{i} ~= 2
        error(['Error']);
    end
end

% =============================
% Add current node
% =============================
if Ibus_1st<=N_Bus

fprintf('Adding current node...\n')

% Get the inner loop parameters
% Notes: Assume all inverters are same
kp_i = DevicePara{Ibus_1st}.kp_i_dq;  
ki_i = DevicePara{Ibus_1st}.ki_i_dq;
Lf   = DevicePara{Ibus_1st}.L;
Rf   = DevicePara{Ibus_1st}.R;

kp_pll = DevicePara{Ibus_1st}.kp_pll;
ki_pll = DevicePara{Ibus_1st}.ki_pll;
PI_pll = kp_pll + ki_pll/s;

Gdel = 1;               % Assume no control delay
Z_PIi = (kp_i + ki_i/(s-1i*W0))*Gdel;
Z_inv = Z_PIi + s*Lf+Rf;
% Y_inv = 1/Z_inv;
Y_inv = (s-1i*W0)/((kp_i + s*Lf+Rf)*(s-1i*W0) + ki_i);

% Convert Y_inv to double
Y_inv_ = double(subs(Y_inv,'s',1i*(W0+dW)));
Y_inv = double(subs(Y_inv,'s',1i*W0));

% Add Y_inv to nodal admittance matrix
if 0
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        % Self branch
        Ybus(i,i) = Ybus(i,i)+Y_inv;
        Ybus_(i,i) = Ybus_(i,i) + Y_inv_;
    end
end
end

end

% =============================
% Eliminate floating node
% =============================
% Notes:
% The floating bus (i.e., no device bus) is assumed as zero-current bus,
% and eliminated here after converting the Y matrix to Y-Z hybrid matrix.
if 1

if Fbus_1st<=N_Bus
    
fprintf('Eliminating floating node...\n')

Ybus = HybridYZMatrix(Ybus,Fbus_1st);
Ybus_ = HybridYZMatrix(Ybus_,Fbus_1st);

if Ibus_1st>N_Bus
   Ibus_1st = Fbus_1st;     % Update Ibus_1st
end
N_Bus = Fbus_1st-1;

Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

end

end

% =============================
% Add voltage node
% =============================
fprintf('Adding voltage node...\n')

J = DevicePara{1}.J;
D = DevicePara{1}.D;
Lsg = DevicePara{1}.L;
Rsg = DevicePara{1}.R;
Zsg = s*Lsg + Rsg;
Y_sg = 1/Zsg;

% Convert Ysg to double
Y_sg_ = double(subs(Y_sg,'s',1i*(Wbase+dW)));
Y_sg = double(subs(Y_sg,'s',1i*Wbase));

J = J*Wbase;
D = D*Wbase;
% Notes: Adding '*Wbase' is because the power swing equation rather than
% the torque swing equation is used, and P=T*w0 if w is not in per unit
% system.
    
if 1

% Prepare star-delta conversion by adding new buses
Ybus = PrepareDY(Ybus,Ibus_1st,N_Bus,Y_sg);
Ybus_ = PrepareDY(Ybus_,Ibus_1st,N_Bus,Y_sg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = HybridYZMatrix(Ybus,N_Bus+1);
Ybus_ = HybridYZMatrix(Ybus_,N_Bus+1);

% Eliminate the zero current bus
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

end

%% Network matrix
fprintf('Calculating network matrix...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = HybridYZMatrix(Ybus,Ibus_1st);
Gbus_ = HybridYZMatrix(Ybus_,Ibus_1st);

% Notes:
% It should be ensured that the buses are listed in the form like this:
% [Vbus1, Vbus2, Vbus3, Ibus4, Ibus5]
    
Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.
Gbus_ = -Gbus_;
    
% Get G_prime
Gbus_prime = (Gbus_ - Gbus)/(1i*dW);

% Nnotes:
% Gbus should statisfy: Output = -Gbus*Input

% Update input and output to network admittance/impedance matrix
Input = [V(1:Ibus_1st-1);
         I(Ibus_1st:end)];
Output = [I(1:Ibus_1st-1);
          V(Ibus_1st:end)];
   
% Get S matrix
S = conj(Input)*transpose(Input);

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
% K = vpa(K);
K = double(K);
K = -K;             % For negative feedback
[K11,K12,K21,K22] = BlockMatrix(K,Ibus_1st-1,Ibus_1st-1);

% Get Gamma matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_Gamma(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        Gamma(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_Gamma(m,n));
    end
end
% Gamma = vpa(Gamma);
Gamma = double(Gamma);
Gamma = -Gamma;    % For negative feedback

%% Apparatus Matrix
fprintf('Calculating apparatus matrix...\n')

% Transfer function form:
% omega = 1/(D + J*s) * W;
T_V = 1/(D + J*s);
F_V = -s/T_V;

% State space form:
% dx/dt = [domega]/dt = [-D/J,0]*[omega] + [1/J]*[W];
%         [dtheta]      [1   ,0] [theta]   [0  ]
% y = [omega] = [1,0]*[omega] + [0]*[W]
%     [theta]   [0,1] [theta]   [0]
Av = [-D/J,0;1,0];
Bv = [1/J;0];
Cv = [1,0;0,1];
Dv = [0;0];
T_V_ss = ss(Av,Bv,Cv,Dv);

% Auto conversion of state space
% T_V_ss = [SimplexPS.sym2ss(T_V);
%           SimplexPS.sym2ss(T_V/s)];

if Ibus_1st<=N_Bus
    
% Transfer function form
% omega = (kp_pll + ki_pll/s) * W
T_I = PI_pll;
F_I = -s/T_I;

% State space form
% dx/dt = [dlamda]/dt = [0 ,0]*[lambda] + [1 ]*[W]
%         [dtheta]      [ki,0] [theta ]   [kp]
% y = [omega] = [ki,0]*[lambda] + [kp]*[W]
%     [theta]   [1 ,0] [theta ]   [0 ]
Ai = [0,0;ki_pll,0];
Bi = [1;kp_pll];
Ci = [ki_pll,0;1,0];
Di = [kp_pll;0];
T_I_ss = ss(Ai,Bi,Ci,Di);

% Auto conversion of state space
% T_I_ss = [SimplexPS.sym2ss(T_I);
%           SimplexPS.sym2ss(T_I/s)];

end

Hinv = eye(N_Bus);      % Let Hinv be identity matrix

%% Stability criterion
fprintf('Calculating stability criterion...\n')
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = BlockMatrix(KH,Ibus_1st-1,Ibus_1st-1);
[phi,xi] = eig(KH);
Gamma_Hphi = inv(phi)*Hinv*Gamma*phi;
[GH11,GH12,GH21,GH22] = BlockMatrix(Gamma_Hphi,Ibus_1st-1,Ibus_1st-1);
[~,sigma,~] = svd(Gamma_Hphi);
sigma_max = max(max(sigma));
% if min( min( real(xi) ) )<-1e-4
%     error(['Error: Min(Real(xi)) = ' num2str(min(min(real(xi)))) ' < 0.']);
% end

% Notes:
% We analyze the current node and voltage node seperately because their
% transfer functions T are different. By asuming the current node is much
% "faster" than the voltage node, we can analyze the current node first
% with assume no voltage node. Then, we can analyze the voltage node by
% adding current node back into K and gamma of the voltage node.

% Current node
if Ibus_1st<=N_Bus
    KH_I = K22;
    GH_I = GH22;
    [phi_I,xi_I] = eig(KH_I);
    [~,sigma_I,~] = svd(GH_I);
    sigma_I_max = max(max(sigma_I));
    [zeta_m_I,w_min_I] = CalcZeta(T_I,diag(xi_I));
  	if min(min(real(xi_I)))<-1e-4
        error(['Error: xi_I_min = ' num2str(min(min(real(xi_I)))) ' < 0.']);
    end
else
    xi_I = [];
    sigma_I = [];
    zeta_m_I = [];
    w_min_I = [];
end

% Voltage node
if isempty(KH22)
    KH_V = KH;
    GH_V = Gamma_Hphi;
else
    KH_V = KH11-KH12*inv(KH22)*KH21;
    GH_V = GH11-KH12*inv(KH22)*GH21;
end
[phi_V,xi_V] = eig(KH_V);
[~,sigma_V,~] = svd(GH_V);
sigma_V_max = max(max(sigma_V));
[zeta_m_V,w_min_V] = CalcZeta(T_V,diag(xi_V));
if min(min(real(xi_V)))<-1e-4
    error(['Error: xi_V_min = ' num2str(min(min(real(xi_V)))) ' < 0.']);
end

%% Calculating the state space representation
% Get whole system Tss
Tss = [[],[],[],[]];
for i = 1:N_Bus
    if DeviceSourceType{i} == 1
        Tss = append(Tss,T_V_ss);
    elseif DeviceSourceType{i} == 2
        Tss = append(Tss,T_I_ss);
    else
        error(['Error']);
    end
end
feedin = [1:N_Bus];
feedout_L1 = [1:N_Bus]*2;
feedout_L2 = [1:N_Bus]*2-1;
Tss_ = feedback(Tss,KH,feedin,feedout_L1);
Tss__ = feedback(Tss_,Gamma_Hphi,feedin,feedout_L2);

pole_sys = pole(Tss__)/2/pi;

%% Plot: symbolic

figure_n = 0;

w_p = logspace(-1,1,500)*2*pi;
w_pn = [-flip(w_p),w_p];
s_pn = 1i*w_pn;

figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_V,s_pn);
scatter(real(diag(xi_V)),imag(diag(xi_V)),'x','LineWidth',1.5); hold on; grid on;

if Ibus_1st<=N_Bus
figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_I,s_pn);
scatter(real(diag(xi_I)),imag(diag(xi_I)),'x','LineWidth',1.5); hold on; grid on;
end

%% Plot: state space

figure_n = 1000;

figure_n = figure_n+1;
figure(figure_n)
scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;

%% Output

xi
sigma

xi_I
sigma_I
zeta_m_I

xi_V
sigma_V
zeta_m_V
Dpu_crit = max(max(sigma_V)*W0)     % Critical damping for SGs

%% Check the stability
% Criterion
fprintf('Check the stability by the proposed criterion:\n')
if min(zeta_m_V)>sigma_V_max
    StableVoltageNode = 1;
else
    StableVoltageNode = 0;
end
StableCurrentNode = 1;
if ~isempty(zeta_m_I)
    if min(zeta_m_I)>sigma_I_max
        StableCurrentNode = 1;
    else
        StableCurrentNode = 0;
    end
end
if StableVoltageNode==1 && StableCurrentNode==1
    fprintf('Stable!\n')
else
 	if StableVoltageNode==0
        fprintf('Unstable voltage node!\n')
    end
    if StableCurrentNode==0
        fprintf('Unstable current node!\n')
    end
end

% Pole
fprintf('Check the stability by poles:\n')
if isempty(find(real(pole_sys)>1e-9, 1))
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end