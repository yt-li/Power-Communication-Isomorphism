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
% Pure SG system
% UserData = 'ExamplePowerSystem_v2.xlsx';      % pure SG system, 5 machines
% UserData = 'ExamplePowerSystem_v2_1.xlsx';   	% pure SG system, 2 machine
% UserData = 'ExamplePowerSystem_v2_2.xlsx';   	% pure SG system, 3 machine
% UserData = 'ExamplePowerSystem_v2_3.xlsx';   	% pure SG system, random J,D
% UserData = 'ExamplePowerSystem_v2_5.xlsx';   	% pure SG system, 3 machine, large Lflux and small Lline
% UserData = 'ExamplePowerSystem_v2_6.xlsx';    % pure SG system, 5 machines, random all parameters with Zsg=0
% UserData = 'ExamplePowerSystem_v2_7.xlsx';    % pure SG system, 5 machines, random all parameters
% UserData = 'ExamplePowerSystem_v2_8.xlsx';   	% pure SG system, 3 machine, take Lflux to line
% UserData = 'ExamplePowerSystem_v2_9.xlsx';   	% equivalent to v2_8 but with Lflux inside SG
% UserData = 'ExamplePowerSystem_v2_10.xlsx';     % pure SG system, 5 machines, very small line impedance
% UserData = 'ExamplePowerSystem_v2_11.xlsx';    % pure SG system, 5 machines, random all parameters, with floating buses
% UserData = 'ExamplePowerSystem_v2_12.xlsx';    % pure SG system, 5 machines, random all parameters, with floating buses, with floating PL
% UserData = 'ExamplePowerSystem_v2_13.xlsx';    % pure SG system, 5 machines, random all parameters, with high-order floating bus

% Hybrid SG-IBR system
UserData = 'ExamplePowerSystem_v3.xlsx';    	% SG + IBR
% UserData = 'ExamplePowerSystem_v3_1.xlsx';    % SG + IBR + floating bus
% UserData = 'ExamplePowerSystem_v3_2.xlsx';  	% SG + floating bus
% UserData = 'ExamplePowerSystem_v3_3.xlsx';   	% One SG + One IBR + No load
% UserData = 'ExamplePowerSystem_v3_4.xlsx';  	% One SG + One IBR

% Nature example
% UserData = 'Nature_NETS_NYPS_68Bus_original.xlsx';        % Standard
% UserData = 'Nature_NETS_NYPS_68Bus_NoR';    	% No resistance in Lline and Lflux     
% UserData = 'Nature_NETS_NYPS_68Bus_IBR';
% UserData = 'Nature_NETS_NYPS_68Bus_NoTrans';
% UserData = 'Nature_NETS_NYPS_68Bus_NoQL';
% UserData = 'Nature_NETS_NYPS_68Bus_NegQLi';

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

% For printting later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);
ListPowerFlow_ = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

% Update V and I
[V,I] = PowerFlowUpdateVI(PowerFlowNew);

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
Ybus = YbusCalc_s_sym(ListLineNew,W0,'albe');

% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Convert Ybus to double
dW = 1e-5*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(W0+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*W0));

% Ybus_ = double(subs(Ybus,'s',1i*(0+dW)));   % Used for calculating derivative numerically
% Ybus = double(subs(Ybus,'s',1i*0));

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
    if DeviceType{i}==1
      	DeviceSourceType{i} = 1;    % Voltage node
    elseif DeviceType{i}==11
    	DeviceSourceType{i} = 2;    % Current node
    elseif DeviceType{i}==100       
        DeviceSourceType{i} = 3;    % Floating node     
    else
     	error(['Error']);
    end
end

% Find the index of the first floating bus/node
n_Fbus_1st = N_Bus+1;  % Default: no floating bus
for i = 1:N_Bus
    if DeviceSourceType{i} == 3
        n_Fbus_1st = i;
        break;
    end
end

% Check if the following buses are all floating buses
for i = n_Fbus_1st:N_Bus
    if DeviceSourceType{i} ~= 3
        error(['Error']);
    end
end

% Find the index of the first current bus/node
n_Ibus_1st = N_Bus+1;   % Default: no current bus
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        n_Ibus_1st = i;
        break;
    end
end

% Check if the following buses are all current buses
for i = n_Ibus_1st:(n_Fbus_1st-1)
    if DeviceSourceType{i} ~= 2
        error(['Error']);
    end
end

% =============================
% Add current node
% =============================
if n_Ibus_1st>N_Bus
    fprintf('Warning: The system has no current node.\n')
else

fprintf('Adding current node...\n')

% Get the inner loop parameters
% Notes: Assume all inverters are same
for i = n_Ibus_1st:(n_Fbus_1st-1)
kp_i = DevicePara{i}.kp_i_dq;  
ki_i = DevicePara{i}.ki_i_dq;
Lf   = DevicePara{i}.L;
Rf   = DevicePara{i}.R;

kp_pll{i} = DevicePara{i}.kp_pll;
ki_pll{i} = DevicePara{i}.ki_pll;
PI_pll{i} = kp_pll{i} + ki_pll{i}/s;
end

Gdel = 1;               % Assume no control delay
Z_PIi = (kp_i + ki_i/(s-1i*W0))*Gdel;
Z_inv = Z_PIi + s*Lf+Rf;
% Y_inv = 1/Z_inv;
Y_inv = (s-1i*W0)/((kp_i + s*Lf+Rf)*(s-1i*W0) + ki_i);

% Convert Y_inv to double
Y_inv_ = double(subs(Y_inv,'s',1i*(W0+dW)));
Y_inv = double(subs(Y_inv,'s',1i*W0));
% Notes:
% When current controller is very fast, Y_inv -> 0 and can be ignored.

% Add Y_inv to nodal admittance matrix
if 0                                                                            % ??? 
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
if n_Fbus_1st>N_Bus
    fprintf('Warning: The system has no floating node.\n')
else
    
fprintf('Eliminating floating node...\n')

Ybus = HybridYZMatrix(Ybus,n_Fbus_1st);
Ybus_ = HybridYZMatrix(Ybus_,n_Fbus_1st);

if n_Ibus_1st>N_Bus
   n_Ibus_1st = n_Fbus_1st;     % Update Ibus_1st
end
N_Bus = n_Fbus_1st-1;

Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

end

% =============================
% Add voltage node
% =============================
fprintf('Adding voltage node...\n')

for i = 1:(n_Ibus_1st-1)
J{i} = DevicePara{i}.J;
D{i} = DevicePara{i}.D;
J{i} = J{i}*Wbase;
D{i} = D{i}*Wbase;
% Notes: Adding '*Wbase' is because the power swing equation rather than
% the torque swing equation is used, and P=T*w0 if w is not in per unit
% system.

Lsg = DevicePara{i}.L;
Rsg = DevicePara{i}.R;
Zsg = s*Lsg + Rsg;
Y_sg{i} = 1/Zsg;

% Convert Ysg to double
Y_sg_{i} = double(subs(Y_sg{i},'s',1i*(W0+dW)));
Y_sg{i} = double(subs(Y_sg{i},'s',1i*W0));
end    

if 1                                                                            % ??? 

% Prepare star-delta conversion by adding new buses
Ybus = PrepareDY(Ybus,n_Ibus_1st,N_Bus,Y_sg);
Ybus_ = PrepareDY(Ybus_,n_Ibus_1st,N_Bus,Y_sg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = HybridYZMatrix(Ybus,N_Bus+1);
Ybus_ = HybridYZMatrix(Ybus_,N_Bus+1);

% Eliminate the old voltage bus, i.e., zero current bus
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

% Update V and I
% Notes: The steady-state voltage at voltage buses are changed if we split
% the inductor outside the apparatus. The following calculation is
% equivalent to "V = V+I/Y_sg".
I = I(1:N_Bus,:);
V = inv(Ybus)*I;

end

%% Network matrix
fprintf('Calculating network matrix...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = HybridYZMatrix(Ybus,n_Ibus_1st);
Gbus_ = HybridYZMatrix(Ybus_,n_Ibus_1st);

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
Input = [V(1:n_Ibus_1st-1);
         I(n_Ibus_1st:end)];
Output = [I(1:n_Ibus_1st-1);
          V(n_Ibus_1st:end)];
   
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
[K11,K12,K21,K22] = BlockMatrix(K,n_Ibus_1st-1,n_Ibus_1st-1);

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
[Gamma11,Gamma12,Gamma21,Gamma22] = BlockMatrix(Gamma,n_Ibus_1st-1,n_Ibus_1st-1);

%% Apparatus Matrix
fprintf('Calculating apparatus matrix...\n')

% Initialize Hinv
Hinv = eye(N_Bus);      % Let Hinv be identity matrix initially

% Choose the reference node
n_i_ref = n_Ibus_1st;        % Select the first current node as the reference
n_v_ref = 1;               % Select the first voltage node as the reference

% Symbolic transfer function form:
% omega = 1/(D + J*s) * W;
for i = 1:(n_Ibus_1st-1)
T_V{i} = 1/(D{i} + J{i}*s);
end
F_V = -s/T_V{n_v_ref};

% Update Hinv
for i = 1:(n_Ibus_1st-1)
    Hinv(i,i) = 1/(J{i}/J{n_v_ref});
    Hinv(i,i) = double(Hinv(i,i));
end

% State space form:
% dx/dt = [domega]/dt = [-D/J,0]*[omega] + [1/J]*[W];
%         [dtheta]      [1   ,0] [theta]   [0  ]
% y = [omega] = [1,0]*[omega] + [0]*[W]
%     [theta]   [0,1] [theta]   [0]
Av = [-D{n_v_ref}/J{n_v_ref}, 0;
      1,                      0];
Bv = [1/J{n_v_ref};
      0];
Cv = [1,0;
      0,1];
Dv = [0;
      0];
T_V_ss = ss(Av,Bv,Cv,Dv);

% Auto conversion of state space
% T_V_ss = [SimplexPS.sym2ss(T_V);
%           SimplexPS.sym2ss(T_V/s)];

if n_Ibus_1st<=N_Bus
    
% Symbolic transfer function form
% omega = (kp_pll{i} + ki_pll{i}/s) * W
for i = n_Ibus_1st:N_Bus
    T_I{i} = PI_pll{i};
end
F_I = -s/T_I{n_i_ref};

% Update Hinv for current node
for i = n_Ibus_1st:N_Bus
    Hinv(i,i) = kp_pll{i}/kp_pll{n_i_ref};
    Hinv(i,i) = double(Hinv(i,i));
end

% State space form
% dx/dt = [dlamda]/dt = [0 ,0]*[lambda] + [1 ]*[W]
%         [dtheta]      [ki,0] [theta ]   [kp]
% y = [omega] = [ki,0]*[lambda] + [kp]*[W]
%     [theta]   [1 ,0] [theta ]   [0 ]
Ai = [0,               0;
      ki_pll{n_i_ref}, 0];
Bi = [1;
      kp_pll{n_i_ref}];
Ci = [ki_pll{n_i_ref}, 0;
      1,               0];
Di = [kp_pll{n_i_ref};
      0];
T_I_ss = ss(Ai,Bi,Ci,Di);

% Auto conversion of state space
% T_I_ss = [SimplexPS.sym2ss(T_I);
%           SimplexPS.sym2ss(T_I/s)];

end

%% Stability criterion
fprintf('Calculating stability criterion...\n')
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = BlockMatrix(KH,n_Ibus_1st-1,n_Ibus_1st-1);
[phi,xi] = eig(KH);
Gamma_Hphi = inv(phi)*Hinv*Gamma*phi;
[GH11,GH12,GH21,GH22] = BlockMatrix(Gamma_Hphi,n_Ibus_1st-1,n_Ibus_1st-1);
[~,sigma,~] = svd(Gamma_Hphi);
sigma_max = max(max(sigma));
if min( min( real(xi) ) )<-1e-4
    fprintf(['Error: Min(Real(xi)) = ' num2str(min(min(real(xi)))) ' < 0.\n']);
end

% Notes:
%
% We analyze the current node and voltage node seperately because their
% transfer functions T are different. By asuming the current node is much
% "faster" than the voltage node, we can analyze the current node first
% with assume no voltage node. Then, we can analyze the voltage node by
% adding current node back into K and gamma of the voltage node.
%
% If the power system is in the loop topology. KH and Gamma_Hphi should be
% semi-diagonalized and should have some relation. The values of KH will
% rotate in each row.

% Current node
% This calculation might be wrong!                                                      ???  
if n_Ibus_1st<=N_Bus
    % KH_I = K22;
    % Gamma_Hphi_I = GH22;
    K_I = K22;
    Gamma_I = Gamma22;
    Hinv_I = Hinv(n_Ibus_1st:N_Bus,n_Ibus_1st:N_Bus);
    KH_I = Hinv_I*K_I;
    [phi_I,xi_I] = eig(KH_I);
    Gamma_Hphi_I = inv(phi_I)*Hinv_I*Gamma_I*phi_I;
    [~,sigma_I,~] = svd(Gamma_Hphi_I);
    sigma_I_max = max(max(sigma_I));
    [zeta_m_I,w_min_I] = CalcZeta(T_I{n_i_ref},diag(xi_I));
  	if min(min(real(xi_I)))<-1e-4
        fprintf(['Error: xi_I_min = ' num2str(min(min(real(xi_I)))) ' < 0.\n']);
    end
else
    xi_I = [];
    sigma_I = [];
    zeta_m_I = [];
    w_min_I = [];
end

% Voltage node
% This calculation might also be wrong! Because H is not eye matrix now!        	???  
if isempty(KH22)
    KH_V = KH;
    [phi_V,xi_V] = eig(KH_V);
    Gamma_Hphi_V = Gamma_Hphi;
    [~,sigma_V,~] = svd(Gamma_Hphi_V);
else
    K_V = K11 - K12*inv(K22)*K21;
    Gamma_V =  Gamma11 - K12*inv(K22)*Gamma21;
    Hinv_V = Hinv(1:(n_Ibus_1st-1),1:(n_Ibus_1st-1));
    KH_V = Hinv_V*K_V;
    [phi_V,xi_V] = eig(KH_V);
    Gamma_Hphi_V = inv(phi_V)*Hinv_V*Gamma_V*phi_V;
    [~,sigma_V,~] = svd(Gamma_Hphi_V);
end
sigma_V_max = max(max(sigma_V));
[zeta_m_V,w_min_V] = CalcZeta(T_V{n_v_ref},diag(xi_V));
if min(min(real(xi_V)))<-1e-4
    fprintf(['Error: xi_V_min = ' num2str(min(min(real(xi_V)))) ' < 0.']);
end

%% Calculate the state space representation
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

% Calculate the whole-system closed loop state space model
feedin = [1:N_Bus];
feedout_L1 = [1:N_Bus]*2;
feedout_L2 = [1:N_Bus]*2-1;
T1cl = feedback(Tss*Hinv,K,feedin,feedout_L1);
T12cl = feedback(T1cl,Gamma,feedin,feedout_L2);

% Calculate the whole system pole
pole_sys_T1cl = pole(T1cl)/2/pi;
pole_sys_T12cl = pole(T12cl)/2/pi;

%% Plot: symbolic

if 0

figure_n = 0;

w_p = logspace(-1,1,500)*2*pi;
w_pn = [-flip(w_p),w_p];
s_pn = 1i*w_pn;

figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_V,s_pn);
scatter(real(diag(xi_V)),imag(diag(xi_V)),'x','LineWidth',1.5); hold on; grid on;

if n_Ibus_1st<=N_Bus
figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_I,s_pn);
scatter(real(diag(xi_I)),imag(diag(xi_I)),'x','LineWidth',1.5); hold on; grid on;
end

end

%% Plot: state space

figure_n = 1000;

figure_n = figure_n+1;
figure(figure_n)
scatter(real(pole_sys_T12cl),imag(pole_sys_T12cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_sys_T1cl),imag(pole_sys_T1cl),'x','LineWidth',1.5); hold on; grid on;
legend('L12','L1')

figure_n = figure_n+1;
figure(figure_n)
scatter(real(pole_sys_T12cl),imag(pole_sys_T12cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_sys_T1cl),imag(pole_sys_T1cl),'x','LineWidth',1.5); hold on; grid on;
pole_sys_toolbox = load('pole_sys').pole_sys;
index = find(abs(imag(pole_sys_toolbox))<45);
pole_sys_toolbox = pole_sys_toolbox(index);
index = find(abs(real(pole_sys_toolbox))<1);
pole_sys_toolbox = pole_sys_toolbox(index);
scatter(real(pole_sys_toolbox),imag(pole_sys_toolbox),'x','LineWidth',1.5); hold on; grid on;
legend('Yunjie With F Shift','Yunjie Wihtout F Shift','Toolbox')

%% Check stability
% Criterion
fprintf('Check the stability by the proposed criterion:\n')
if min(zeta_m_V)>=sigma_V_max
    StableVoltageNode = 1;
else
    StableVoltageNode = 0;
end
StableCurrentNode = 1;
if ~isempty(zeta_m_I)
    if min(zeta_m_I)>=sigma_I_max
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
if isempty(find(real(pole_sys_T12cl)>1e-4, 1))
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end