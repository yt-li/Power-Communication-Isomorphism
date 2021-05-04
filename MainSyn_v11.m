% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% PLi and QLi should be combined into the nodal admittance matrix.
%
% Devices should be listed in the same order as the bus.
%
% The flux inductor of SGs should be combined into the line impedances.
%
% The reference direction should be doubled checked: Ybus, Gbus, T,

%% Clear
clear all
clc
close all

%% Select data

% Single SG
% UserData = 'SingleSG.xlsx';         % SG + resistive load

% Single IBR
% UserData = 'SingleIBR.xlsx';    	% IBR + resistive load
% UserData = 'SingleIBR_v1';

% Multi-machine system
% UserData = 'ExamplePowerSystem_v3.xlsx';    	% SG + SG + IBR

% UserData = 'ExamplePowerSystem_v4';         % SG + IBR

% Nature example
UserData = 'Nature_NETS_NYPS_68Bus_original';
% UserData = 'Nature_NETS_NYPS_68Bus_IBR_all';
% UserData = 'Nature_NETS_NYPS_68Bus_1SG_OtherIBR';
% UserData = 'Nature_NETS_NYPS_68Bus_2SG_OtherIBR';
% UserData = 'Nature_NETS_NYPS_68Bus_IBR';
% UserData = 'Nature_NETS_NYPS_68Bus_IBR_v1';

% Test Power Flow
% UserData = '3MachineModel_v1_test';
% UserData = '3MachineModel_v1_3SG';
% UserData = '3MachineModel_v2_SG_IBR_SG';
% UserData = '3MachineModel_v3_SG_IBR_IBR';
% UserData = '3MachineModel_v4_3IBR';
% UserData = '2MachineModel_test';

% UserData = 'IEEE_14Bus';
% UserData = 'IEEE_30Bus';
% UserData = 'IEEE_57Bus';

%% Enable
Enable_Plot_ComparisonToolbox   = 1;
Enable_VoltageNode_InnerLoop    = 1;
Enable_CurrentNode_InnerLoop    = 0;

%% Load data
fprintf('Loading data...\n')
LoadData();

%% Power Flow
fprintf('Doing power flow anlaysis...\n')
PowerFlowCal();

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Set
W0 = Wbase;

% Calculate nodal admittance matrix
fprintf('Calculating nodal admittance matrix...\n')
Ybus = YbusCalc_s_sym(ListLineNew,W0,'albe');

%% Apparatus Data
% Get the device source type:
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

% Notes:
% In the netlist, the device source type has to be listed like this:
% all voltage nodes, all current nodes, all floating bus nodes.
%
% This should be improved. Meanwhile, Ybus and PowerFlow should also be
% updated.

% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Convert Ybus to double
dW = 1e-10*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(W0+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*W0));

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. The good news is that, the
% passive component, and the inner loops of inverters, are indeed symmtric
% in complex dq and alpha/beta frame.

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
% Eliminate floating node
% =============================
% Notes:
% The floating bus (i.e., no device bus) is assumed as zero-current bus,
% and eliminated here after converting the Y matrix to Y-Z hybrid matrix.
if n_Fbus_1st>N_Bus
    fprintf('Warning: The system has no floating node.\n')
else
    
fprintf('Eliminating floating node...\n')

Ybus = HybridMatrixYZ(Ybus,n_Fbus_1st);
Ybus_ = HybridMatrixYZ(Ybus_,n_Fbus_1st);

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
if n_Ibus_1st == 1
    fprintf('Warning: The system has no voltage node.\n')
else
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

% Doing D-Y conversion
if Enable_VoltageNode_InnerLoop                                           	% ??? 

% Prepare star-delta conversion by adding new buses
Ybus = PrepareConvertDY(Ybus,n_Ibus_1st,N_Bus,Y_sg);
Ybus_ = PrepareConvertDY(Ybus_,n_Ibus_1st,N_Bus,Y_sg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = HybridMatrixYZ(Ybus,N_Bus+1);
Ybus_ = HybridMatrixYZ(Ybus_,N_Bus+1);

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

end

Ybus_dwn = Ybus_;
Ybus_dwm = Ybus_;

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
kp_pll{i} = kp_pll{i}/abs(P(i));
ki_pll{i} = ki_pll{i}/abs(P(i));
PI_pll{i} = kp_pll{i} + ki_pll{i}/s;

% Notes:
% Because we use Q-PLL later, we scale PI_pll here by P (approximately
% id) to ensure the actual PLL bandwidth is same to that of a v_q-PLL.

end

wm = sym('wm');

% alpha/beta
Z_PIi = kp_i + ki_i/(s-1i*wm);        
Z_Lf = s*Lf+Rf;
% ki_i = ki_i*5;                                                                 % ??????
Y_inv = (s-1i*wm)/((kp_i + s*Lf+Rf)*(s-1i*wm) + ki_i);

Y_inv_dwn = subs(Y_inv,'wm',W0);
Y_inv_dwn = double(subs(Y_inv_dwn,'s',1i*(W0+dW)));

Y_inv_dwm = subs(Y_inv,'s',1i*W0);
Y_inv_dwm = double(subs(Y_inv_dwm,'wm',(W0+dW)));

Y_inv = subs(Y_inv,'wm',W0);
Y_inv_ = double(subs(Y_inv,'s',1i*(W0+dW)));
Y_inv = double(subs(Y_inv,'s',1i*W0));

% Calculate Y_inv_prime
Y_inv_prime = (Y_inv_ - Y_inv)/(1i*dW);

% Notes:
% When current controller is very fast, Y_inv -> 0 and can be ignored.

% Add Y_inv to nodal admittance matrix
if Enable_CurrentNode_InnerLoop                                                 % ??? 
    
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        % Self branch
        Ybus(i,i) = Ybus(i,i) + Y_inv;
        Ybus_(i,i) = Ybus_(i,i) + Y_inv_;
        Ybus_dwn(i,i) = Ybus_dwn(i,i) + Y_inv_dwn;
        Ybus_dwm(i,i) = Ybus_dwm(i,i) + Y_inv_dwm;
    end
end

% Update I
V = V(1:N_Bus,:);
I = Ybus*V;

end

end

% Notes:
% Doing D-Y conversion first, add Yinv next???

%% Network matrix
fprintf('Calculating network matrix...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = HybridMatrixYZ(Ybus,n_Ibus_1st);
Gbus_ = HybridMatrixYZ(Ybus_,n_Ibus_1st);
Gbus_dwn = HybridMatrixYZ(Ybus_dwn,n_Ibus_1st);
Gbus_dwm = HybridMatrixYZ(Ybus_dwm,n_Ibus_1st);

% Notes:
% It should be ensured that the buses are listed in the form like this:
% [Vbus1, Vbus2, Vbus3, Ibus4, Ibus5, ...]
    
Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.
Gbus_ = -Gbus_;
Gbus_dwn = -Gbus_dwn;
Gbus_dwm = -Gbus_dwm;
    
% Get G_prime
% Notes: It is calculaed by numerical calculation
if 1                                                                            % ???
Gbus_prime = (Gbus_ - Gbus)/(1i*dW);  
else
Gbus_prime = (Gbus_dwn - Gbus)/(1i*dW) + (Gbus_dwm - Gbus)/(1i*dW);
end

% Notes:
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
      	if real(I(i))>0
            epsilon(i) = -epsilon(i);
            test1 = 1
        end
    else
        error(['Error']);
    end
end
% Notes:
% For current node, conventional PLL uses vq rather than Q to achieve the
% synchronization. In load convention, Q = vq*id - vd*iq. If iq = 0, Q =
% vq*id if iq = 0. That means, the power flow direction id will influence
% the sign of Q, which also means the sign of loop gain also has to be
% changed to ensure the system stability. Here, we change the the value of
% epsilon_m depending on the power flow, to ensure xi >=0 and the
% stability. Noting that Q is in load convention, so, id should also be
% load convention. That means, when I>0 which means the IBR is generating
% active power, the sign of epsilon should be changed.

% For test                                                                      % ??????
% epsilon(2) = -epsilon(2)
% test2 = 1

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
% Notes:
% For single-machine-load system, i.e., a single bus system, K = 0, based
% on equation (15) in the paper.

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

% Notes:
% K can be interpreted as the synchronizing torque coefficient, as dS is
% proportional to K*dtheta. Gamma can be interpreted as the damping torque
% coefficient, as dS is proportional to Gamma*dtheta.

%% Apparatus Matrix
fprintf('Calculating apparatus matrix...\n')

% Initialize Hinv
Hinv = eye(N_Bus);      % Let Hinv be identity matrix initially

% Choose the reference node
n_i_ref = n_Ibus_1st;        % Select the first current node as the reference
n_v_ref = 1;               % Select the first voltage node as the reference

% ### Voltage node
if n_Ibus_1st ~= 1
    
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

end

% ### Current node
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
if 1                                                                          % ???
% Without LPF
% dx/dt = [dw_pll_i]/dt = [0 ,0]*[w_pll_i] + [ki]*[W]
%         [dtheta  ]      [1, 0] [theta  ]   [kp]
% y = [omega] = [1, 0]*[w_pll_i] + [kp]*[W]
%     [theta]   [0 ,1] [theta  ]   [0 ]
Ai = [0, 0;
      1, 0];
Bi = [ki_pll{n_i_ref};
      kp_pll{n_i_ref}];
Ci = [1, 0;
      0, 1];
Di = [kp_pll{n_i_ref};
      0];
else
% With LPF
% dx/dt = [dw      ]/dt
%         [dw_pll_i]
%         [dtheta  ]
% y = [omega]
%     [theta]
tau = 1/(2*pi*500);
Ai = [-1/tau,1/tau,0;
      0,0,0;
      1,0,0];
Bi = [kp_pll{n_i_ref};
      ki_pll{n_i_ref};
      0];
Ci = [1,0,0;
      0,0,1];
Di = [0;
      0];
end

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
%
% If analyzing voltage and current nodes seperately, KH and Gamma_Hphi for
% two subloops should be re-calculated based on K, Hinv, Gamma, and can not
% be seperate directly from the whole system KH and Gamma_Hphi.

% ### Voltage node
% This calculation might also be wrong! Because H is not eye matrix now!        	???  
if n_Ibus_1st~= 1
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
else
	xi_V = [];
    sigma_V = [];
    zeta_m_V = [];
    w_min_V = [];
end

% ### Current node
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
T1cl = minreal(T1cl);
T12cl = minreal(T12cl);

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

if Enable_Plot_ComparisonToolbox
figure_n = figure_n+1;
figure(figure_n)
scatter(real(pole_sys_T12cl),imag(pole_sys_T12cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_sys_T1cl),imag(pole_sys_T1cl),'x','LineWidth',1.5); hold on; grid on;
pole_sys_toolbox = load('pole_sys').pole_sys;
index = find(abs(imag(pole_sys_toolbox))<35);
pole_sys_toolbox = pole_sys_toolbox(index);
scatter(real(pole_sys_toolbox),imag(pole_sys_toolbox),'x','LineWidth',1.5); hold on; grid on;
legend('Yunjie With F Shift','Yunjie Wihtout F Shift','Toolbox')
end

%% Check stability
% Criterion
fprintf('Check the stability by the proposed criterion:\n')

% Set the default values
StableVoltageNode = 1;
StableCurrentNode = 1;

% Check voltage nodes
if ~isempty(zeta_m_V)
    if min(zeta_m_V)>=sigma_V_max
        StableVoltageNode = 1;
    else
        StableVoltageNode = 0;
    end
end

% Check current nodes
if ~isempty(zeta_m_I)
    if min(zeta_m_I)>=sigma_I_max
        StableCurrentNode = 1;
    else
        StableCurrentNode = 0;
    end
end

% Output
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