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

%% Clear
clear all
clc
close all

%% Load data
fprintf('Loading data...\n')

% UserData = 'ExamplePowerSystem_v1.xlsx';
% UserData = 'ExamplePowerSystem_v1_1.xlsx';
UserData = 'ExamplePowerSystem_v2.xlsx';

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

%% Power Flow
fprintf('Doing power flow anlaysis...\n')

% ### Power flow analysis
[PowerFlow,~,V,I,Ang,P,Q,Vm] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
% Form of PowerFlow{i}: P, Q, V, xi, w
% for i = 1:N_Bus
%     SynPowerFlow{i} = [S(i),P(i),Q(i),V(i),I(i)];
% end

% For printing later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlow_ = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

% Notes:
% Because V and I rather than PowerFlow are used later, the Load2SelfBranch
% will not work here. We have to make sure PLi and QLi = 0 in the netlist.

%% Network Matrix
fprintf('Calculating network matrix...\n')
% New Ybus
% Ybus = SimplexPS.PowerFlow.YbusCalc(ListLine);
Ybus = YbusCalc_s(ListLine,Wbase,'albe');
% Ybus = YbusCalc_w(ListLine,Wbase);

% Notes:
% If using s-domain Ybus calculation and using vector as input, then, the
% system admittance matrix has to be symmetric. The good news is that, the
% passive component, and the inner loops of inverters, are indeed symmtric
% in complex dq frame.

% Seperate Ybus according to source types
Is_n = [];
for i = 1:N_Bus
    if DeviceSourceType{i} == 2
        Is_n = i;
        break;
    end
end
if isempty(Is_n)
    Is_n = N_Bus+1;
end
Y11 = Ybus(1:Is_n-1,1:Is_n-1);
Y12 = Ybus(1:Is_n-1,Is_n:end);
Y21 = Ybus(Is_n:end,1:Is_n-1);
Y22 = Ybus(Is_n:end,Is_n:end);

% Get Gbus
G11 = Y11 - Y12*inv(Y22)*Y21;
G12 = -Y12*inv(Y22);
G21 = inv(Y22)*Y21;
G22 = inv(Y22);

Gbus = [G11,G12;
        G21,G22];
Gbus = -Gbus;       % To change power to load convention
    
% Get Gbus_prime        % This value is wrong. Should be dGbus/ds
Gbus_prime = diff(Gbus,'s');

% Convert sym to num
Gbus = subs(Gbus,'s',1i*Wbase);
Gbus_prime = subs(Gbus_prime,'s',1i*Wbase);
Gbus = double(Gbus);
Gbus_prime = double(Gbus_prime);
    
% Update input and output
In = [V(1:Is_n-1);
      I(Is_n:end)];
Out = [I(1:Is_n-1);
       V(Is_n:end)];
   
   % The convention should be double-checked ?????????????????
   % In and Out, Gbus, T(s)
   
% Get S matrix
% In(1) = abs(In(1))*exp(1i*0);
% In(2) = abs(In(2))*exp(1i*pi/2/10*1);
% In(3) = abs(In(2))*exp(1i*pi/2/10*2);
% In(4) = abs(In(2))*exp(1i*pi/2/10*3);
% In(5) = abs(In(2))*exp(1i*pi/2/10*4);
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
K = double(K);
K = -K;             % For negative feedback

% Get Gamma matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_Gamma(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        Gamma(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_Gamma(m,n));
    end
end
Gamma = double(Gamma);
Gamma = - Gamma;    % For negative feedback

%% Apparatus Matrix
fprintf('Calculating apparatus matrix...\n')

% Symbolic
s = sym('s');

% Synchronous generaotr
D = 5/Wbase^2;
J = 3.5*2/Wbase^2;

% Droop
m = 1/D;
Tf = J*m;
% increase m is equivalent to decrease J?

% PI
w_pll = 10*2*pi;
kp_pll = w_pll;
ki_pll = w_pll^2/4;
PI_pll = kp_pll + ki_pll/s;

THinv = [];
THinv = sym(THinv);
for i = 1:N_Bus
    if DeviceType{i} == 0
        THinv_i(i) = 1/(D + J*s);
    elseif DeviceType{i} == 11
        if 0
            THinv_i(i) = PI_pll;
        else
            % TH(i) = 1/(1/m + Tf/m*s);
            THinv_i(i) = 1/(D + J*s)*10;
        end
    elseif DeviceType{i} == 20
        THinv_i(i) = m/(1+Tf*s);
    else
        error(['Error']);
    end
    THinv = blkdiag(THinv,THinv_i(i));
end

T = 1/(D+J*s);
Hinv = simplify(THinv/T);
Hinv = double(Hinv);

if ~SimplexPS.is_eye(Hinv)
    error(['Error']);
end
    
F = -s/T;

%% Stability criterion
fprintf('Calculating stability criterion...\n')
K_H = Hinv*K;

[phi,xi] = eig(K_H);
xi_diag = diag(xi); % Get the diagonal elements

Gamma_Hphi = inv(phi)*inv(Hinv)*Gamma*phi;

[~,sigma,~] = svd(Gamma_Hphi);
sigma_max = max(max(sigma));

w = sym('w','real');
T = subs(T,'s',1i*w);

for i = 1:length(xi_diag)
    
    if xi_diag(i) > 0
        
%         zeta{i} = (xi_diag(i) + s/T)/(s);
% 
%         dzeta{i} = diff(zeta{i},s);
%         dzeta{i} = simplify(dzeta{i});
% 
%         ddzeta{i} = diff(dzeta{i},s);
%         ddzeta{i} = simplify(ddzeta{i});
%
%         zeta_{i} = matlabFunction(zeta{i});
%         dzeta_{i} = matlabFunction(dzeta{i});
%         ddzeta_{i} = matlabFunction(ddzeta{i});
%
%         s_equili{i} = fzero(dzeta_{i},100);
%         zeta_m(i) = zeta_{i}(s_equili{i});
%         
%       	if s_equili{i}<=0
%             error(['Error']);
%         end

        zeta0{i} = (xi_diag(i) + s/T)/s;
        zeta0{i} = subs(zeta0{i},'s',1i*w);
        
        zeta{i} = zeta0{i}*conj(zeta0{i});
        zeta{i} = simplify(zeta{i});

%         dzeta{i} = diff(zeta{i},w);
%         dzeta{i} = simplify(dzeta{i});
% 
%         ddzeta{i} = diff(dzeta{i},w);
%         ddzeta{i} = simplify(ddzeta{i});

        % s_equili{i} = roots(dzeta{i});

        zeta0_{i} = matlabFunction(zeta0{i});
        zeta_{i} = matlabFunction(zeta{i});
%         dzeta_{i} = matlabFunction(dzeta{i});
%         ddzeta_{i} = matlabFunction(ddzeta{i});
        
        w_min(i) = fminbnd(zeta_{i},-1e4,1e4)

%         w_min{i} = fzero(dzeta_{i},100);
         zeta_m(i) = abs(zeta0_{i}(w_min(i)))

        stop
        
    elseif xi_diag(i) == 0
    else
        error(['Error']);
    end

end

% 可以用穷举法，画出zeta的值在复平面的变化，取最小值

% 

% 
% stop

%% Plot

fn = 0;

w_p = logspace(-1,2,500)*2*pi;
w_pn = [-flip(w_p),w_p];
s_pn = 1i*w_pn;

fn = fn+1;
figure(fn)
SimplexPS.nyquist_c(F,s_pn);
scatter(real(xi_diag),imag(xi_diag),'x','LineWidth',1.5); hold on; grid on;

if 1

T1 = 1/(D+J*s);
T2 = 1/D;
T3 = 1/(J*s);

F1 = -s/T1;
F2 = -s/T2;
F3 = -s/T3;
    
fn = fn+1;
figure(fn)
[p1,~] = SimplexPS.nyquist_c(F1,s_pn,'Color','b');
[p2,~] = SimplexPS.nyquist_c(F2,s_pn,'Color','r');
[p3,~] = SimplexPS.nyquist_c(F3,s_pn,'Color','y');
scatter(real(xi_diag),imag(xi_diag),'x','LineWidth',1.5); hold on; grid on;
legend([p1,p2,p3],{'1','2','3'})
    
end
% plot(T_Re,T_Im);

%% Output

Hinv

K_H

phi

xi

sigma

sigma_max

ListPowerFlow



stop