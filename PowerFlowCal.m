% ### Power flow analysis
[PowerFlow,~,V_,I_,~,~,~,~] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
        % The Gauss-Seidel method is always used here.
% [PowerFlow] = SimplexPS.PowerFlow.PowerFlowNR(ListBus,ListLine,Wbase);                      % ???
% Form of PowerFlow{i}: P, Q, V, xi, w

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% ### For >90 test                                                                  % ??????
% Load convention
% wL = 0.1
% PowerFlow{1} = [-9, -14.359, 1, 0,       pi*2*50];
% PowerFlow{2} = [9,  -14.359, 1, -2.0218, pi*2*50];
% PowerFlowNew = PowerFlow;
% wL = 1
% PowerFlow{1} = [-0.9, -1.4359, 1, 0,       pi*2*50];
% PowerFlow{2} = [0.9,  -1.4359, 1, -2.0218, pi*2*50];
% PowerFlowNew = PowerFlow;

% For printting later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);
ListPowerFlowNew = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);
P = ListPowerFlowNew(:,2);

% Update V and I
[V,I] = PowerFlowUpdateVI(PowerFlowNew);

% Notes:
% The codes in this part are borrowed from the SimplexPS toolbox. The V and
% I are updated based on the new power flow.

% ### For test
VoltageTheta = ListPowerFlowNew(:,5);
Max_VoltageThetaDiff = CalDiffMax(VoltageTheta);
Max_VoltageThetaDiff = Max_VoltageThetaDiff/pi*180
% fprintf(['Max angle difference between voltages in the whole system: ' num2str(Max_VoltageThetaDiff) ' Degree.\n'])

% AngleDiff = angle(V(1)) - angle(I(2));
% AngleDiff = AngleDiff/pi*180                % In load convention