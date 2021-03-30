% Author(s): Yitong Li

%% 
function [YbusDSS,YbusCell] = YbusCalc_s_dss(ListLine,w0) 

w = w0;

% Get the data
FB = ListLine(:,1);             % From bus number
TB = ListLine(:,2);             % To bus number

Rlist = ListLine(:,3);              % Resistance, R
Xlist = ListLine(:,4);              % Inductance, wL
Blist = ListLine(:,5);              % Capacitance, wC
Glist = ListLine(:,6);              % Conductance, G

Tlist = ListLine(:,7);              % Turns ratio, T

XLlist = ListLine(:,8);
AreaType = ListLine(:,9);       % AC or DC type

% Get number
N_Bus = max(max(FB),max(TB));    % Number of buses
N_Branch = length(FB);

% Calculate y
for n = 1:N_Branch
    
    
        R = Rlist(n);
        X = Xlist(n);
        B = Blist(n);
        G = Glist(n);
        XL = XLlist(n);
        
    
        if ( (R==0) && (X==0) && ( isinf(G) || isinf(B) || (XL==0) ) )
            error(['Error: short circuit, branch from ' num2str(FB(n)) ' to ' num2str(TB(n))]);
        else
            % ### Mutual branch
            if FB(n)~=TB(n)
                if ~( isinf(G) || isinf(B) )
                    error(['Error: Mutual branch ' num2str(FB(n)) num2str(TB(n)) ' contains G or B']);
                end
                if X == 0                       % R branch
                    A_RL = []; B_RL = []; E_RL = []; C_RL = []; D_RL = inv(R);
                else                            % RL or L branch
                    % KVL equations
                    % [v] = {[R] + [sL]}*[i]
                    % State equations
                    % d_i/dt = 1/L*[-R]*[i] + 1/L*[1]*[v]
                    A_RL = 1/(X/w) * [-R];
                    B_RL = 1/(X/w) * [1];
                    E_RL = [1];
                    % Output equations
                    % [id] = [1]*[i] + [0]*[v]
                    C_RL = [1];
                    D_RL = [0];
                end
                % Get the branch model
                Y_RL = dss(A_RL,B_RL,C_RL,D_RL,E_RL); 
                Ybranch{n} = Y_RL;
            % ### Self branch
            else
                % Error check
                if ~( (R==0) && (X==0) )         % GC branch, normally for self branch
                    error(['Error: Self branch ' num2str(FB(n)) num2str(TB(n)) ' contains R or X.']);
                end
                if B == 0                           % G branch
                    A_GC = []; B_GC = []; E_GC = [];
                    C_GC = []; D_GC = inv(G);
                else                                % GC or C branch
                    % KCL equations
                    % [i] = {[G] + [sC]}*[v]
                    % State equation
                    % [d_v]/dt = 1/C*[-G]*[v] + 1/C*[1]*[i]
                    A_GC = 1/(B/w)*[-G];
                    B_GC = 1/(B/w)*[1];
                    E_GC = [1];
                    % Output equation
                    % [v] = [1]*[v] + [0]*[i]
                    C_GC = [1];
                    D_GC = [0];
                end
                % Get the branch model
                Z_GC = dss(A_GC,B_GC,C_GC,D_GC,E_GC);
                Y_GC= SimplexPS.DssSwitchInOut(Z_GC,1);
                Ybranch{n} = Y_GC;
                % For self branch, connect load to it
                if isinf(XL)
                    Ybranch{n} = Ybranch{n};
                elseif (XL)==0
                    error(['Error: The inductive load is short-circuit. Please check QLi settings.']);
                else
                    % KVL equation for XL
                    % [v] = {[sL]}*[i]
                    % => State equation
                    % [d_i]/dt = 1/L*[1]*[v] 
                    A_XL = 0;
                    B_XL = 1/(XL/w)*[1];
                    E_XL = [1];
                    % => Output equation
                    % [i] = [1]*[i] + [0]*[v]
                    C_XL = [1];
                    D_XL = [0];

                    Y_XL = dss(A_XL,B_XL,C_XL,D_XL,E_XL);
                    Ybranch{n} = SimplexPS.DssSum(Ybranch{n},Y_XL);
                end
            end
        end
    
    Ybranch{n};
    
end

% Initialise YBus
Y0 = ss([],[],[],[0]);
for m = 1:N_Bus
    for n = 1:N_Bus
        YbusCell{m,n} = Y0;
    end
end

for k = 1:N_Branch     
    if(FB(k) ~= TB(k))    
        % Formation of the Off Diagonal Elements...
        YbusCell{FB(k),TB(k)} = SimplexPS.DssSum(YbusCell{FB(k),TB(k)}, - Ybranch{k}/Tlist(k));      % Mutual admittance is negative
        YbusCell{TB(k),FB(k)} = YbusCell{FB(k),TB(k)};
        
        % Formation of the Diagonal Elements...
        YbusCell{FB(k),FB(k)} = SimplexPS.DssSum(YbusCell{FB(k),FB(k)}, Ybranch{k}/(Tlist(k)^2));  % Self admittance is positive
        YbusCell{TB(k),TB(k)} = SimplexPS.DssSum(YbusCell{TB(k),TB(k)}, Ybranch{k});
    else
        % Formation of the Diagonal Elements...
        YbusCell{FB(k),TB(k)} = SimplexPS.DssSum(YbusCell{FB(k),TB(k)}, Ybranch{k});
    end 
    
end

YbusDSS = SimplexPS.DssArrange(YbusCell);      	% Whole state space form
 
end