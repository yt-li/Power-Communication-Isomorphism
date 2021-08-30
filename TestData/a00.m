clear all
clc
close all

Gamma_SimData_SG = load('Gamma_SimData_SG').Gamma_SimData_SG;

Time = Gamma_SimData_SG.time;
Index = find(Time == 7)
length(Time)

w_SG_ = Gamma_SimData_SG.signals(5).values;
for i = 1:length(w_SG_)
w_SG(i) = w_SG_(:,:,i);
end

Gamma_SimData_SG.time = Gamma_SimData_SG.time(Index:end);
Gamma_SimData_SG.signals(1).values = Gamma_SimData_SG.signals(1).values(Index:end,:);
Gamma_SimData_SG.signals(2).values = Gamma_SimData_SG.signals(2).values(Index:end,:);
Gamma_SimData_SG.signals(3).values = Gamma_SimData_SG.signals(3).values(Index:end,:);
Gamma_SimData_SG.signals(4).values = Gamma_SimData_SG.signals(4).values(Index:end,:);
Gamma_SimData_SG.signals(5).values = w_SG;
Gamma_SimData_SG.signals(6).values = 0;
Gamma_SimData_SG.signals(7).values = Gamma_SimData_SG.signals(7).values(Index:end,:);
Gamma_SimData_SG.signals(8).values = Gamma_SimData_SG.signals(8).values(Index:end,:);
Gamma_SimData_SG.signals(9).values = Gamma_SimData_SG.signals(9).values(Index:end,:);
Gamma_SimData_SG.signals(10).values = Gamma_SimData_SG.signals(10).values(Index:end,:);
