
GbusVIF_diag = diag(GbusVIF);
GbusVIF_diag = GbusVIF_diag(Order_New2Old);
Index = find(imag(GbusVIF_diag)<0);

GbusVIF_diag
Index

% CapIndex = find(ListLineNew(:,5)>1);
% BusIndex = ListLineNew(CapIndex,[1,5])