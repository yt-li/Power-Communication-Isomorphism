if Enable_Plot_Symbolic

w_p = logspace(-1,1,500)*2*pi;
w_pn = [-flip(w_p),w_p];
s_pn = 1i*w_pn;

if Exist_Vbus == 1
figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_V_sym,s_pn);
scatter(real(diag(xi_V)),imag(diag(xi_V)),'x','LineWidth',1.5); hold on; grid on;
end
if Exist_Ibus == 1
figure_n = figure_n+1;
figure(figure_n)
SimplexPS.nyquist_c(F_I_sym,s_pn);
scatter(real(diag(xi_I)),imag(diag(xi_I)),'x','LineWidth',1.5); hold on; grid on;
end

end