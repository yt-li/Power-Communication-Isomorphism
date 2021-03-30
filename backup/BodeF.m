% bode plot for complex transfer function in lnear form

% Authors(s): Yitong Li

%%
function Xw = BodeF(X,f_p,varargin)

    f_pn = [-flip(f_p),f_p];
    w_pn = f_pn*2*pi;
    s_pn = 1i*w_pn;
    f_p_log = log10(f_p);
    f_pn_log = [-flip(f_p_log),f_p_log];

    funcX = matlabFunction(X);
 
 	Xw = zeros(1,length(s_pn));
    for n = 1:length(s_pn)
        Xw(n) = funcX(s_pn(n));
    end
        
    arg_Xw = angle(Xw)/pi*180;    % Degree
    mag_Xw = abs(Xw);
    mag_Xw_log = 20*log10(mag_Xw);
    
    subplot(2,1,1)
    plot(f_pn_log,mag_Xw_log,'linewidth',1); grid on; hold on;
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(f_pn_log,arg_Xw,'linewidth',1); grid on; hold on;
    ylabel('Phase (Degree)')
    xlabel('Log(Frequency) Hz')
    
end
    
    


