% bode plot for complex transfer function in lnear form

% Authors(s): Yitong Li

%%
function Xw = bode_c_linear(X,sbd,varargin)

    funcX = matlabFunction(X);
 
 	Xw = zeros(1,length(sbd));
    for n = 1:length(sbd)
        Xw(n) = funcX(sbd(n));
    end
        
    fbd = imag(sbd/2/pi);
    arg_Xw = angle(Xw)/pi*180;    % Degree
    mag_Xw = abs(Xw);
    
    subplot(2,1,1)
    plot(fbd,mag_Xw,'linewidth',1); grid on; hold on;
    ylabel('Magnitude')
    subplot(2,1,2)
    plot(fbd,arg_Xw,'linewidth',1); grid on; hold on;
    ylabel('Phase (Degree)')
    xlabel('Frequency')
    
end
    
    


