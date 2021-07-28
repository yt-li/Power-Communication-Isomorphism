% This function is used to plot the wb-lambda contour

% Author(s): Yitong Li, Yunjie Gu

%%
clear all
clc
close all

%% Basic
f0 = 50;
w0 = 2*pi*f0;

fb = 50;
wb = 2*pi*fb;

x = sym('x','real');    % Real part of lambda
y = sym('y','real');    % Imaginary part of lambda

lambda_f = x + 1i*y;
lambda_w = 2*pi*lambda_f;

F = (1i*w0 - lambda_w)/(1i*(w0 + wb) - lambda_w);   % Low pass filter

%% Get the function
Real_F = real(F);
Imag_F = imag(F);

Real_Fun = matlabFunction(Real_F);
Imag_Fun = matlabFunction(Imag_F);

Mag_F_sqr = Real_F^2 + Imag_F^2;
Mag_F_sqr = Mag_F_sqr - 1/2;
Mag_Fun = matlabFunction(Mag_F_sqr);

Ang_F = atan2(Imag_F,Real_F);
Ang_Fun = matlabFunction(Ang_F);
Ang_F1 = Ang_F - pi/2/2;
Ang_F2 = -pi/2/2 - Ang_F;
Ang_Fun1 = matlabFunction(Ang_F1);
Ang_Fun2 = matlabFunction(Ang_F2);

%% Check the creterion
N_Point = 1000;
x = linspace(-150,150,N_Point);
y = linspace(-150,150,N_Point);
y = flip(y);

xlength = length(x);
ylength = length(y);

for n = 1:ylength
    for m = 1:xlength
        % Get the x,y with a very small perturbation, to avoid the NaN for
        % the case such as 0/0.
        x_ = x(m) + eps(x(m));
        y_ = y(n) + eps(y(n));
        
        % Calculate magnitude
        Mag_F_Value(n,m) = Mag_Fun(x_,y_);
        Mag_F_Sign(n,m) = sign(Mag_F_Value(n,m));
        
        % Calculate phase angle
        Ang_F_Value1(n,m) = Ang_Fun1(x_,y_);
        Ang_F_Value2(n,m) = Ang_Fun2(x_,y_);
        if Ang_F_Value1(n,m)>0 || Ang_F_Value2(n,m)>0
            Ang_F_Sign(n,m) = -1;
        else
            Ang_F_Sign(n,m) = 1;
        end
        
        % Calculate the whole forbidden area
        % Notes: -1 is the forbidden region, i.e., the area of lambda that
        % does not satisfy the requirement of bandwidth wb.
        if Mag_F_Sign(n,m)<0 || Ang_F_Sign(n,m)<0
            F_Sign(n,m) = -1;
        else
            F_Sign(n,m) = 1;
        end
    end
end

% Get F_Plot
F_Plot = zeros(ylength,xlength);
for n = 1:ylength
    for m = 1:xlength
        if F_Sign(n,m) == 1
            Position = [n,m];
            value = 1;
            mode = 2;
            Check = CheckNeighbor(F_Sign,Position,value,mode);
            if Check ~= 1
                F_Plot(n,m) = 1;
            end
        end
    end
end

% Get position
% for n = 1:ylength
%     Find = find(F_Plot(n,:)==1);
%     for i = 1:length(Find)
%         
%     end
% end

figure(1)
contour(x,y,F_Plot,[1,1])

stop

%% Output
Mag_F_Value;
Mag_F_Sign

Ang_F_Value1;
Ang_F_Value2;
Ang_F_Sign


F_Sign
F_Plot





% mag_fun(0,0)

% fsolve(fun_mag,0)


% Given y, i.e., the frequency, solve x.

% Mag_F = abs(F);
% Ang_F = angle(F);
% 
% Mag_F = simplify(Mag_F)

%% Backup

%figure(1)
%fcontour(fun_mag)

%stop
