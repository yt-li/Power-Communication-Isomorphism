% clear all
% clc
% close all
% 
% K1_ = [3,-10,7;
%        3,0,-3;
%        2,-2,0];
% [phi1_,xi1_] = eig(K1_)
% 
% K2_ = [5,-2,-3;
%        -4,10,-6;
%        -14,-6,20];
% [phi2_,xi2_] = eig(K2_)


% v1 = phi(1,:)
% v2 = phi(2,:)
% 
% t1 = v2*transpose(v1)
% t2 = v1*transpose(v1)
% 
% v1 = phi(:,1);
% v2 = phi(:,2);
% 
% t1 = transpose(v1)*v1

H_ = diag([10000,1,1]);
Hinv_ = inv(H_);
KH_ = Hinv_*K
[~,xi_] = eig(KH_)
% KH__ = KH_ - diag(diag(KH_))
% [~,xi__] = eig(KH__)
