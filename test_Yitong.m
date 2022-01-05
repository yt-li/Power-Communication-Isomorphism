clear all
clc
close all

s = sym('s','real');
wr = sym('wr','real');
K = sym('K','real');

T = [1,1i;
     1,-1i];
M = [s,-wr;
     wr,s];
Gdq = M*M;
  
Gdq_ = T*Gdq*inv(T);
Gdq_ = simplify(Gdq_)
  
GW = Gdq + eye(2)*K^2;
% GW = Gdq_ + eye(2)*K^2;
  
F = GW(1,1)*GW(2,2) - GW(2,1)*GW(1,2);
   
F1 = subs(F,'s',1i*wr+1i*K);
F1 = simplify(F1)

F2 = subs(F,'s',1i*wr-1i*K);
F2 = simplify(F2)

F3 = subs(F,'s',-1i*wr+1i*K);
F3 = simplify(F3)

F4 = subs(F,'s',-1i*wr-1i*K);
F4 = simplify(F4)