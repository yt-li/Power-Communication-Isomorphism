% clear all
% clc
% close all
% 
% s = [1 1 1 1 1 1 2 3 4 5 6 7 7 7 7 8 9 10 11 8 6];
% t = [2 3 4 5 6 7 3 4 5 6 2 8 9 10 11 10 10 11 8 1 11];
% G = graph(s,t);
% h = plot(G)
% 
% [T,p] = minspantree(G);
% highlight(h,T,'EdgeColor','r','LineWidth',1.5)

% KH_ = transpose(KH)*KH
% [phi_,xi_] = eig(KH_)


% clear all
% clc
% close all
% 
% KH = [5,-2,-3;
%      -2,2,0;
%      -3,0,3];
 
% [phi,xi] = eig(KH)
% 
% KH_GraphData = digraph(KH,'omitselfloops');
% 
% KH_ = laplacian(KH_GraphData);
% [phi_,xi_] = eigs(KH)
% % [phi_,xi_] = eigs(A,2,'smallestabs')
% 
% % [phi_,xi_] = eigs(A,2,'smallestabs')
% 
% A = [1,1;
%      1,0];
% 
% [phi,xi] = eig(A)

t1 = [1,2,3];

t2 = t1.*t1