clear all
clc
close all

s = [1 1 1 1 1 1 2 3 4 5 6 7 7 7 7 8 9 10 11 8 6];
t = [2 3 4 5 6 7 3 4 5 6 2 8 9 10 11 10 10 11 8 1 11];
G = graph(s,t);
h = plot(G)

[T,p] = minspantree(G);
highlight(h,T,'EdgeColor','r','LineWidth',1.5)