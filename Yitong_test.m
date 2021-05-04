clear all
clc
close all

A = [3,2,1,4,5,6,100,10000];
A = transpose(A);
[Diff,Pos] = CalDiffMax(A)