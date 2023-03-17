% % Messner & Penot
clear ;
close all ; 
clc ;

% % Initialisation des paramètres
fe = 1e4 ; % Fréquence d ’ échantillonnage
M = 4; % Nombre de symboles dans la modulation
n_b = log2 ( M ) ; % Nombre de bits par symboles
% ... autres paramètres
Eg = % Energie du filtre de mise en forme
Ds = 1e3;
Ts = 0.001 ;
f0 = 2.5e3 ;
alpha = 0.5;
Tg = 4*Ts;
Fse = Ts*fe;
Ns = 5000;
sigmacarre = No/4;

