clear all;
clc;

% cell volume parameters
B = [0.1882, 1.6028];
Vol = @(b,x)  b(1)*exp(b(2)*x);

% growth rate
% MOPS complete, MOPS minimal, MOPS complete w/ methionine
GR = log(2)./([21.5, 56.3, 26.5]/60);

% OxyR
OxyR_mn = [1867, 726, 1332];
OxyR_conc = OxyR_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

% KatE
KatE_mn = [26	1967	15];
KatE_conc = KatE_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

% KatG
KatG_mn = [2814	 1908	1550];
KatG_conc = KatG_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

% Gor
Gor_mn = [3255	1055	1664];
Gor_conc = Gor_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

% GrxA
GrxA_mn = [1247	431	807];
GrxA_conc = GrxA_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

% trxB (thioredoxin reductase)
TrxR_mn = [8599	2980	4301];
TrxR_conc = TrxR_mn ./ 6.02e23 ./ (Vol(B, GR) * 1e-15) * 1e6; % uM

