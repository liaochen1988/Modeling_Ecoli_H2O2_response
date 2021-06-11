clear all;
clc;

syms Km_ahp Ho k_ahp PA
f = sqrt((Km_ahp-Ho+k_ahp/PA)^2+4*Km_ahp*Ho);
taylor(f,Ho)