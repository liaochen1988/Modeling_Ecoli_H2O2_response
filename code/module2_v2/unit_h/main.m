clear;
clc;

hostparams = readParameters();
[nutr,growthRate] = run_nutrient_limitation(hostparams);
pchip(nutr,growthRate,0.12*3600)