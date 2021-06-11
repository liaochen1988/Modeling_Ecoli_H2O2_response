function dH = calcHin(x,k_met,Ho,PA,k_kat,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp,n_ahp)

dH = k_met+(Ho-x)*PA-k_kat*x/(Km_kat+x)/(1+x/Ki_katG_h2o2*0)-k_ahp*x^n_ahp/(Km_ahp^n_ahp+x^n_ahp);

end