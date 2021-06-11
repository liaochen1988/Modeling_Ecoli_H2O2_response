function dH = calcHin(x,k_met,Ho,PA,k_cat,k_ahp,Km_ahp,n_ahp)

dH = k_met+(Ho-x)*PA-k_cat*x-k_ahp*x^n_ahp/(Km_ahp^n_ahp+x^n_ahp);

end