function initR = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)

initR = k1p*k2p*Peptide.*Gsh.^2.*Grx./(k1p*Peptide+k2p*Gsh.^2+k2n*Gssg);

end

