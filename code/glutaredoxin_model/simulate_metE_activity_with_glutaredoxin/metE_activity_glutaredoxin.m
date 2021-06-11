function dy = metE_activity_glutaredoxin(t,x,kon,koff,k1p,k1n,k2p,k2n,Grx)

metE_SH = x(1);
GSSG = x(2);
metE_SSG = x(3);
GSH = x(4);

dy = zeros(4,1);
dy(1) = -kon*metE_SH*GSSG+koff*metE_SSG*GSH+...
       (k1p*k2p*metE_SSG*GSH^2-k1n*k2n*metE_SH*GSH*GSSG)/(k1p*metE_SSG+k2n*GSSG+k1n*metE_SH*GSH+k2p*GSH^2)*Grx;
dy(2) = dy(1);
dy(3) = -dy(1);
dy(4) = -dy(1);

end

