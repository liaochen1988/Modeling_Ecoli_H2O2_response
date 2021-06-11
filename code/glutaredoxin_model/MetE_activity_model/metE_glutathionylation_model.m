function dy = metE_glutathionylation_model(t,x,kon,koff)

metE_SH = x(1);
GSSG = x(2);
metE_SSG = x(3);
GSH = x(4);

dy = zeros(4,1);
dy(1) = -kon*metE_SH*GSSG+koff*metE_SSG*GSH;
dy(2) = dy(1);
dy(3) = -dy(1);
dy(4) = -dy(1);

end

