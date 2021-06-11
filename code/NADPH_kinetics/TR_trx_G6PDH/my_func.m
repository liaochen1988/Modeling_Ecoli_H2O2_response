function dx = my_func(t,x,GSSG,Gor)

NADPH = x(1);

v_gor = (267 * GSSG  + 6.55e5 * GSSG ^ 2) * NADPH * Gor ...
        / (2.22e-5 * GSSG + 9.7e-5 * NADPH + GSSG * NADPH + 0.022 * GSSG ^ 2 + 3.9e3 * GSSG ^ 2 * NADPH);

dx = zeros(1,1);
dx(1) = -v_gor;

end

