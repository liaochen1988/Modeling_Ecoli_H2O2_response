function dy = mse_h2o2_removal_model_R(p,timeOD,h2o2OD,timeEM,h2o2EM,...
    k_o,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,k_diff,Vin,Vout,r_g0,r_mc,beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia)

y0 = [0, 3.44, 3/12, 1.1e7, 0, 1, 0, p(5)];
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,8),'NonNegative',1:8);
[t_sim,y_sim] = ode15s(@h2o2_removal_model_R,[0:0.01:20],y0,options,...
    p(4),k_o,p(2),p(1),Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,p(3),k_diff,Vin,Vout,r_g0,r_mc,...
    beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);

% dy1 = (pchip(t_sim,y_sim(:,2),timeEM)-h2o2EM)/max(h2o2EM);
% dy2 = (pchip(t_sim,y_sim(:,4),timeOD)-h2o2OD)/max(h2o2OD);
dy1 = (pchip(timeEM,h2o2EM,t_sim)-y_sim(:,2))/max(h2o2EM);
dy2 = (pchip(timeOD,h2o2OD,t_sim)-y_sim(:,4))/max(h2o2OD);
dy = [dy1;dy2];

end

