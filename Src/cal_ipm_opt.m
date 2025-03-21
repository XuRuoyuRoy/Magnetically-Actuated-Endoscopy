%% ========================================================================
%  Script Name:    cal_ipm_opt.m
%  Description:    Calculate inverse dynamics of ipm using optimization
%
%  Author:         Ruoyu Xu
%  Date:           2025-03-12
%  Version:        1.0
%
%  Example:        
%      x_ipm = cal_ipm_opt([0;0;0.5; eul2rotm([0 15/57.3 0])*[-680;0;0]],[0;0;0.3;10;0;0])
%
% Notes:             x0-- ini state of ipm 6x1, x_epm 6x1 state of epm      
%                   z is constant p=-0.1, p=p_ipm-p_epm，
%                  p_hat = p/norm(p); PM_D = 3*(p_hat*p_hat.')-eye(3);
%                  headc=PM_D*m_epm/norm(PM_D*m_epm);
%
%  Dependencies:   na
%
%  Copyright:      (c) Ruoyu Xu. All rights reserved.
% ========================================================================  
function x_ipm = cal_ipm_opt(x_epm,x0) 
    t_sim=tic;
    options = optimoptions('fmincon','Display','none','Algorithm','interior-point','MaxIterations',300,'OptimalityTolerance',1e-5);
    m_c=norm(x0(4:6));
    x_cost = @(optx)f_cost(optx,x_epm); 
    nonlcon = @(optx)normcon(optx,m_c);
    lb=[x0(1:2)-0.05;x0(3);x0(4:6)-m_c];
    ub=[x0(1:2)+0.05;x0(3);x0(4:6)+m_c];
    [x_opt,cost_val] = fmincon(x_cost,x0,[],[],[],[],lb,ub,nonlcon,options);
    x_ipm=[x_opt(1:3);x_opt(4:6)];
%     disp(['Cost: ', num2str(cost_val), ' after ', num2str(1000*toc(t_sim)), 'ms'])
end

%% subfun
function [c,ceq] = normcon(optx,m_c)
    c=[];
    ceq=norm(optx(4:6))-m_c;
end

function [PM_F,PM_tau] = ipm_force(p,m_epm,m_ipm)
    u0= 4*pi*10^(-7);
    m_epm_norm = norm(m_epm);
    m_ipm_norm = norm(m_ipm);
    m_epm_hat = m_epm/m_epm_norm;
    m_ipm_hat = m_ipm/m_ipm_norm;
    p_norm = norm(p);
    p_hat = p/p_norm;
    PM_D = 3*(p_hat*p_hat.')-eye(3);
    PM_Z = eye(3)-5*(p_hat*p_hat.');
    PM_B = u0*m_epm_norm/(4*pi*p_norm^3)*PM_D*m_epm_hat;
    PM_F = 3*u0*m_epm_norm*m_ipm_norm/(4*pi*p_norm^4)*(m_epm_hat*m_ipm_hat'+m_ipm_hat*m_epm_hat'+(m_ipm_hat'*PM_Z*m_epm_hat)*eye(3))*p_hat;
    PM_tau=cross(m_ipm,PM_B);
end

function x_cost = f_cost(optx,x_epm)
    p=optx(1:3)-x_epm(1:3);
    m_ipm=optx(4:6);
    m_epm=x_epm(4:6);
    [PM_F,PM_tau] = ipm_force(p,m_epm,m_ipm);
    x_cost=norm([PM_F(1:2);PM_tau]);
end

