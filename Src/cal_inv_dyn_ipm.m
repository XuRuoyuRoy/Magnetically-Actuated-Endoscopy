%% ========================================================================
%  Script Name:    cal_inv_dyn_ipm.m
%  Description:    Calculate inverse dynamics of ipm
%
%  Author:         Ruoyu Xu
%  Date:           2025-03-11
%  Version:        1.0
%
%  Example:        
%      x_ipm = cal_inv_dyn_ipm([0;0;0.5; eul2rotm([0 15/57.3 0])*[-680;0;0]],[0;0;0.4;10;0;0])
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
function x_ipm = cal_inv_dyn_ipm(x_epm,x0) 
    iter_max=1000;% max iteration
    tol=1e-2;% tolerence
    f_d=[0;0;0]; % stable force
    tau_d=[0;0;0]; % stable torque
    p=x0(1:3)-x_epm(1:3);
    m_ipm=x0(4:6);
    m_epm=x_epm(4:6);
%     p=[0.0;0.0;-0.1];
%     m_ipm=10*[1;0;0];
    m_ipm_norm = norm(m_ipm);
    m_epm_norm = norm(m_epm);
    [f,tau]=ipm_force(p,m_epm,m_ipm);
    df=f-f_d;
    dtau=tau-tau_d;
    err=norm([df;dtau]);
    k=0;
    while err>=tol
        k=k+1;
        df=f-f_d;
        df(3)=0; % fix the z direction
        dtau=tau-tau_d;
        J_F_=jacob_ipm(p,m_epm,m_ipm);
        dx=-0.3*pinv(J_F_)*[df;dtau];
        
        dx(3)=0;
        dx=LimitJointState(dx,[0.05;0.05;0.05;0.05;0.05;0.05]);

        p(1:2)=p(1:2)+dx(1:2); % only update x and y of IPM

        
        m_ipm=m_ipm+dx(4:6);
        m_ipm=m_ipm_norm*m_ipm/norm(m_ipm);
%         p_hat=p/norm(p);
%         mtemp=(3*(p_hat*p_hat.')-eye(3))*m_epm;
%         m_ipm=m_ipm_norm*mtemp/norm(mtemp);
        
        err=norm([df;dtau]);
        [f,tau]=ipm_force(p,m_epm,m_ipm);
        if k>=iter_max
            disp('max iteration')
            break;
        end

    end
    x_ipm=[p+x_epm(1:3);m_ipm];
    disp(['Accuracy: ', num2str(err), ' after ', num2str(k), ' iterations'])
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

function jacob_out = jacob_ipm(p,m_epm,m_ipm)
    u0= 4*pi*10^(-7);
    epm_norm = (m_epm'*m_epm)^0.5;
    ipm_norm = (m_ipm'*m_ipm)^0.5;
    p_norm = (p'*p)^0.5;
    m_epm_hat = m_epm/epm_norm;
    m_ipm_hat = m_ipm/ipm_norm;
    p_hat = p/p_norm;
    PM_D = 3*(p_hat*p_hat')-eye(3);
    PM_Z = eye(3)-5*(p_hat*p_hat');
    PM_G = eye(3)- p_hat*p_hat';
    F_p = 3*u0*epm_norm*ipm_norm/(4*pi*p_norm^5)*(m_epm_hat*m_ipm_hat'*PM_Z+m_ipm_hat*m_epm_hat'*PM_Z ...
        +m_ipm_hat'*m_epm_hat*PM_Z-5*(p_hat*p_hat')*m_epm_hat*m_ipm_hat'*PM_G ...
        -5*(p_hat*p_hat')*m_ipm_hat*m_epm_hat'*PM_G-5*m_ipm_hat'*(p_hat*p_hat')*m_epm_hat*PM_Z);
    F_mc = 3*u0*epm_norm*ipm_norm/(4*pi*p_norm^4)*(m_epm_hat'*p_hat*eye(3)+m_epm_hat*p_hat'+p_hat*m_epm_hat'*PM_Z);
    T_p = 3*u0*epm_norm*ipm_norm/(4*pi)*(skew(m_ipm_hat/p_norm^3)*(p_hat*m_epm_hat'*(PM_G/p_norm) ...
        +(PM_G/p_norm)*(p_hat'*m_epm_hat))+skew(PM_D*m_epm_hat)*(m_ipm_hat*p_hat'/p_norm^4));
    T_mc = -u0*epm_norm*ipm_norm/(4*pi*p_norm^3)*skew(PM_D*m_epm_hat);
    jacob_out=[F_p F_mc;T_p T_mc];
end

function S = skew(V)
	S = [
   	    0    -V(3)    V(2)
        V(3)      0    -V(1)
       -V(2)    V(1)      0
   ];
end
function  dq_out  = LimitJointState(dq_in,dqm)
    sca=abs(dq_in./dqm);
    sca_max=max(sca);
    if sca_max>1
        dq_out=dq_in/sca_max;
    else
        dq_out=dq_in;
    end
end
