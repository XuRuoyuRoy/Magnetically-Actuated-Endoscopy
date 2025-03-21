function [PM_F,PM_tau,PM_D,PM_Z] = cal_ipm_force(p,m_epm,m_ipm)
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

        