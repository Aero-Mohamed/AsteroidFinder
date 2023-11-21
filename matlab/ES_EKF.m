% ErrorState-ExtendedKalmanFilter

function [corrected_state, p_cov_corrected] = ES_EKF(varGyro, P_cov_est, state_est, y_k)
    global H_jac
    
    R = varGyro*eye(3, 3);
    K = (P_cov_est*H_jac')*inv(H_jac*P_cov_est*H_jac' + R);
    
    errState = K*(y_k' - H_jac*state_est');
    
    delta_q = errState(1:3);
    delta_q = quatnormalize([sqrt(1 - sum(delta_q.^2)), delta_q']);
    q_est = quatnormalize([sqrt(1 - sum(state_est(1:3).^2)), state_est(1:3)]);
    
    delta_pqr = errState(4:6);
    delta_T = errState(7:9);
    
    pqr_hat = state_est(4:6) + delta_pqr';
    T_hat = state_est(7:9) + delta_T';
    q_hat = quatmultiply(delta_q, q_est);
    
    corrected_state = [q_hat, pqr_hat, T_hat];
    
    % Propagate the uncertainty
    p_cov_corrected = (eye(9, 9) - K*H_jac)*P_cov_est;

end