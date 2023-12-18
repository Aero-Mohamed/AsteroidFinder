function [w123_rw_dot, error, integral]   = reactionWheelController(dt, state, previousError, integral)
    global acc_rw_max maxSpeed;

    Ki_rw = eye(3)*1;
    Kp_rw = eye(3)*100;
    Kd_rw = eye(3)*1;
    
    
    pqr = state(11:13)';
    w123_rw = state(14:16);
    %w123_rw_dot = [0.01; 0; 0]; % rad/s^2
    
    pqr_command = [0;0;0];
    error = pqr - pqr_command;
    
    integral = integral + error * dt;
    derivative = (error - previousError) / dt;
    
    w123_rw_dot = Kp_rw * error + Ki_rw * integral + Kd_rw * derivative;
    w123_rw_dot = w123_rw_dot';
    
    
    for idx = 1:3
       if abs(w123_rw(idx)) >= maxSpeed
           w123_rw_dot(idx) = 0;
       else
           if abs(w123_rw_dot(idx)) >= acc_rw_max
               w123_rw_dot(idx) = sign(w123_rw_dot(idx))*acc_rw_max;
           end
       end
    end
    
end