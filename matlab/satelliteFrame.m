function [dF, dM] = satelliteFrame(w123_rw_dot)

    reactionWheelParams;
    
    global Inertia_wr1_B Inertia_wr2_B Inertia_wr3_B;
    
    % negative sign due to the fact that the satellite move opposite to the
    % reaction wheel rotation direction.
    dM = -(Inertia_wr1_B*w123_rw_dot(1)*n1 + Inertia_wr2_B*w123_rw_dot(2)*n2 + Inertia_wr3_B*w123_rw_dot(3)*n3);
    dF = [0, 0, 0]';
end