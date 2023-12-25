function dState_dT = Satellite(~, state, forces, moments)
    global BB; % Magnitic Field Body Frame
    global magFieldTimeStep currentTime lastMagFieldTime;
    global Inertia_rw1_cg Inertia_rw2_cg Inertia_rw3_cg;
    global w123_dot_rw_current; % Reaction Wheel Command Signal From the Controller
    
    planet;             % Planet Variables
    satellite_inertia;  % Mass & Inertia
    reactionWheelParams;% 3 Reaction wheels parameters
    
    % Extract States
    x = state(1); y = state(2); z = state(3);
    p = state(11); q = state(12); r = state(13);
    pqr = state(11:13);
    q_0123 = state(7:10);
    q0 = q_0123(1); q1 = q_0123(2); q2 = q_0123(3); q3 = q_0123(4);
    w123_rw = state(14:16); % reaction wheel anglure velocity
    
    
    %%% Kinamatics
    vel = state(4:6);
    q_dot = 0.5*[ 0 -p -q -r;
                  p  0  r -q;
                  q -r  0  p;
                  r  q -p  0]*q_0123;
    
    r = state(1:3);
    rho = norm(r);
    r_hat = r/rho;
    
    %%% Dynamics
    F_grav = -(Mu*m/rho^2)*r_hat;
    F = F_grav + forces;
    accel = F/m;
    
    LMN_magtorquers = [0;0;0];
    Moments = LMN_magtorquers + moments;
    
    % angular momentum
    H = Inertia_satellite * pqr + Inertia_rw1_cg*w123_rw(1)*n1 + Inertia_rw2_cg*w123_rw(2)*n2 + Inertia_rw3_cg*w123_rw(3)*n3;
    pqrdot = invI*(Moments - cross(...
        pqr, H...
    ));
    
    %% Magnitic Field | 
    % sample this reading every n steps for fast simulation execution.
    if(currentTime >= lastMagFieldTime)
        %%% Magnetic Field 
        phiE = 0;
        thetaE = acos(z/rho);
        psiE = atan2(y, x);
        latitude = 90 - thetaE*180/pi;
        longitude = psiE*180/pi;
        altitude_km = (rho)/1000; %km

        [BN, BE, BD] = igrf('01-Jan-2020', latitude, longitude, altitude_km, 'geocentric');
        B_NED = [BN; BE; -BD]; % North_East_Down Frame
        [S, C, ~] = SCT([phiE, thetaE+pi, psiE]);

        % MagField in Inertial Frame
        BI = [
            C.theta*C.epsi, (S.phi*S.theta*C.epsi - C.phi*S.epsi), (C.phi*S.theta*C.epsi + S.phi*S.epsi);
            C.theta*S.epsi, (S.phi*S.theta*S.epsi + C.phi*C.epsi), (C.phi*S.theta*S.epsi - S.phi*C.epsi);
            -S.theta, S.phi*C.theta, C.phi*C.theta
        ] * B_NED;
                
        Rotation_B2I = [(q0^2+q1^2-q2^2-q3^2) 2*(q1*q2-q0*q3)  2*(q0*q2+q1*q3);
              2*(q1*q2+q0*q3) (q0^2-q1^2+q2^2-q3^2) 2*(q2*q3-q0*q1);
              2*(q1*q3-q0*q2) 2*(q0*q1+q2*q3) (q0^2-q1^2-q2^2+q3^2)];
        BB = Rotation_B2I'*BI; % MagField in Body Frame
        BB = BB*(1e-9); % In Tesla
        lastMagFieldTime = currentTime + magFieldTimeStep;
    end
    
    
    
    %% return State    
    dState_dT = [vel; accel; q_dot; pqrdot; w123_dot_rw_current];
end