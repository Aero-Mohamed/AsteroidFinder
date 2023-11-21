function ptp = Quaternions2EulerAngles(q0123)
    q0 = q0123(:,1);
    q1 = q0123(:,2);
    q2 = q0123(:,3);
    q3 = q0123(:,4);

    ptp(:,1) = (atan2(2.*(q0.*q1 + q2.*q3),1-2.*(q1.^2 + q2.^2))); %phi
    ptp(:,2) = asin(2.*(q0.*q2-q3.*q1)); %theta
    ptp(:,3) = atan2(2.*(q0.*q3 + q1.*q2),1-2.*(q2.^2 + q3.^2)); %psi

    % q0 = cos(phi/2).*cos(theta/2).*cos(psi/2) + sin(phi/2).*sin(theta/2).*sin(psi/2)
    % q1 = sin(phi/2).*cos(theta/2).*cos(psi/2) - cos(phi/2).*sin(theta/2).*sin(psi/2)
    % q2= cos(phi/2).*sin(theta/2).*cos(psi/2) + sin(phi/2).*cos(theta/2).*sin(psi/2)
    % q3 = cos(phi/2).*cos(theta/2).*sin(psi/2) - sin(phi/2).*sin(theta/2).*cos(psi/2)

    ptp = real(ptp);
end