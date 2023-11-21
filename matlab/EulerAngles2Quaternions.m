function q = EulerAngles2Quaternions( angles )

    phi = angles(:,1);
    theta = angles(:,2);
    psi = angles(:,3);

    q0 = cos(phi/2).*cos(theta/2).*cos(psi/2) + sin(phi/2).*sin(theta/2).*sin(psi/2);
    q1 = sin(phi/2).*cos(theta/2).*cos(psi/2) - cos(phi/2).*sin(theta/2).*sin(psi/2);
    q2= cos(phi/2).*sin(theta/2).*cos(psi/2) + sin(phi/2).*cos(theta/2).*sin(psi/2);
    q3 = cos(phi/2).*cos(theta/2).*sin(psi/2) - sin(phi/2).*sin(theta/2).*cos(psi/2);

    q = [q0,q1,q2,q3]';
end