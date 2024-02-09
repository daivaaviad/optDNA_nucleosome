function R = quaternion_to_rot(quat)

% Computes the rotation matrix for a given quaternion
   
    quat = quat/norm(quat);

    q1q1 = quat(1)*quat(1);
    q1q2 = quat(1)*quat(2);
    q1q3 = quat(1)*quat(3);
    q1q4 = quat(1)*quat(4);

    q2q2 = quat(2)*quat(2);
    q2q3 = quat(2)*quat(3);
    q2q4 = quat(2)*quat(4);

    q3q3 = quat(3)*quat(3);
    q3q4 = quat(3)*quat(4);

    q4q4 = quat(4)*quat(4);


    d1 = [q1q1 - q2q2 - q3q3 + q4q4
        2*(q1q2 + q3q4)
        2*(q1q3 - q2q4)];

    d2 = [2*(q1q2 - q3q4)
        -q1q1 + q2q2 - q3q3 + q4q4
        2*(q2q3 + q1q4)];

    d3 = [2*(q1q3 + q2q4)
        2*(q2q3 - q1q4)
        -q1q1 - q2q2 + q3q3 + q4q4];

    R = [d1 d2 d3];
end