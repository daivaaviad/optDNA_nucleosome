function quatval = rot_to_quaternion(Q)

% Computes the quaternion vector from a rotation matrix Q

tr = Q(1,1) + Q(2,2) + Q(3,3);

if (tr > 0)
    S = sqrt(tr+1.0) * 2;
    q0 = 0.25 * S;
    q1 = (Q(3,2) - Q(2,3)) / S;
    q2 = (Q(1,3) - Q(3,1)) / S;
    q3 = (Q(2,1) - Q(1,2)) / S;
elseif ((Q(1,1) > Q(2,2)) && (Q(1,1)> Q(3,3)))
    S = sqrt(1.0 + Q(1,1) - Q(2,2) - Q(3,3)) * 2;
    q0 = (Q(3,2) - Q(2,3)) / S;
    q1 = 0.25 * S;
    q2 = (Q(1,2) + Q(2,1)) / S;
    q3 = (Q(1,3) + Q(3,1)) / S;
elseif (Q(2,2) > Q(3,3))
    S = sqrt(1.0 + Q(2,2) - Q(1,1) - Q(3,3)) * 2;
    q0 = (Q(1,3) - Q(3,1)) / S;
    q1 = (Q(1,2) + Q(2,1)) / S;
    q2 = 0.25 * S;
    q3 = (Q(2,3) + Q(3,2)) / S;
else
    S = sqrt(1.0 + Q(3,3) - Q(1,1) - Q(2,2)) * 2;
    q0 = (Q(2,1) - Q(1,2)) / S;
    q1 = (Q(1,3) + Q(3,1)) / S;
    q2 = (Q(2,3) + Q(3,2)) / S;
    q3 = 0.25 * S;
end

q1 = q1;
q2 = q2;
q3 = q3;
q4 = q0;

quatval = [q1;q2;q3;q4];

end

