%Test functions for Quaternion

for i = 1:100
    quat1_vec = rand([1,4]);
    quat2_vec = rand([1,4]);
    quat1 = Quaternion(quat1_vec);
    quat2 = Quaternion(quat2_vec(1), quat2_vec(2), quat2_vec(3), quat2_vec(4));

    %Testing constructor
    assert(approx(quat1.to_array(), quat1_vec));
    assert(approx(quat2.to_array(), quat2_vec));

    %Testing generalized multiplication
    assert(approx(quat1.quaternion_multiply(quat2).to_array(), ...
        quatmultiply(quat1_vec, quat2_vec)));

    %Testing generalized inversion
    assert(approx(quat1.quaternion_inverse().to_array(), quatinv(quat1_vec)));

    %Testing normalization & magnitude
    assert(approx(quat1.quaternion_normalize().to_array(), quat1_vec / norm(quat1_vec)));
    assert(approx(quat2.quaternion_normalize().to_array(), quat2_vec / norm(quat2_vec)));
    quat1 = quat1.quaternion_normalize();
    quat2 = quat2.quaternion_normalize();
    quat1_vec = quat1_vec / norm(quat1_vec);
    quat2_vec = quat2_vec / norm(quat2_vec);
    
    %Testing conjugation
    assert(approx(quat1.quaternion_conjugate().to_array(), quatconj(quat1_vec)));

    %Testing q2rotvec
    rot = quat1.quaternion2rotation_vec();
    assert(approx(rot, rotvec(quaternion(quat1_vec))));

    %Testing rotvec2q
    matlab_q = quaternion(rot, 'rotvec');
    [matlab_qw, matlab_qi, matlab_qj, matlab_qk] = parts(matlab_q);
    assert(approx(Quaternion.rotation_vec2quaternion(rot).to_array(), [matlab_qw, matlab_qi, matlab_qj, matlab_qk]));

    %Testing rotation applications
    rotating_vector = rand([1, 3]);
    assert(approx(quat1.apply_rotation(rotating_vector), quatrotate(quatinv(quat1_vec), rotating_vector)));
end

function r = approx(v1, v2)
    r = all(isapprox(v1, v2) == 1, 'all');
end