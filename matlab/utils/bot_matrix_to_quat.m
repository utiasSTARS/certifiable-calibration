function q = bot_matrix_to_quat(R)

    q = zeros(4,1);
    
    m00 = R(1,1); m01 = R(1,2); m02 = R(1,3);
    m10 = R(2,1); m11 = R(2,2); m12 = R(2,3);
    m20 = R(3,1); m21 = R(3,2); m22 = R(3,3);
    
    
    q(1) = sqrt( max( 0, 1 + m00 + m11 + m22 ) ) / 2;
    q(2) = sqrt( max( 0, 1 + m00 - m11 - m22 ) ) / 2;
    q(3) = sqrt( max( 0, 1 - m00 + m11 - m22 ) ) / 2;
    q(4) = sqrt( max( 0, 1 - m00 - m11 + m22 ) ) / 2;
    q(2) = copysign( q(2), m21 - m12 );
    q(3) = copysign( q(3), m02 - m20 );
    q(4) = copysign( q(4), m10 - m01 );

function a = copysign(a,b)
    a = sign(b) * abs(a);