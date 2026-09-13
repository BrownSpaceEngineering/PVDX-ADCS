% INPUTS:
% r_eci - (3 x 1): ECI coordinates of satellite position
% v_eci - (3 x 1): ECI velocity of satellite position
% q_b2eci - (1 x 4): quaternion of current orientation in body to eci,ACTIVE convention
% providence_eci - (1 x 3): ECI coordinates of Providence, Rhode Island
%
% RETURN:
% q_tgtb - (1 x 4): error quaternion, with axis in body coordinates
%
function [q_tgtb, z_want] = pointing_error(r_eci, v_eci, q_b2eci, providence_eci)
    persistent q_want_prev 
    if isempty(q_want_prev)
        q_want_prev = [1 0 0 0];   % initial guess, any unit quaternion
    end
    z_want = providence_eci - r_eci;
    z_want = z_want / norm(z_want);

    y_want = cross(z_want, v_eci);
    y_want = y_want / norm(y_want);

    x_want = cross(y_want, z_want);
    x_want = x_want / norm(x_want);

    R_want_eci = [x_want, y_want, z_want];

    q_want_eci = rotm2quat(R_want_eci);
    if dot(q_want_eci, q_want_prev) < 0
        q_want_eci = -q_want_eci;
    end
    q_want_prev = q_want_eci;   % store for next timestep
   
    q_tgtb = quatmultiply(quatinv(q_b2eci), q_want_eci);
  
end






       
    


    
