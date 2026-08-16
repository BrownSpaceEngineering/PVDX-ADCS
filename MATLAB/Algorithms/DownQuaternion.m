% INPUTS:
%   position_eci   - (1 x 3): ECI coordinates of satellite position
%   q_current.     - (1 x 4): quaternion of current orientation in Body to ECI frame 
%   providence_eci - (1 x 3): ECI coordinates of Providence, Rhode Island
%
% RETURN:
%   q_want - (1 x 4): desired down pointing quaternion in Body to ECI frame
% 
% All tests pass
function q_want = DownQuaternion(position_eci, q_current, providence_eci)

if nargin == 0
    run_tests();
    return
end

%subtract the current Providence from current ECI to get the relative vector from satellite to Rhode island
nadir = providence_eci - position_eci;

if(nadir == 0)
    error("Nadir is 0")
end

%normalize this vector
nadir = nadir / norm(nadir);


%get the rotation matrix of body to ECI quaternion, in row major form
body_to_ECI_rotation_matrix = quat2rotm(q_current);
body_to_ECI_rotation_matrix = body_to_ECI_rotation_matrix.';

%CHECK: z_direction corresponds to sattelite radar
current_z = body_to_ECI_rotation_matrix(3,:);

%if current_z points along the nadir, then return the current q
if all(isapprox(current_z, nadir, "loose") == 1)
    q_want = q_current;
    return;
end

axis_q = cross(current_z, nadir);

%antiparallel case

if(all(isapprox(axis_q, 0, "loose")))
    %then just use x-axis since perp to both vectors
    axis_q = body_to_ECI_rotation_matrix(1,:);
end

axis_q = axis_q / norm(axis_q);

%clamping for floating point error
angle_q = acos(max(min(dot(current_z, nadir), 1), -1));

%get imaginary and real parts of quaternion
imaginary = (sin(angle_q / 2) * axis_q);
real = cos(angle_q / 2);

%get error quaternion
q_error = [real, imaginary];

%convert error quaternion to wanted quaternion
q_want = quatmultiply(q_error, q_current);

end

%TEST:
function run_tests()

providence_eci = [0, 0, 0];

positions = [
    1000,    0,    0;
       0, 1000,    0;
       0,    0, 1000;
     600,  600,  600;
    -500,  800, -300;
    7000, 1200, -900
];

initial_rotms = {
    eye(3), ...
    rotz_deg(30) * roty_deg(15), ...
    rotx_deg(90), ...
    roty_deg(45) * rotx_deg(60), ...
    rotz_deg(-70) * rotx_deg(20) * roty_deg(10)
};

for p = 1:size(positions, 1)
    position_eci = positions(p, :);
    expected_nadir = (providence_eci - position_eci) / norm(providence_eci - position_eci);

    for r = 1:numel(initial_rotms)
        q_current = rotm2quat(initial_rotms{r});

        q_want = DownQuaternion(position_eci, q_current, providence_eci);
       
        %quaternion should stay normalized
        assert(approx(norm(q_want), 1));

        %z-axis of resulting rotation should point at nadir
        rotm_want = quat2rotm(q_want);
        z_axis = rotm_want(:, 3)';
        assert(approx(z_axis, expected_nadir));

        %resulting rotation should be the *minimal* rotation from q_current:
        %rotation angle from q_current to q_want should equal the angle
        %between the old z-axis and the nadir direction (no extra twist)
        rotm_current = quat2rotm(q_current);
        old_z = rotm_current(:, 3)';
        expected_angle = acos(max(min(dot(old_z, expected_nadir), 1), -1));

        q_delta = quatmultiply(q_want, quatinv(q_current));
        delta_angle = 2 * acos(max(min(abs(q_delta(1)), 1), -1));

        assert(approx(delta_angle, expected_angle));
    end
end

disp('All DownQuaternion tests passed.');

end

function r = approx(v1, v2)
    r = all(isapprox(v1, v2, "loose") == 1, 'all');
end

function R = rotx_deg(t)
    R = [1, 0, 0; 0, cosd(t), -sind(t); 0, sind(t), cosd(t)];
end

function R = roty_deg(t)
    R = [cosd(t), 0, sind(t); 0, 1, 0; -sind(t), 0, cosd(t)];
end

function R = rotz_deg(t)
    R = [cosd(t), -sind(t), 0; sind(t), cosd(t), 0; 0, 0, 1];
end
