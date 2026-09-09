classdef Quaternion 
    %#codegen
    %Quaternions are thought of in the scalar-first form! [qw, qx, qy, qz]
    properties
        qw % Scalar part of the quaternion
        qx % First component of the vector part
        qy % Second component of the vector part
        qz % Third component of the vector part
    end

    methods
        function quaternion = Quaternion(varargin)
            %Constructor
            if nargin == 1
                args = varargin{1};
                quaternion.qw = args(1);
                quaternion.qx = args(2);
                quaternion.qy = args(3);
                quaternion.qz = args(4);
            elseif nargin == 4
                quaternion.qw = varargin{1};
                quaternion.qx = varargin{2};
                quaternion.qy = varargin{3};
                quaternion.qz = varargin{4};
            else
                quaternion.qw = 1;
                quaternion.qx = 0;
                quaternion.qy = 0;
                quaternion.qz = 0;
            end
               
        end
        function result = quaternion_multiply(self, quat2)
            %Performs quaternion multiplication between two quaternions
            w = self.qw*quat2.qw - self.qx*quat2.qx - self.qy*quat2.qy - self.qz*quat2.qz;
            x = self.qw*quat2.qx + self.qx*quat2.qw + self.qy*quat2.qz - self.qz*quat2.qy;
            y = self.qw*quat2.qy - self.qx*quat2.qz + self.qy*quat2.qw + self.qz*quat2.qx;
            z = self.qw*quat2.qz + self.qx*quat2.qy - self.qy*quat2.qx + self.qz*quat2.qw;
            result = Quaternion(w,x,y,z);
        end

        function result = quaternion_magnitude(self)
            %Returns the magnitude of the quaternion
            result = sqrt(self.qw^2 + self.qy^2 + self.qz^2 + self.qx^2);
        end

        function result = quaternion_conjugate(self)
            %Takes the conjugate of the quaternion
            result = Quaternion(self.qw, -self.qx, -self.qy, -self.qz);
        end

        function result = quaternion_normalize(self)
            %Normalizes the quaternion to be a unit quaternion
            magnitude = self.quaternion_magnitude();
            result = Quaternion(self.qw/magnitude, self.qx/magnitude, self.qy/magnitude, ... 
                self.qz/magnitude);
        end

        function result = quaternion_inverse(self)
            %Inverts the quaternion. Identical to conjugation for unit
            %quaternions
            magnitude = self.quaternion_magnitude();
            conjugate = self.quaternion_conjugate();
            result = Quaternion(conjugate.qw / magnitude^2, conjugate.qx / magnitude^2,... 
                conjugate.qy / magnitude^2, conjugate.qz / magnitude^2);
        end

        function result = apply_rotation(self, vec)
            %Uses quaternion multiplication to apply a rotation to a vector
            q = self.quaternion_normalize();
            pure_q = Quaternion(0, vec(1), vec(2), vec(3));
            pure_result = q.quaternion_multiply(pure_q.quaternion_multiply(q.quaternion_inverse()));
            result = [pure_result.qx, pure_result.qy, pure_result.qz];
        end

        function result = quaternion2rotation_vec(self)
            res = self.to_array();
            if(self.qw < 0)
                res = -res;
            end
            if(self.qw >= 1)
                result = zeros([1,3]);
                return
            end

            theta = 2 * acos(res(1));
            if(theta == 0)
                result = zeros([1,3]);
                return
            end
            result = theta * res(2:4) / sqrt(1 - res(1)^2);
            
        end

        function result = to_array(self)
            %Turns qw, qx, qy, qz into an array
            result = [self.qw, self.qx, self.qy, self.qz];
        end
        function s = string(self)
            %override to_string method
            s = self.to_array();
        end
    end

    methods(Static)
        function result = rotation_vec2quaternion(vec)
            %Returns the quaternion quivalent for a rotation vector
            angle = sqrt(vec(1)^2 + vec(2)^2 + vec(3)^2);
            axis = [vec(1)/angle, vec(2)/angle, vec(3)/angle];
            if angle > 0
                result = Quaternion([cos(angle/2), sin(angle/2) * axis(1), ...
                    sin(angle/2)*axis(2), sin(angle/2)*axis(3)]);
            else
                result = Quaternion([1, 0, 0, 0]); 
            end
        end
        function diff = quat_diff(first_quat, second_quat)
            %Returns the rotation needed to get from first_quat to
            %second_quat

            qd = second_quat.quaternion_multiply(first_quat.quaternion_inverse());

            diff = qd.quaternion2rotation_vec();
        end
    end
end