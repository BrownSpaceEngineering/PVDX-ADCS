function I = moment2current3axis(m, Imax)
    % moment2current3axis
    % Converts desired moment into currents for 3 orthogonal magnetorquers
    %
    % Inputs:
    %   m       : 3x1 magnetic moment [mx my mz] [A m^2]
    %   Imax    : scalar or 3x1 current limit [A] (optional)
    %
    % Output:
    %   I       : 3x1 coil currents [Ix Iy Iz] [A]
    
    %Hard-coded parameters (to be measured/calculated)
    %   n       : 3x1 turns per coil
    %   A       : 3x1 coil areas [m^2]
    %   L,r     : 3x1 rod geometry (ignored if mu_r=1)
    %   mu_r    : 3x1 permeability (1 for air-core)
    n = [1 1 1];
    A = [1 1 1]; 
    L = [1 1 1];
    r = [1 1 1];

    G = ones(3,1);
    for i = 1:3
        if mu_r(i) > 1
                %Demagnetization factor
                ratio = L(i)./r(i);
                Nd = 4*(log(ratio)-1) / (ratio^2 - 4*log(ratio));
                
                %Core Amplification
                G(i) = 1 + (mu_r(i)-1) / (1 + (mu_r(i)-1)*Nd);
        end 
    end
    
    % Convert to currents
    I = m ./ (n(:).*A(:).*G);
    
    % Saturation
    if nargin == 2
        I = max(min(I, Imax), -Imax);
    end

end
