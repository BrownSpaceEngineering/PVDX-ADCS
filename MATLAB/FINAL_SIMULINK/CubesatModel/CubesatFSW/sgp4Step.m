function [r_gcrf, v_gcrf, oe_osc, errCode] = sgp4Step(oe_epoch, bstar, epochJD, dt)
%#codegen
% sgp4Step  Simulink MATLAB Function block entry point for SGP4 (10 Hz)
%
% Hold the last uplinked elements on oe_epoch/bstar/epochJD and drive dt
% with (current UTC - element epoch) in seconds. sgp4Init reruns only
% when the held elements change (i.e. once per uplink); every other step
% only runs sgp4Propagate.
%
% See sgp4Init / sgp4Propagate for units.

    persistent sat lastKey
    key = [oe_epoch(:); bstar; epochJD];
    if isempty(sat) || ~isequal(key, lastKey)
        sat = sgp4Init(oe_epoch, bstar, epochJD);
        lastKey = key;
    end

    [r_gcrf, v_gcrf, oe_osc, ~, ~, ~, errCode] = sgp4Propagate(sat, dt);
end
