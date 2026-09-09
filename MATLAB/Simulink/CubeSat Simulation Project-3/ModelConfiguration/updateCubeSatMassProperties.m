function updateCubeSatMassProperties()
%updateCubeSatMassProperties Set CubeSat Simulation Model Mass Properties

%  Copyright 2019-2024 The MathWorks, Inc.

    d = Simulink.data.connect('asbCubeSatModelData.sldd');
    initVehicle = d.get('vehicle');

    initMass = num2str(initVehicle.mass);
    moi = num2str(initVehicle.inertia);

    initMOI = regexprep(['[' moi(1,:) ';' moi(2,:) ';' moi(3,:) ']'],' +',' ');

    prompt = {'CubeSat Mass (kg):','CubeSat Moments of Inertia:'};
    title = 'CubeSat Mass Properties';
    dims = [1 50];
    definput = {initMass,initMOI};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    userInput = inputdlg(prompt,title,dims,definput,opts);

    if ~isempty(userInput)
        vehicle.mass = eval(userInput{1});
        vehicle.inertia = eval(userInput{2});
        d.set('vehicle', vehicle);
    end
end
