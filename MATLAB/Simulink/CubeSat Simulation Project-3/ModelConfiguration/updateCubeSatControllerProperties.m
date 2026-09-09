function updateCubeSatControllerProperties()
% updateCubeSatControllerProperties Set CubeSat Simulation Model Controller
% Properties

%  Copyright 2019-2024 The MathWorks, Inc.

    d = Simulink.data.connect('asbCubeSatModelData.sldd');
    initGains = d.get('gains');

    initKp = num2str(initGains.Kp);
    initKi = num2str(initGains.Ki);
    initKd = num2str(initGains.Kd);

    prompt = {'CubeSat Proportional Controller Gain (Kp):','CubeSat Integral Controller Gain (Ki):',...
              'CubeSat Derivative Controller Gain (Kd):'};
    title = 'CubeSat Controller Properties';
    dims = [1 50];
    definput = {initKp,initKi,initKd};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    userInput = inputdlg(prompt,title,dims,definput, opts);

    if ~isempty(userInput)
        gains.Kp = eval(userInput{1});
        gains.Ki = eval(userInput{2});
        gains.Kd = eval(userInput{3});
        d.set('gains', gains);
    end

end
