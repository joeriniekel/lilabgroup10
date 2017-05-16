function result = lightsensor(controller)

    % Define the pin that contains the lightsensor
    pin = 'A0';
    % Return a predicate to the model trace that contains the value of the
    % lightsensor
    val = readVoltage(controller, pin);
    disp(val)
    %controller
    
    result = predicate('lightvalue' ,val);
end