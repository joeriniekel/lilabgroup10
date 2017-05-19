function result = bsensor(controller)

    % Define the pin that contains the lightsensor
    pin = 'A1';
    % Return a predicate to the model trace that contains the value of the
    % lightsensor
    result = predicate('b2' ,readVoltage(controller, pin));
end