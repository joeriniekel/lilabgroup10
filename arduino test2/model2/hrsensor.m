function result = hrsensor(controller)

    % Define the pin that contains the lightsensor
    pin = 'A0';
    % Return a predicate to the model trace that contains the value of the
    % lightsensor
    result = predicate('hr2' ,readVoltage(controller, pin));
end