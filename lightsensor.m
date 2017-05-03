function result = lightsensor(controller)

    % Define the pin that contains the lightsensor
    pin = 'A0';
    % Return a predicate to the model trace that contains the value of the
    % lightsensor
    result = predicate('lightvalue' ,readVoltage(controller, pin));
end