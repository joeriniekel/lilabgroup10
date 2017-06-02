function result = breathingsensor(controller)

    % Define the pin that contains the lightsensor
    % pin = 'A1';
    % Return a predicate to the model trace that contains the value of the
    % lightsensor
    val = readVoltage(controller, 'A1');
    % disp(val)

    result = predicate('breathingvalue' ,val);
end
