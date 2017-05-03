function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

%ADD RULES BELOW
function result = ddr1( model, trace, params, t )
    erl_new = trace(t).erl + (1-model.parameters.default.beta) * (( (1-model.parameters.default.w) * model.parameters.default.v1 + model.parameters.default.w * trace(t).v2 ) - trace(t).erl) * model.parameters.default.delta_t;
    result = {t+1, 'erl', erl_new};
end
function result = ddr2( model, trace, params, t )
    v2_new = trace(t).v2 - ( model.parameters.default.a2 * trace(t).d / model.parameters.default.d_max ) * model.parameters.default.delta_t;
    result = {t+1, 'v2', v2_new};
end
function result = ddr3( model, trace, params, t )
    d_new = trace(t).erl - model.parameters.default.erl_norm;
    result = {t+1, 'd', d_new};
end