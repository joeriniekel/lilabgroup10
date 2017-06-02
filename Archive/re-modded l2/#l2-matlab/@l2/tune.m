function result = tune(this, params, t_max, varargin)
    %TUNE Tune parameters with simulations of t timesteps
    %   Using en evolutionary algorithm, tune the parameters given
    %   on simulations of length t. Parameters to tune are given in
    %   a cell array, containing pairs of parameternames and
    %   possible values. Real numbered intervals are given with the
    %   minimal and maximal value, a specific set of numbers is
    %   provided as an array, boolean parameters or parameters that
    %   have a sort as possible values are paired with that sort's
    %   name as a string value. Optionally a specific scenario
    %   can be given, or a cell array containing multiple scenarios
    %   to tune parameters on. When multiple scenarios are given,
    %   the average fitness will be used as a fitness score.
    %   Moreover, as a last argument a different parameter set may
    %   be provided as the base parameters.
    %
    %   model.TUNE({'param1', [0 1], 'param2', [1:10], 'param3', 'SORT1'}, 5)
    %   model.TUNE({..}, 5, 'alternative')
    %   model.TUNE({..}, 5, {'default' 'alternative'})
    %   model.TUNE({..}, 5, 'default', 'extreme')
    %
    % See also: SIMULATE, EVALUATE
    
    %% Settings
    ng=30; %30             % number of generations
    np=100; %100              % population size, NOT directly used by aga but useful here
    nm=2;               % elite size
    nr=floor(np*0.2);   % number of mutants
    nn=floor(np*0.75);   % number of normal sons (the remaining are newcommers)
    goal=0;             % goal (if aga finds a value below goal, iterations will be stopped)
    info=1;            % Verbosity level
    label=10000;        % Just a label (an arbitray integer number) 
    rng('shuffle')      % we don't want repetability in the GA  
    
    %% Functions
    
    switch nargin
        case 3
            scenarioSet = {'default'};
            parameterSet = 'default';
        case 4
            scenarioSet = varargin{1};
            parameterSet = 'default';
        case 5
            scenarioSet = varargin{1};
            parameterSet = varargin{2};
    end
    
    % discard identical individuals
    funique=@(pop) pop; % ignored

    % fitness
    fitfun=@(x) score(x, this, t_max, params, scenarioSet, parameterSet);

    % mutation
    mutfun=@(x,f) mutate( this, params, x, f );

    % reproduction
    reproduccio=@(x,y) reproduce(x, y, params); 

    % random individual
    ranfun=@() birth(this, params);

    % print individuals
    prifun=@(x) printChild(x);
    
    %% Setup

    % initial population
    pop={};
    for i=1:np
        pop{i}=ranfun(); 
    end
    
    %% Run algorithm

    [popFinal,best,nite] =  aga(info,label,pop,...
                            ng,nm,nr,nn,goal,...
                            funique,fitfun,mutfun,reproduccio,ranfun,prifun);
    
    result = popFinal{1};
end

function fitness = score(x, model, t_max, params, scenarioSet, parameterSet)
    fitness = 0;
    
    % Set a complete set of parameters
    parameters = model.parameters.(parameterSet);
    for j=1:length(x)
        parameters.(params{j*2-1}) = x{j};
    end
    
    % Evaluate for each scenario and sum up all the scores
    for i=1:length(scenarioSet)
        scenario = scenarioSet{1};
        fitness = fitness + model.evaluate(t_max, scenario, parameters);
    end

    % Calculate average fitness
    fitness = fitness / length(scenarioSet);
end

function mutation = mutate( model, params, x, ~ )

    mutationSize = 0.1; % the size of numerical mutations
                
    mutation = cell(size(x));

    % Go through each parameter of the individual
    for i=1:length(x)
        p = x{i};           % the parameter
        v = params{i*2};    % the possible values
        
        % if v is a string, and thus a sort
        if ischar(v)
            if strcmp(v, 'BOOLEAN')
                mutation{i} = logical(round(rand()));
            else
                v = model.sorts.(v);
                mutation{i} = v{randi(length(v))};
            end
        
        % if v describes the possible range    
        else
            % if v descibes min and max value
            if length(v) == 2
                mutRange = (v(2) - v(1)) * mutationSize;
                mutation{i} =  x{i} + (2*mutRange*rand() - mutRange);
                mutation{i} = max(v(1), min(v(2), mutation{i}));
            else
                mutRange = ceil(length(v) * mutationSize);
                [~, iRange] = min(abs(v-x{i}));
                iMutation =  iRange + randi([-mutRange,mutRange]);
                iMutation = max(1, min(length(v), iMutation));
                mutation{i} = v(iMutation);
            end
        end
    end
end

function child = reproduce(x,y,params)
    child = cell(size(x));
    for i=1:length(x)
        if ischar(params{i*2})
            % pick value of single parent
            switch randi(2)
                case 1
                    child{i} = x{i};
                case 2
                    child{i} = y{i};
            end
        else
            childV = (x{i}+y{i})/2; % take the average
            
            %check if the value is ok
            %if length v == 2, range is given and average is ok
            %else, check for closest match in given values
            v = params{i*2};    % the possible values
            if length(v)>2
                vDif = abs(v-childV);
                vMin = min(vDif);
                iMin = find(vDif==min(vDif));
                childV = v(iMin(randi(length(iMin))));
            end
            child{i} = childV;
        end
    end
end

function child = birth(model, params)
    child = cell(1,length(params)/2);
    for i=1:length(params)/2;
        v = params{i*2};    % the possible values
        
        % if v is a string, and thus a sort
        if ischar(v)
            if strcmp(v, 'BOOLEAN')
                child{i} = logical(round(rand()));
            else
                v = model.sorts.(v);
                child{i} = v{randi(length(v))};
            end
        
        % if v describes the possible range    
        else
            % if v descibes min and max value
            if length(v) == 2
                child{i} = v(1) + (v(2)-v(1)).*rand();
            else
                child{i} = v(randi(length(v)));
            end
        end
    end
end

function printChild(x)
   x 
end