classdef l2<handle
    % L2 Class for simulating l2-matlab models (v1.1.1)
    %   model = L2('Folder') creates an l2 object of the files in 'Folder'.
    %   model = L2('Folder', struct(..)) creates an l2 object of the files 
    %   in 'Folder', with filenames as given in the structure.
    %   model = L2('Folder', {.., ..}) creates an l2 object of the files 
    %   in 'Folder', with filenames changed as given in the cell array.
    %
    % l2 Properties:
    %   model - Folder name containing the model
    %   trace - Simulation results
    %   sorts - Sorts as defined in sorts.l2
    %   predicates - Predicates as defined in predicates.l2
    %   parameters - Parameters as defined in parameters.l2
    %   scenarios - Scenarios as defined in scenarios.l2
    %   files - Structure containing the filenames used
    % 
    % l2 Methods:
    %   reset - Reset the model, optionally with new model/filenames
    %   simulate - Simulate the model for t steps
    %   plot - Plot the current trace
    %   export - Plot and save the current trace to file
    %   validate - Check if predicate conforms to model definition
    %   evaluate - Calculate fitness of a trace
    %   tune - Tune given parameters using an evolutionary algorithm
    % 
    % l2 Methods (Static):
    %   exists - Check if a particular value of predicate holds true
    %   getall - Retrieve all matching values of a predicate
    % 
    % For more information, see the <a href="matlab: 
    % web('http://www.few.vu.nl/~jmn300/l2-matlab/')">L2-Matlab</a> website.
    %
    % See also PREDICATE, PRED2STR, STR2PRED, PREDCMP
    
    
    properties(GetAccess = 'public', SetAccess = 'private')
        model;          % folder containing the model
        trace;          % simulation results
        sorts;          % default sorts.l2
        predicates;     % default predicates.l2
        parameters;     % default parameters.l2
        scenarios;      % default scenarios.l2
    end
    
    properties(GetAccess = 'public', SetAccess = 'public')
        %FILES Structure containing the filenames used by l2. Can be set
        %   when creating or reseting the model by specifying all filenames
        %   in a structure or any number of file-filename combinations in a
        %   cell array.
        %
        %   rules - Matlab file containing the dynamic rules
        %   sorts - Text file (.l2) containing sort definitions
        %   predicates - Text file (.l2) containing predicate definitions
        %   scenarios - Text file (.l2) containing scenario definitions
        %   parameters - Text file (.l2) containing parameter definitions
        %   fitness - Matlab file containing the fitness function
        %
        %   Example:
        %   Creating a model with an entire set of different filenames:
        %   model = new l2('folder', struct('rules','rls.m','sorts','srts.l2','predicates','prdcts.l2','scenarios','scnrs.l2','parameters','prmtrs.l2');
        %   Resetting the model with a new parameter and rules file:
        %   model.reset({'rules' 'rls.m' 'parameters' 'prmtrs.l2'});
        %
        %   See also RESET
        files = struct('rules','rules.m','sorts','sorts.l2','predicates','predicates.l2','scenarios','scenarios.l2','parameters','parameters.l2', 'fitness', 'fitness.m');
    end
    
    methods
        %% Constructor
        function obj = l2(varargin)
            warning off backtrace

            switch nargin
                case 0 % empty model
                    obj.model = '';
                    obj.trace = struct([]);
                    obj.sorts = struct([]);
                    obj.predicates = struct([]);
                    obj.parameters = struct([]);
                    obj.scenarios = struct([]);
                    
                case 1 % default files
                    try
                        curpath = cd(varargin{1});
                        obj.model = pwd;
                        cd(curpath)
                    catch
                        throw(MException('L2:FolderMissing', 'Model folder ''%s'' not found.', varargin{1}));
                    end
                    obj.trace = struct([]);
                    obj = obj.reset();
                    
                otherwise % custom files
                    try
                        curpath = cd(varargin{1});
                        obj.model = pwd;
                        cd(curpath)
                    catch
                        throw(MException('L2:FolderMissing', 'Model folder ''%s'' not found.', varargin{1}));
                    end
                    obj.trace = struct([]);
                    obj = obj.reset(varargin{2:end});
            end
        end
        
        %% Change/reset model
        function result = reset(this, varargin)
            %RESET Function to reset the model. Optionally, a new model
            %   folder or file names can be added as arguments. Afterwards,
            %   the trace will be cleared and external files have been
            %   reloaded, thus providing a fresh l2-model.
            %
            %   model.RESET - Keep existing model and filenames
            %   model.RESET('model', 'folder') - Change the model folder
            %   model.RESET(file, name, ..., fileX, nameX) - Change any of the filenames.
            %   Model folder and filenames can be changed simultaneously.
            % 
            %   See also FILES
            if nargin > 1 % extra arguments given, process these first
                if mod(length(varargin),2) ~= 0 % odd number of arguments, cannot be correct
                    throwAsCaller(MException('L2:RESET', 'Missing inputs, specify both file/model and name.'));
                end
                try
                    for i=1:2:nargin-1
                        if strcmp(varargin{i}, 'model')
                            if isdir([pwd filesep varargin{1}])
                                this.model = [pwd filesep varargin{1}];
                            elseif isdir(varargin{1})
                                this.model = varargin{1};
                            else
                                throwAsCaller(MException('L2:FolderMissing', 'Model folder ''%s'' not found.', varargin{1}));
                            end
                        else
                            this.files.(varargin{i}) = varargin{i+1};
                        end
                    end
                catch err
                    throwAsCaller(err)
                end
            end
            this.trace = struct([]);
            try
                this = this.readFiles();    
            catch err
                throwAsCaller(err);
            end
            result = this;
        end
        
        %% Global commands
        function result = simulate(this, t_max, varargin)
            %SIMULATE Simulate the model for t timesteps. Optionally, a
            %   specific scenario can be provided as a second parameter or a
            %   particular parameter set as a third parameters.
            %
            %   Example:
            %   model.SIMULATE(10) - Simulate the default scenario using
            %   default parameters for 10 timesteps
            %   model.SIMULATE(10, 'alternative') - Simulate the
            %   'alternative' scenario for 10 timesteps using default
            %   parameter set
            %   model.SIMULATE(10, 'default', 'extreme') - Simulate the
            %   default scenario for 10 timesteps using the 'extreme'
            %   parameter set
            
            scnToSet = 'default';
            paramToSet = 'default';
            if nargin>=3
                scnToSet = varargin{1};
            end
            if nargin>=4
                paramToSet = varargin{2};
            end
            
            % initialise trace and parameters
            if isstruct(paramToSet)
                params = paramToSet;
            elseif ~isempty(this.parameters) && ~isempty(this.parameters.(paramToSet))
                params = this.parameters.(paramToSet);
            else
                params = {};
            end
            
            % Initialise the scenario
            display('Setting the scenario...')
            this.trace = this.resetTrace(t_max);
%##### UNCOMMENT THE LINE BELOW AND LINE 10 IN step.m TO REVERT TO THE
%ORIGINAL VERSION
%             this.trace = setScenario(this, this.trace, scnToSet);
            
            % Get function handle to rules file
            prevDir = pwd;
            clean__ = onCleanup(@(~)cd(prevDir));
            cd(this.model);
            [~, rules, ~] = fileparts(this.files.rules);
            rules = str2func(['@' rules]);
            delete(clean__);
            
            % Construct trace by simulating the model with current settings
            h=waitbar(0,'Simulating..');
            clean__ = onCleanup(@()close(h));
            for t=1:t_max-1
                waitbar(t/t_max)             
                %apply rules
                this.trace = step( this, this.trace, params, rules, t, scnToSet );
            end
            delete(clean__);
            % Remove timepoints after t
            this.trace = this.trace(1:t_max);
            result = this.trace(1:t_max);
        end
        
        function p = plot(this, varargin)
            %PLOT Plot the current trace file. By default, all predicates
            %   are plotted, but optionally a cell array can be provided
            %   containing a specific set of predicates to be plotted. 
            %
            %   model.PLOT
            %   model.PLOT({'predicate2' 'predicate4'})
            %
            %   See also EXPORT
            
            % Get the predicates to plot
            predicates = fieldnames(this.trace);
            if nargin == 2
                if iscell(varargin{1})
                    predicates = varargin{1};
                else
                    throw(MException('L2:Plot', 'Unknown argument ''%s''.', varargin{1}));
                end
            end
            
            % Plot the trace
            [p, height] = this.plot_panel(predicates, 'Screen');
        end
        
        function export(this, filename, varargin)
            %EXPORT Export a plot of the current trace to a file. 
            %   The first argument contains the filename (filetype extension 
            %   is optional) and optionally a specific set of predicates to 
            %   be plotted can be provided as the last argument. The 
            %   resulting file will be stored in the model folder in a 
            %   separate directory named 'plots'.
            %
            %   model.EXPORT('filename')
            %   model.EXPORT('filename.pdf' {'predicate2' 'predicate4'})
            %
            %   See also PLOT
            
            % Get the predicates to plot
            predicates = fieldnames(this.trace);
            if nargin == 3
                if iscell(varargin{1})
                    predicates = varargin{1};
                else
                    throw(MException('L2:Plot', 'Unknown argument ''%s''.', varargin{1}));
                end
            end
            
            % Plot the trace
            [p, height] = this.plot_panel(predicates, 'File');
                        
            % Get (and create) the file directory
            [path, filename, ext] = fileparts(filename);
            filedir = [this.model filesep 'plots' filesep];
            if ~isempty(path)
                filedir = [filedir filesep path filesep];
            end
            if ~isdir(filedir)
                mkdir(filedir)
            end
            
            % Save the figure(s)
            wb=waitbar(0,'Saving..');
            clean__ = onCleanup(@(~)close(wb));
            if length(p) == 1
                    h = ['-h' num2str(height(1))];
                    p{1}.export([filedir filename ext], h);
                    close(p{1}.figure)
                    waitbar(1)
            else
                for i=1:length(p)
                    h = ['-h' num2str(height(i))];
                    p{i}.export([filedir filename '_' num2str(i) ext], h);
                    close(p{i}.figure)
                    waitbar(i/length(p))
                end
            end
            
            % Close the windows
            delete(clean__);
        end
        
        function [result, err] = validate(this, pred)
            %VALIDATE Check if a predicate conforms to model definition.
            %   Therefore, the name is checked with defined predicates and
            %   each of its values is matched with possible values for the
            %   corresponding sort as defined in the predicate.
            
            % empty predicate is ok
            if isempty(pred)
                result = true;
                err = [];
                return
            end
            
            % check if predicate is defined
            if ~any(strcmp(pred.name, fieldnames(this.predicates)))
                result = false;
                err = MException('L2:VALIDATE:PredNotDefined', 'Predicate ''%s'' not defined in ''%s''.', pred.name, this.model);
                return
            end
            
            % check for required number of arguments in predicate
            sorts = this.predicates.(pred.name);
            if ~iscell(sorts)
                sorts = {sorts};
            end
            if length(sorts) ~= length(pred.arg)
                result = false;
                err = MException('L2:VALIDATE:NrArg', 'Predicate ''%s'' should contain %i instead of %i arguments.', pred.name, length(this.predicates.(pred.name)), length(pred.arg));
                return
            end
            
            % check each value
            for i=1:length(pred.arg)
                value = pred.arg{i};
                sort = sorts{i};
                
                if isa(value, 'predicate')
                    nestedPred = value;
                    value = ['<' nestedPred.name '>'];
                    [validated, err] = this.validateValue(value, sort);
                    if ~validated
                        result = false;
                        causeErr = MException('L2:VALIDATE', 'Predicate ''%s''.', pred.name);
                        err = addCause(err, causeErr);
                        return
                    else
                        [validated, err] = this.validate(nestedPred);
                        if ~validated
                            result = false;
                            return
                        end
                    end
                else
                    [validated, err] = this.validateValue(value, sort);
                    if ~validated
                        result = false;
                        causeErr = MException('L2:VALIDATE', 'Predicate ''%s''.', pred.name);
                        err = addCause(err, causeErr);
                        return
                    end
                end
            end
            
            % Otherwise, the predicate is ok!
            result = true;
            err = [];
        end   
        
        function [result, err] = validateValue(this, value, sort)  
            err = [];
            if strcmp(sort, 'REAL')
                result = isnumeric(value);
            elseif strcmp(sort, 'BOOLEAN')
                result = islogical(value);
            elseif strfind(value, 'rtinput_')
                result = true;
            else
                % check if predicate is defined
                if ~any(strcmp(sort, fieldnames(this.sorts)))
                    result = false;
                    err = MException('L2:VALIDATE:SortNotDefined', 'Sort ''%s'' not defined in ''%s''.', sort, this.model);
                    return
                end
                result = any(strcmp(value, this.sorts.(sort)));
            end
            if ~result
                err = MException('L2:VALIDATE:Value', 'Value ''%s'' is not of sort ''%s''.', value, sort);
            end
        end
            
        
        function score = evaluate(this, varargin)
            %EVALUATE Evaluate the current model using fitness function
            %   When calling the function with no parameters, the current
            %   trace is evaluated using a custom fitness function in the
            %   model folder. If arguments are provided, the model is first
            %   simulated using those arguments, after which the resulting
            %   trace is evaluated.
            %
            %   model.EVALUATE() - Evaluate the current trace
            %   model.EVALUATE(varargin) - Simulate model using varargin
            %   and evaluate the resulting trace
            %
            % See also: SIMULATE, TUNE
            
            % Get function handle to fitness function
            prevDir = pwd;
            clean__ = onCleanup(@(~)cd(prevDir));
            cd(this.model);
            [~, fitness, ~] = fileparts(this.files.fitness);
            fitness = str2func(['@' fitness]);
            delete(clean__);
            
            % If no arguments are given, evaluate current trace
            if nargin == 1
                % A trace needs to be present to evaluate
                if isempty(this.trace)
                    throwAsCaller(MException('L2:FITNESS', 'Trace is empty, first simulate the model before calculating fitness.'));
                end
                   
                % Get the fitness score
                score = feval(fitness, this);
                
            % If only one argument is given    
            else
                % Simulate using default scenario and get fitness
                this.simulate(varargin{:});
                score = feval(fitness, this);
            end
        end
        
        result = tune(this, params, t_max, varargin)
        
        function result = getSort( this, predicate )
            result = this.predicates.(predicate);
        end
        

    end
    
    %% Static methods
    methods (Static)
        [exists, values] = exists( trace, t, varargin )
        values = getall( trace, t, varargin )
        result = pred2str( predicate )
        result = str2pred( string )
        result = predcmp( pred1, pred2 )
    end

    %% Hidden methods
    methods (Hidden)
        
        function result = readFiles(this)    
            try
                this.sorts = getfield(readL2([this.model filesep this.files.sorts], 'sorts'), 'default');
            catch err
                warning('Error reading sorts-file, using only REAL and BOOLEAN.')
                this.sorts = struct([]);
            end
            try
                this.predicates = getfield(readL2([this.model filesep this.files.predicates], 'predicates'), 'default');
            catch err
                throwAsCaller(err);
            end
            try
                this.parameters = readL2([this.model filesep 'parameters.l2'], 'parameters');
            catch err
                warning('Error reading parameters-file, using no parameters.')
                this.parameters = struct([]);
            end
            try
                this.scenarios = readL2([this.model filesep this.files.scenarios], 'scenarios');
            catch err
                throwAsCaller(err);
            end
            result = this;
        end
        
        function trace = resetTrace(this, t_max)
            trace=struct;
            fields = fieldnames(this.predicates);
            for i=1:size(fields)
            	trace(t_max).(fields{i}) = [];
            end           
        end
    end

    %% Hidden, static methods
    methods (Hidden, Static)
        fncs = getRules()
        structure = flattenStruct(input)
        [lop,fop,nite]=aga(ninfo,label, ... 
                           pop, ... 
                           ng,nm,nr,nn, goal, ... 
                           funique,fitfun,mutfun,reproduccio,ranfun,prifun)  
    end
    
    %% Hide methods from handle class, just for clearity
    methods (Hidden)
        function addlistener(varargin)
            addlistener@addlistener(varargin{:})
        end
        function eq(varargin)
            eq@eq(varargin{:})
        end
        function findobj(varargin)
            findobj@findobj(varargin{:})
        end
        function findprop(varargin)
            findprop@findprop(varargin{:})
        end
        function ge(varargin)
            ge@ge(varargin{:})
        end
        function gt(varargin)
            gt@gt(varargin{:})
        end
        function le(varargin)
            le@le(varargin{:})
        end
        function lt(varargin)
            lt@lt(varargin{:})
        end
        function ne(varargin)
            ne@ne(varargin{:})
        end
        function notify(varargin)
            notify@notify(varargin{:})
        end
    end
    
end

