function [ fncs ] = rules()
    % DO NOT EDIT
    fncs = l2.getRules();
    for i=1:length(fncs)
        fncs{i} = str2func(fncs{i});
    end
end

%ADD RULES BELOW

%DDR1a Effect of training on regulation ability
function result = ddr1a(trace, params, t)
    result = {};
    
    %If training takes place at the current time point increase the
    %emotion regulation ability
    training = trace(t).training.arg{1};
    if training
        for emotion_regulation_ability = l2.getall(trace, t, 'emotion_regulation_ability', {NaN})
            new_ability = emotion_regulation_ability.arg{1} + params.ability_training;
                
            result = {result{:} {t+1, 'emotion_regulation_ability', {new_ability}}};
        end
        
    %If no training occurs maintain the same emotion regulation ability
    else
        for emotion_regulation_ability = l2.getall(trace, t, 'emotion_regulation_ability', {NaN})
            new_ability = emotion_regulation_ability.arg{1};
            
            result = {result{:} {t+1, 'emotion_regulation_ability', {new_ability}}};
        end
    end
end

%DDR1b Effect of training on expectation about others
function result = ddr1b(trace, params, t)
    result = {};
    
    %If training takes place at the current time point decrease  your
    %expectatoin about others
    training = trace(t).training.arg{1};
    if training
        for expectation_of_others = l2.getall(trace, t, 'expectation_about_others', {NaN})
            new_expectation = expectation_of_others.arg{1} - params.expectation_training;

            result = {result{:} {t+1, 'expectation_about_others', {new_expectation}}};
        end
        
    %If no training occurs maintain the same expectation about others
    else
        for expectation_of_others = l2.getall(trace, t, 'expectation_about_others', {NaN})
            new_expectation = expectation_of_others.arg{1};

            result = {result{:} {t+1, 'expectation_about_others', {new_expectation}}};
        end
    end
end

%DDR2 Effect of activity engagement on social activities
function result = ddr2(trace, params, t)
    result = {};
    
    %Increase the amount of social activities enjoyed by the amount of
    %social activities engaged in.
    for activity_engagement = l2.getall(trace, t, 'activity_engagement', {NaN})
        engagement = activity_engagement.arg{1};
        
        for social_activities = l2.getall(trace, t, 'social_activities', {NaN})
            new_activities = social_activities.arg{1} + engagement;
            
            result = {result{:} {t+1, 'social_activities', {new_activities}}};
        end
    end
end

%DDR3 Effect of level of vulnerability, emotion regulation ability, social
%activities and expectation about others on feeling of loneliness
function result = ddr3(trace, params, t)
    result = {};
    
    for level_of_vulnerability = l2.getall(trace, t, 'level_of_vulnerability', {NaN})
        vulnerability = level_of_vulnerability.arg{1};
        
        for emotion_regulation_ability = l2.getall(trace, t, 'emotion_regulation_ability', {NaN})
            emotion_regulation = emotion_regulation_ability.arg{1};
            
            for expectation_about_others = l2.getall(trace, t, 'expectation_about_others', {NaN})
                expectation = expectation_about_others.arg{1};
                
                for social_activities = l2.getall(trace, t, 'social_activities', {NaN})
                    activities = social_activities.arg{1};
                        
                    %Initialize influence variable to convert the
                    %non-numerical value of level_of_vulnerability to a
                    %numerical value to calculate new_feeling with.
                    influence = 0;
                    
                    if strcmp(vulnerability, 'high')
                        influence = params.vulnerability_high;
                    else
                        influence = params.vulnerability_low;
                    end

                    %Sum value new_feeling is determined by the following
                    %equation.
                    new_feeling = influence + expectation - emotion_regulation - activities;

                    %Initialize loneliness variable (this value is not
                    %used).
                    loneliness = false;

                    %If the sum variable exceeds a certain value set the
                    %feeling of loneliness accordingly.
                    if new_feeling < params.loneliness_threshold
                        loneliness = false;
                    else
                        loneliness = true;
                    end

                    result = {result{:} {t+1, 'feeling_of_loneliness', {loneliness}}};
                end
            end
        end
    end
end

%DDR4a Effect of feeling of loneliness on alcoholism, performance and
%cardiovascular diseases
function result = ddr4a(trace, params, t)
    result = {};
    
    for feeling_of_loneliness = l2.getall(trace, t, 'feeling_of_loneliness', {NaN})
        feeling = feeling_of_loneliness.arg{1};
        
        for performance_measure = l2.getall(trace, t, 'performance', {NaN})
            performance = performance_measure.arg{1};
            
            %Initialize new_performance variable (this value is not used).
            new_performance = 'high';
            
            %If there is a feeling of loneliness the new work performance
            %will be low, otherwise the work performance will be high
            if feeling
                new_performance = 'low';
            else
                new_performace = 'high';
            end              

            result = {result{:} {t+1, 'performance', {new_performance}}};
        end
    end
end

%DDR4b Effect of feeling of loneliness on alcoholism
function result = ddr4b(trace, params, t)
    result = {};
    
    for feeling_of_loneliness = l2.getall(trace, t, 'feeling_of_loneliness', {NaN})
        feeling = feeling_of_loneliness.arg{1};

        for alcoholism_level = l2.getall(trace, t, 'alcoholism', {NaN})
            alcoholism = alcoholism_level.arg{1};
            
            %Initialize the new_alcoholism variable (this value is not actually used)
            new_alcoholism = alcoholism;

            %If there is a feeling of loneliness increase the amount of
            %alcohol that is consumed, otherwise decrease it.
            if feeling
                new_alcoholism = alcoholism * params.alcoholism_increase;
            else
                new_alcoholism = alcoholism * params.alcoholism_decrease;
            end                  

            result = {result{:} {t+1, 'alcoholism', {new_alcoholism}}};
        end
    end
end

%DDR4c Effect of feeling of loneliness on cardiovascular diseases
function result = ddr4c(trace, params, t)
    result = {};
    
    for feeling_of_loneliness = l2.getall(trace, t, 'feeling_of_loneliness', {NaN})
        feeling = feeling_of_loneliness.arg{1};

        diseases = trace(t).cardiovascular_diseases.arg{1};
        
        %Initialize the new_diseases variable (This value is used in case 
        %none of the conditions of the if-statement are true).
        new_diseases = false;

        %If cardiovascular diseases are already present they stay present
        %(not recoverable), otherwise check if there is a feeling of
        %loneliness. If there is a feeling of loneliness cardiovascular
        %diseases are present.
        if diseases
            new_diseases = diseases;
        elseif feeling
            new_diseases = true;
        end              

        result = {result{:} {t+1, 'cardiovascular_diseases', {new_diseases}}};
    end
end

%--------------------------------------------------------------------------

%ADR1 Observation of cardiovascular diseases
function result = adr1(trace, params, t)
    result = {};
    
    for cardiovascular_diseases = l2.getall(trace, t, 'cardiovascular_diseases', {NaN})
        diseases = cardiovascular_diseases.arg{1};
        
        result = {result{:} {t+1, 'observation', predicate('cardiovascular_diseases', diseases)}};
    end
end

%ADR2 Observation of performance
function result = adr2(trace, params, t)
    result = {};
    
    for performance = l2.getall(trace, t, 'performance', {NaN})
        level = performance.arg{1};
        
        result = {result{:} {t+1, 'observation', predicate('performance', level)}};
    end
end

%ADR3 belief of performance
function result = adr3(trace, params, t)
    result = {};
    
    for observation = l2.getall(trace, t, 'observation', {predicate('performance', NaN)})
        belief = observation.arg{1};
        
        result = {result{:} {t+1, 'belief', belief}};
    end
end

%ADR4 belief of cardiovascular diseases
function result = adr4(trace, params, t)
    result = {};
    
    for observation = l2.getall(trace, t, 'observation', {predicate('cardiovascular_diseases', NaN)})
        belief = observation.arg{1};
        
        result = {result{:} {t+1, 'belief', belief}};
    end
end

%ADR5 belief of feeling of loneliness
function result = adr5(trace, params, t)
    result = {};
    
    for cardio_belief = l2.getall(trace, t, 'belief', predicate('cardiovascular_diseases', NaN))
        cardio = cardio_belief.arg{1}.arg{1};
        for performance_belief = l2.getall(trace, t, 'belief', predicate('performance', NaN))
            performance = performance_belief.arg{1}.arg{1};
            
            %Initialize feeling_of_loneliness variavble (this value is not
            %used)
            feeling_of_loneliness = false;
            
            %If there is a belief of cardiovascular diseases that are
            %present and there is a belief of a high work performance then
            %there is a belief of a feeling of loneliness
            
            %If there is a belief of cardiovascular diseases that are 
            %present and there is a belief of a low work performance then 
            %there is a belief of no feeling of loneliness
            if cardio
                if strcmp(performance, 'low')
                    feeling_of_loneliness = true;
                else
                    feeling_of_loneliness = false;
                end
            else
                feeling_of_loneliness = false;                
            end
            
            result = {result{:} {t+1, 'belief', predicate('feeling_of_loneliness', feeling_of_loneliness)}};
            
        end
    end
end

%ADR6 assessment of feeling of loneliness
function result = adr6(trace, params, t)
    result = {};
    
    for loneliness_belief = l2.getall(trace, t, 'belief', predicate('feeling_of_loneliness', NaN))
        belief = loneliness_belief.arg{1}.arg{1};
        
        %for loneliness_desire = l2.getall(trace, t, 'desire', {NaN}) 
        for loneliness_desire = l2.getall(trace, t, 'desire', predicate('feeling_of_loneliness', NaN))
            desire = loneliness_desire.arg{1}.arg{1};
            
            %Initialize assessment variable (this value is not used)
            assessment = 0;
            
            %If the belief of the value of feeling of loneliness is not the
            %same as the desire for the feeling of loneliness an assessment
            %is made of a feeling of loneliness.
            if belief ~= desire
                assessment = predicate('feeling_of_loneliness', true);
            else
                assessment = predicate('feeling_of_loneliness', false);
            end

            result = {result{:} {t+1, 'assessment', assessment}};
        end
    end
end

%--------------------------------------------------------------------------

% from observe( emotion_regulation_ability ) 	to belief( emotion_regulation_ability )

function result = belief_emotion_regulation_ability(trace, params, t)
    result = {};
    
      for observe_emotion_reg_ability = l2.getall(trace, t, 'emotion_regulation_ability', {NaN})
        ability = observe_emotion_reg_ability.arg{1};
        
        result = {result{:} {t+1, 'belief', predicate('emotion_regulation_ability', ability)}};
    end
end

% from observe( expectation_about_others ) 	to belief( expectation_about_others )

function result = belief_expectation_about_others(trace, params, t)
    result = {};
    
      for observe_expectations_about_others = l2.getall(trace, t, 'expectation_about_others', {NaN})
        expectations = observe_expectations_about_others.arg{1};
        
        result = {result{:} {t+1, 'belief', predicate('expectation_about_others', expectations)}};
    end
end
    
function result = belief_level_of_vulnerability(trace, params, t)
    result = {};
    
      for observe_level_of_vulnerability = l2.getall(trace, t, 'level_of_vulnerability', {NaN})
        vulnerability = observe_level_of_vulnerability.arg{1};
        
        result = {result{:} {t+1, 'belief', predicate('level_of_vulnerability', vulnerability)}};
      end
end

function result = belief_social_activities(trace, params, t)
    result = {};
    
      for observe_social_activities = l2.getall(trace, t, 'social_activities', {NaN})
        activity = observe_social_activities.arg{1};
        
        result = {result{:} {t+1, 'belief', predicate('social_activities', activity)}};
      end
end



  
function result = belief_loneliness(trace, params, t)
    result = {};  
    
    for belief_level_of_vulnerability = l2.getall(trace, t, 'belief' , predicate('level_of_vulnerability', {NaN}))
        vulnerability = belief_level_of_vulnerability.arg{1};
        if strcmp(vulnerability, 'high');
            influence = params.vulnerability_high;
        else
            influence = params.vulnerability_low;
        end
        
        for belief_social_activities = l2.getall(trace, t, 'belief' , predicate('social_activities', {NaN}))
             social_activities = belief_social_activities.arg{1};
                  
            for belief_emotion_regulation_ability = l2.getall(trace, t, 'belief' , predicate('emotion_regulation_ability', {NaN}))
                  emotion_regulation_ability = belief_emotion_regulation_ability.arg{1};
                      
                for belief_expectations_about_others = l2.getall(trace, t, 'belief' , predicate('expectation_about_others', {NaN}))
                    expectations_about_others = belief_expectations_about_others.arg{1};

                    loneliness = ((influence - emotion_regulation_ability) + (expectations_about_others - social_activities));
                    result = {result{:} {t+1, 'belief', predicate('feeling', loneliness)}}; 
                end
            end
        end
    end
end
   
%from assessment( feeling_of_loneliness) to desire( feeling )


function result = desire_loneliness(trace, params, t)      
    result = {};  
    for assessment_loneliness = l2.getall(trace, t,  'assessment' , predicate('feeling_of_loneliness', {NaN}))    
        assessment = assessment_loneliness.arg{1}.arg{1};
        
        for belief_loneliness = l2.getall(trace, t,  'belief' , predicate('feeling', {NaN}))
            belief = belief_loneliness.arg{1}.arg{1};
         
            if assessment %disp('assessment is true');
                desired_feeling = belief;
            else    %disp('assessment is false');
                desired_feeling = belief - params.y1;
            end
            result = {result{:} {t+1, 'desire', predicate('feeling' , desired_feeling)}};
        end
    end 
end
             
%from desire( feeling ), belief( social_activities ) to desire( social_activities )

function result = desire_social_activities(trace, params, t)      
   result = {}; 

   for desire_feeling = l2.getall(trace, t,  'desire' , predicate('feeling', {NaN}))    
        desire_f = desire_feeling.arg{1}.arg{1};
        
         for belief_feeling = l2.getall(trace, t,  'belief' , predicate('feeling', {NaN}))    
             belief_f = desire_feeling.arg{1}.arg{1};
             
             for  belief_social_activities = l2.getall(trace, t,  'belief' , predicate('social_activities', {NaN}))    
                  belief_s = desire_feeling.arg{1}.arg{1};
             
                  difference_feeling = (desire_f - belief_f);
                  desire_s = ((belief_s + difference_feeling) * params.y4);
                  
                  result = {result{:} {t+1, 'desire', predicate('social_activities' , desire_s)}};
             end
        end
    end  
end

   

function result = desire_activity_engagement(trace, params, t)      
        result = {}; 
        
        for desire_s_activities = l2.getall(trace, t,  'desire' , predicate('social_activities', {NaN}))    
            desire_social_activities = desire_s_activities.arg{1}.arg{1};
            
            for belief_s_activities = l2.getall(trace, t,  'belief' , predicate('social_activities', {NaN}))    
                belief_social_activities = belief_s_activities.arg{1}.arg{1};
                
                desire_engagement = (desire_social_activities - belief_social_activities);
                
                result = {result{:} {t+1, 'desire', predicate('activity_engagement' , desire_engagement)}};
                
            end
        end
end
           
%from desire( activity_engagement ) to desire( msg_to_h_activity_engagement )

function result = desire_msg_to_h_activity_engagement(trace, params, t)      
    result = {}; 
    
    for desire_a_engagement = l2.getall(trace, t,  'desire' , predicate('activity_engagement', {NaN}))    
        desire_activity_engagement = desire_a_engagement.arg{1}.arg{1};
        if desire_activity_engagement <= params.activity_lvl_1;
            propose_message = 'stay';
        elseif  desire_activity_engagement < params.activity_lvl_2;
            propose_message = 'raise';
        elseif  desire_activity_engagement < params.activity_lvl_3;
            propose_message = 'raise_more';
        else
            propose_message = 'raise_max';
        end
        result = {result{:} {t+1, 'desire', predicate('message_to_h_activity_engagement' , propose_message)}};
    end
end                  
             
%from desire( msg_to_h_activity_engagement ) to propose( msg_to_h_activity_engagement )

function result = propose_msg_to_h_activity_engagement(trace, params, t)      
	result = {}; 
    for des_message_engagement = l2.getall(trace, t,  'desire' , predicate('message_to_h_activity_engagement', {NaN}))    
        msg_to_h = des_message_engagement.arg{1}.arg{1};        
        result = {result{:} {t+1, 'propose', predicate('message_to_h_activity_engagement' , msg_to_h)}};           
    end
end

