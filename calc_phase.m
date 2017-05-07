function [ phase ] = calc_phase(prev_chest_pos,prev_chest_c,avg_chest_c,amplitude)

  % phi in range[-0.5; 0.5]
  % phase = positive number
  % reduced/instantaneous phase = in range[0,1]
  % phase = mod(n,1);     %the phase of the cycle: a number between 0 and 1

  % phase of breathing cycle   chest_c
  % 0 and 0.5                  'zero-line' = avg_chest_c
  % 0.25                       max
  % 0.75                       min
  
  if strcmp(prev_chest_pos, '1 in') && prev_chest_c > avg_chest_c
    %phase is between 0 and 0.25
    relative_chest_c = (prev_chest_c - avg_chest_c) / amplitude;
    phase = relative_chest_c / 4;
  elseif strcmp(prev_chest_pos, '1 in') && prev_chest_c < avg_chest_c
    %phase is between 0.75 and 0
    relative_chest_c = (avg_chest_c - prev_chest_c) / amplitude;
    phase = (relative_chest_c / 4) + 0.75;
  elseif strcmp(prev_chest_pos, '3 out') && prev_chest_c > avg_chest_c
    %phase is between 0.25 and 0.5
    relative_chest_c = (prev_chest_c - avg_chest_c) / amplitude;
    phase = (relative_chest_c / 4) + 0.25;
  else
    % strcmp(prev_chest_pos, '3 out') && prev_chest_c < avg_chest_c
    %phase is between 0.25 and 0.5
    relative_chest_c = (avg_chest_c - prev_chest_c) / amplitude;
    phase = (relative_chest_c / 4) + 0.5;
  end
end
