avg_chest_c = 50;

sig = 4;    %significance
  prev_f = 2;
  prev_t = 1;
  prev_phi = 0.1;
  prev_A = 10;

  % x1 = x(prev_t) = A * sin(2*pi * prev_f * prev_t + 0);
  x1 = prev_A * sin(2*pi * prev_f * prev_t + 0);
  x1 = round(x1,sig);
  if x1< 0.001, x1 = 0; end;

  %but: fnew
  f = 3;
  t = 2;
  % x2 = x(prev_t) = A * sin(2*pi * f * prev_t + phi);
  % Thus: x1 = x2 = x(prev_t)

  % however, phi is unknown
  % subtract phi:
  % new formula: omschrijven:
  % x2 = x1 = A * sin(2*pi*new_f*prev_t + new_phi);
  % 2*pi*new_f*prev_t + new_phi = asin( x(prev_t) / A) )
  phi = asin( x1 / A) - 2*pi * f * prev_t;

  % controle:
  % x2 = x(prev_t) = A * sin(2*pi*new_f*prev_t + new_phi);
  x2 = A * sin(2*pi* f * prev_t + phi);
  x2 = round(x2,sig);
  if x2< 0.001, x2 = 0; end;

  % Thus:
  % the new formula to calculate chest_c is x2
  % x2(t) = A * sin(2*pi* f * t + phi);
  x2 = A * sin(2*pi* f * t + phi);
  if x2< 0.001, x2 = 0; end;
  chest_c = avg_chest_c + x2;