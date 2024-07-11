alpha = 2;
s = 1/16;

q1 = s .* (1/2 - 1/8 * alpha .* s)
q2 = s .* (1/2 + 1/8 * alpha .* s)

L = q2 - q1;
s_squared = s.^2;

L_again = 4 * L / s_squared
