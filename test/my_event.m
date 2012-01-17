function [value, isterminal, direction] = my_event(t,x)
% Event functions
value = [x(2) - 50;
         x(1) - 5];
isterminal = [0;0];
direction = [1;0];
end