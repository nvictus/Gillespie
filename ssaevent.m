function [status, te_out, xe_out, ie_out] = ssaevent( flag, t, y, event_fcn )
%SSAEVENT Used to determine the location of user-defined events in an ssa simulation
%   Used internally by ssa function when an events function is passed as
%   the 'EventFcn' property. Records the time and system state each time at
%   least one of a list of evaluated event expressions crosses zero.
%
%   The events function (event_fcn) is of the form:
%
%       [value, isterminal, direction] = events(t,y)
%   
%   just as for the MATLAB ODE solvers (see documentation for odeset).
%
%   value(i) = the value of the ith evaluated event expression
%   isterminal(i) = 1 then terminate simulation when value(i) crosses 0
%   direction(i)  = 0 then all zero crossings count, 
%                  +1 then value(i) must be increasing to trigger event, 
%                  -1 then value(i) must be decreasing to trigger event
%
% See also SSA, ODESET


persistent te xe ie value_prev

switch flag
    case 'update'
        % Evaluate events
        [value, is_terminal, direction] = event_fcn(t,y);
        
        % Event is triggered if value crosses zero in the right direction
        is_triggered = (direction>=0 & value>=0 & value_prev<0) |...
                       (direction<=0 & value<=0 & value_prev>0);           
        if any(is_triggered)
            te(end+1,:) = t;
            xe(end+1,:) = y.';
            ie(end+1,:) = find(is_triggered);
            if any(is_triggered & is_terminal)
                status = 1;
                return;
            end
        end        
        value_prev = value;
        
    case 'init'
        te = [];
        xe = [];
        ie = [];
        [value, is_terminal, direction] = event_fcn(t,y);
        
        % Non-directional event may already be triggered
        is_triggered = (value == 0 & direction == 0);
        if any(is_triggered)
            te = t;
            xe = y.';
            ie = find(is_triggered);
            if any(is_triggered & is_terminal)
                status = 1;
                return;
            end
        end
        value_prev = value;
        
    case 'get'
        te_out = te;
        xe_out = xe;
        ie_out = ie;
        
    case 'done' %clean up
        te = [];
        xe = [];
        ie = [];
        value_prev = [];
end

status = 0;
end

