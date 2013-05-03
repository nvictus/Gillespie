function status = ssaprogbar( flag, ti, ~ )
% SSAPROGBAR Used to display a progress bar window during an ssa simulation
%   Used by the ssa function when ssaprogbar is passed to ssa as the
%   'OutputFcn' property. Pushing the stop button terminates the simulation
%   and ssa returns the time course data it collected before the button was
%   pushed.
%
%   See also SSA.

persistent HFIGURE

switch flag
    case 'update'
        ud = get(HFIGURE, 'UserData');
        if ud.stop
            status = 1;
            return;
        end
        waitbar(ti/ud.time_interval, HFIGURE);
    case 'init'
        h = waitbar(0, 'Simulation in progress...');
        pos = get(0,'DefaultUicontrolPosition');
        pos(1) = pos(1) - 15;            
        pos(2) = pos(2) - 20;            
        uicontrol( h, 'Style','pushbutton', ...
                   'String','Stop', ...
                   'Position',pos, ...
                   'Callback',@StopButtonCallback, ...
                   'Tag','stop');        
        ud.stop = 0;
        ud.time_interval = ti(2)-ti(1);
        set(h, 'UserData', ud);        
        HFIGURE = h;        
    case 'done'
        h = HFIGURE;
        HFIGURE = [];
        close(h);
end

status = 0;

end

function StopButtonCallback(~,~)
ud = get(gcbf, 'UserData');
ud.stop = 1;
set(gcbf, 'UserData', ud);
end