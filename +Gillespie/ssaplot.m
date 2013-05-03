function status = ssaplot( flag, t, y, indices )
%SSAPLOT Used to update a time course plot at every time step during an ssa simulation
%   Used by ssa function when the function ssaplot is passed to ssa as the
%   'OutputFcn' property. To plot only some species, pass in their indices
%   as the 'OutputSel' property. ssaplot works similarly to the odeplot
%   function used with the MATLAB ODE solvers and is based on that function.
%
%   See also SSA, ODEPLOT

persistent HFIGURE HAXES

chunk = 128;

switch flag
    
    case 'update'
        if ( isempty(HFIGURE) || isempty(HAXES) )
            error('SSAPLOT:NoInit',...
                  'SSAPLOT not initialized.');
        elseif ( ishghandle(HFIGURE) && ishghandle(HAXES) )

            try
                % Append new data
                userdata = get(HFIGURE, 'UserData');
                y = y(userdata.select);
                ntimepts = 2;
                old_i = userdata.index;
                new_i = old_i + ntimepts;
                chunk = max(chunk, ntimepts);
                [nrows, ncols] = size(userdata.y);
                if new_i > nrows
                    userdata.t = [userdata.t; zeros(chunk,1)];
                    userdata.y = [userdata.y; zeros(chunk,ncols)];
                end
                userdata.t(old_i+1:new_i) = [t; t];
                userdata.y(old_i+1:new_i,:) = [userdata.y(userdata.index,:); y.'];
                userdata.index = new_i;
                set(HFIGURE, 'UserData', userdata);
                
                % Plot new data
                ylim = get(HAXES, 'YLim');
                if ( old_i == 1 || min(y(:)) < ylim(1) || ylim(2) < max(y(:)) )
                    %axes out of range or just initialized -> replot everything
                    for j = 1:ncols
                        set(userdata.lines(j), 'XData', userdata.t(1:new_i),...
                                                'YData', userdata.y(1:new_i,j) );
                    end
                else
                    %extend curves
                    for j = 1:ncols
                        set(userdata.line(j), 'XData', userdata.t(old_i:new_i),...
                                               'YData', userdata.y(old_i:new_i, j));
                    end
                end
                drawnow;

            catch err
                error('SSAPLOT:ErrorUpdatingWindow',...
                      'Error updating plot window: %s', err.message);
            end
        end

    case 'init'
        userdata.select = indices;
        nselected = length(indices);       
        userdata.t = zeros(chunk,1);
        userdata.y = zeros(chunk, nselected);
        userdata.index = 1;
        userdata.t(1) = t(1);
        userdata.y(1,:) = y(indices).';
        
        f = figure(gcf);
        HFIGURE = f;
        HAXES = gca;
        
        t_init = userdata.t(1);
        y_init = userdata.y(1,:);
        
        if ~ishold
            userdata.lines = plot(t_init, y_init, '-');
            hold on;
            userdata.line = plot(t_init, y_init, '-', 'EraseMode', 'none');
            hold off
            set(HAXES, 'XLim', [min(t), max(t)]);
        else
            userdata.lines = plot(t_init, y_init, '-');
            userdata.line = plot(t_init, y_init, '-', 'EraseMode', 'none');
        end
        
        %stop button
        hButton = findobj(f,'Tag','stop');
        if isempty(hButton)
            userdata.stop = 0;
            pos = get(0,'DefaultUicontrolPosition');
            pos(1) = pos(1) - 15;            
            pos(2) = pos(2) - 15;            
            uicontrol( 'Style','pushbutton', ...
                       'String','Stop', ...
                       'Position',pos, ...
                       'Callback',@StopButtonCallback, ...
                       'Tag','stop');
        else
            set(hButton,'Visible','on');     % make sure it's visible
            if ishold
                old_ud = get(f,'UserData');
                userdata.stop = old_ud.stop; % don't change old ud.stop status
            else
                userdata.stop = 0;
            end
        end
        
        set(f, 'UserData', userdata);
        
    case 'done'
        f = HFIGURE;
        ax = HAXES;
        HFIGURE = [];
        HAXES = [];
        
        if ishghandle(f)
            userdata = get(f, 'UserData');
            userdata.t = userdata.t(1:userdata.index);
            userdata.y = userdata.y(1:userdata.index,:);
            set(f, 'UserData', userdata);
            cols = size(userdata.y, 2);
            for j = 1:cols
                set(userdata.lines(j), 'Xdata', userdata.t,...
                                        'YData', userdata.y(:,j));
            end
            if ~ishold
                set(findobj(f,'Tag','stop'),'Visible','off'); %hide stop button
                set(ax, 'XLimMode', 'auto');
                refresh;
            end
        end
        
end
                
if userdata.stop     
    status = 1;
else
    status = 0;
end

end


function StopButtonCallback(~, ~)
ud = get(gcbf, 'UserData');
ud.stop = 1;
set(gcbf, 'UserData', ud);
end
