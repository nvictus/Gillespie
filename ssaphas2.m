function status = ssaphas2( flag, t, y, indices )
%SSAPHAS2 Used to update a 2D phase diagram at every time step during an ssa simulation
%   Used by ssa function when the function ssaphas2 is passed to ssa as the
%   'OutputFcn' property. To choose which two species to plot, pass in
%   their indices as the 'OutputSel' property. ssaphas2 works similarly to
%   the odephas2 function used with the MATLAB ODE solvers and is based on
%   that function.
%
%   See also SSA, ODEPHAS2

persistent HFIGURE HAXES

chunk = 128;

switch flag
    
    case 'update'
        if ( isempty(HFIGURE) || isempty(HAXES) )
            error('SSAPHAS2:NoInit',...
                  'SSAPHAS2 not initialized.');
        elseif ( ishghandle(HFIGURE) && ishghandle(HAXES) )

            try
                userdata = get(HFIGURE, 'UserData');
                
                % Append new data
                y = y(userdata.select);
                old_i = userdata.index;
                new_i = userdata.index + 1;
                nrows = size(userdata.y, 1);
                if new_i > nrows
                    userdata.y = [userdata.y; zeros(chunk,2)];
                end
                
                userdata.y(old_i+1:new_i,:) = y(userdata.select).';
                userdata.index = new_i;
                set(HFIGURE, 'UserData', userdata);
                
                % Plot new data
                ylim = get(HAXES, 'YLim');
                xlim = get(HAXES, 'XLim');
                if ( old_i == 1 || ...
                       ( y(1) < xlim(1) || xlim(2) < y(1) ) || ...
                       ( y(2) < ylim(1) || ylim(2) < y(2) ) );
                    %axes out of range or just initialized -> replot everything
                    set(userdata.lines, 'XData', userdata.y(1:new_i, 1),...
                                        'YData', userdata.y(1:new_i, 2));
                    set(userdata.line,  'XData', userdata.y(old_i:new_i, 1),...
                                        'YData', userdata.y(old_i:new_i, 2));
                else
                    %plot only new data
                    color_order = get(HAXES, 'ColorOrder');
                    %unhighlight old leading segment
                    set(userdata.line, 'Color', color_order(1,:));
                    %add and highlight new leading segment
                    set(userdata.line, 'XData', userdata.y(old_i:new_i, 1),... 
                                       'YData', userdata.y(old_i:new_i, 2),...
                                       'Color', color_order(2,:));
%                     set(userdata.marker, 'XData', userdata.y(new_i, 1),...
%                                          'YData', userdata.y(new_i, 2));
                end
                drawnow;

            catch err
                error('SSAPLOT:ErrorUpdatingWindow',...
                      'Error updating plot window: %s', err.message);
            end
        end

    case 'init'
        nselected = length(indices); 
        if nselected < 2
            error('SSAPHAS2:TooFewSpeciesSelected',...
                  'Two species must be selected for odephas2 phase diagram.');
        end
        
        f = figure(gcf);
        HFIGURE = f;
        HAXES = gca;
        
        userdata.select = indices;
        userdata.y = zeros(chunk, nselected);
        userdata.index = 1;
        userdata.y(1,:) = y(userdata.select).';
        
        color_order = get(HAXES, 'ColorOrder');
        y_init = userdata.y(1,:);
        if ~ishold
            userdata.lines = plot(y_init(1), y_init(2), '-',...
                                                        'Marker', '.');
            hold on;
            userdata.line = plot(y_init(1), y_init(2), '-', ...
                                                       'Marker', '.',...
                                                       'Color', color_order(2,:),...
                                                       'EraseMode', 'none');
%             userdata.marker = plot(y_init(1), y_init(2), 'Marker', 'o',...
%                                                          'Color', color_order(2,:),...
%                                                          'EraseMode', 'background');                                       
            hold off
            
        else
            userdata.lines = plot(y_init(1), y_init(2), '-',...
                                                        'Marker', '.');
            userdata.line = plot(y_init(1), y_init(2), '-', ...
                                                       'Marker', '.',...
                                                       'Color', color_order(2,:),...
                                                       'EraseMode', 'none');
%             userdata.marker = plot(y_init(1), y_init(2), 'Marker', 'o',...
%                                                          'Color', color_order(2,:));
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
        HFIGURE = [];
        HAXES = [];
        
        if ishghandle(f)
            userdata = get(f, 'UserData');
            userdata.y = userdata.y(1:userdata.index,:);
            set(f, 'UserData', userdata);
            set(userdata.lines, 'Xdata', userdata.y(:,1),...
                                'YData', userdata.y(:,2));
            if ~ishold
                set(findobj(f,'Tag','stop'),'Visible','off'); %hide stop button
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
