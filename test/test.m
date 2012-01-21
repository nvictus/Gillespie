%% Create test object to run simulations
t = Test_2StateModel();

%%
fprintf('Testing two-state model, direct method...\n');
t.init();
t.run_ssa();
t.plot_output('Direct');
pause;

%%
clf;
fprintf('Testing two-state model, first-reaction method...\n');
t.init();
t.run_ssa();
t.plot_output('First-Reaction'); 
pause;

%%
clf;
fprintf('Change parameters...\n');
t.init();
t.param.kR = 1;
t.param.gR = 0.05;
t.run_ssa();
t.plot_output();
pause;

%%
clf;
fprintf('Test ssaplot -- real-time plot...\n');
t.init();
t.options.OutputFcn = @ssaplot;
t.run_ssa();
t.plot_output();
pause;

%%
clf;
fprintf('Test ssaphas2 -- real-time phase diagram...\n');
t.init();
t.options.OutputFcn = @ssaphas2;
t.run_ssa();
t.plot_output();
pause;

%%
clf;
fprintf('Test ssaprogbar -- slow-ass waitbar...\n');
t.init();
t.options.OutputFcn = @ssaprogbar;
t.run_ssa();
t.plot_output();
pause;

%%
clf;
fprintf('Test using an event function -- see my_event.m...\n');
t.init();
t.options.EventFcn = @my_event;
t.run_ssa();
t.plot_output('With events!');
t.plot_events();
pause;

