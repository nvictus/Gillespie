%%
fprintf('Testing two-state model...\n');
t = Test_2StateModel();

t.init();
fprintf('    1) no options...');
t.run_ssa();
fprintf('Y\n');

%%
t.init();
t.options.Method = 'FIRST';
fprintf('    2) Method = ''FIRST''...');
t.run_ssa();
t.plot_output('Test 2'); 
pause;
fprintf('Y\n');

%%
t.init();
t.param.kR = 0.01;
t.param.gR = 1;
fprintf('    3) Change rate params...');
t.run_ssa();
t.plot_output('Test 3');
pause;
fprintf('Y\n');

%%
t.init();
t.options.OutputFcn = @ssaplot;
fprintf('    4) OutputFcn = ''ssaplot''...');
t.run_ssa();
t.plot_output('Test 4');
pause;
fprintf('Y\n');

%%
t.init();
t.options.OutputFcn = 'ProgressBar';
fprintf('    5) OutputFcn = ''ProgressBar''...');
t.run_ssa();
t.plot_output('Test 5');
pause;   
fprintf('Y\n');

