function varargout = ssa( stoich_matrix, propensity_fcn, tspan, x0, rate_params, varargin )
%SSA Runs a chemical kinetics simulation using the Gillespie Stochastic Simulation Algorithm
%   Usage:
%       [t, x] = ssa(stoich_matrix, prop_fcn, tspan, x0, rate_params=[])
%       [t, x] = ssa(stoich_matrix, prop_fcn, tspan, x0, rate_params=[], options)
%       [t, x] = ssa(stoich_matrix, prop_fcn, tspan, x0, rate_params=[], 'opt1', val1, 'opt2', val2, ...)
%       [t, x, te, xe, ie] = ssa(stoich_matrix, prop_fcn, tspan, x0, rate_params=[], ...)
%
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)
%
%       For located non-reaction events (see 'EventFcn' option):
%       te:             event times          
%       xe:             state at event times
%       ie:             event index          
%
%   Required:
%       tspan:          Initial and final times, [t_initial, t_final].
%       x0:             Initial species values, [S1_0, S2_0, ... ].
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%       prop_fcn:       Function handle to function that calculates
%                       reaction propensities.
%                       Target function should be of the form
%                           ac = f( xc, rate_params ),
%                       where xc is the current state [S1, S2, ...], and
%                       rate_params is the user-defined rate parameters.
%                       The function should return vector ac of 
%                       propensities (Nreactions x 1) in the same order as
%                       the reactions given in stoich_matrix.
%       
%   Optional:
%       rate_params:    Arbitrary user-defined data passed to prop_fcn.
%       options:        Simulation config options -- ssa accepts name/value
%                       pairs or an options struct. See options below.
%
%       Options
%       -------
%       MaxNumEvents:   Preallocated memory for simulation output (t,x). 
%                       (Default=1000000)
%
%       Method:         Simulation algorithm.
%                       ({'DIRECT'} | 'FIRST-REACTION')
%
%       OutputFcn:      Accepts an output function similar to the MATLAB 
%                       ode library. Use SSAPLOT or SSAPHAS2 to visualize
%                       the simulation trajectory in real time or
%                       SSAPROGBAR to display a graphical progress bar.
%                       (Default=none)
%
%       OutputSel:      Determines which species to display when SSAPLOT or
%                       SSAPHAS2 is used as the output function.
%                       (Default=1:Nspecies)
%
%       EventFcn:       Accepts an event function similar to the MATLAB ode
%                       library --- ssa will return the times, state and
%                       indices of any triggered events. See SSAEVENT for
%                       details.
%                       (Default=none)
%
%       Validate:       Whether to validate required input arguments.
%                       (Default=true)
%
%   See also SSA.directMethod, SSA.firstReactionMethod, SSAPLOT, SSAEVENT
%
%   Author:     Nezar Abdennur
%   Copyright:  (c) Nezar Abdennur 2012
%   Revision:   12.01.19

% If no input, show usage
if nargin == 0
    help(mfilename);
    return;
end

% Check input args
if nargin < 4
    error('SSA:InputError',...
          'Not enough input arguments. Type help(''ssa'').');
end

if ~exist('rate_params','var')
    rate_params = []; %no rate parameters provided
end

% Optional args accepted as a struct or name-value pairs
if ~isempty(varargin)
    if isstruct(varargin{1})
        varargin = reshape([fieldnames(varargin{1}),...
                            struct2cell(varargin{1})]', 1,[]);
    elseif rem(length(varargin),2)~=0
        error('SSA:InputError:NameValueMismatch',...
              'Optional arguments must be given as name-value pairs.');
    end
end

% A struct to store optional input args
opt_names = { 'Method',...
              'OutputFcn',...
              'OutputSel',...
              'EventFcn',...
              'MaxNumEvents',...
              'Validate' };
opt_init = [opt_names; {[], [], [], [], [], []}];
options = struct(opt_init{:});
options.Validate = true; % validate input by default

% Parse optional input args (name-value pairs)
for i = 1:2:length(varargin)
    opt_name = validatestring(varargin{i}, opt_names);
    options.(opt_name) = varargin{i+1};
end

% Validate input model
if options.Validate
    assert(isnumeric(stoich_matrix) && ~isempty(stoich_matrix),...
           'SSA:InputTypeError',...
           'Stoichiometry matrix provided is invalid.');
    num_rxns = size(stoich_matrix, 1);
    num_species = size(stoich_matrix, 2);
    
    assert(isnumeric(x0) && isvector(x0),...
           'SSA:InputTypeError',...
           'Initial value vector provided is invalid');
    assert(length(x0) == num_species,...
           'SSA:InputSizeError',...
           'Size of the initial value vector is inconsistent with the stoichiometry matrix.');
    x0 = x0(:)';
    
    assert(isnumeric(tspan) && length(tspan) == 2 && issorted(tspan),...
           'SSA:InputError',...
           'Invalid tspan vector specified');
    assert(isa(propensity_fcn, 'function_handle'),...
           'SSA:InputTypeError',...
           'propensity_fcn must be a function handle.');
    assert(nargin(propensity_fcn) >= 2 && nargout(propensity_fcn) > 0,...
           'SSA:InputError',...
           'propensity_fcn must take at least 2 input arguments and 1 output argument.');
    try
        a_test = propensity_fcn( zeros(1, num_species) , rate_params);
    catch err
        error('SSA:InputError',...
              'Propensity function failed to evaluate:\n%s', err.message);
    end
    assert(isnumeric(a_test) && isvector(a_test),...
           'SSA:TypeError',...
           'Propensity function must return a numeric vector.');
    assert(length(a_test) == num_rxns,...
           'SSA:SizeError',...
           'Size of vector returned by propensity_fcn is inconsistent with the stoichiometry matrix.'); 
end

% Set user-defined options / default values
if isempty(options.MaxNumEvents)
    max_num_events = 1000000;
elseif options.MaxNumEvents > 0
    max_num_events = round(options.MaxNumEvents);
else
    error('SSA:OptionsError',...
          'MaxNumEvents option must be a positive integer.');
end

if isempty(options.Method)
    update_method = 'DIRECT';
else
    update_method = validatestring(options.Method, {'DIRECT','FIRST-REACTION'});
end

if isempty(options.OutputSel)
    options.OutputSel = 1:num_species;
else
    num_selected = length(options.OutputSel);
    for i = 1:num_selected
        if ~ismember(options.OutputSel(i), 1:num_species)
            error('SSA:OptionsError',...
                  'OutputSel indices must be in range 1 to Nspecies.');
        end
    end
end

special_fcns = {'ssaprogbar','ssaplot','ssaphas2'};
if isempty(options.OutputFcn) && isempty(options.EventFcn)
    output_fcn = [];
    
elseif ~isempty(options.EventFcn)
    if ischar(options.EventFcn)
        options.EventFcn = str2func(options.EventFcn);
    end
    if isa(options.EventFcn, 'function_handle')
        ssaevent('init', tspan, x0(:), @options.EventFcn);
        output_fcn = @(ti,xi) ssaevent('update', ti,xi, @options.EventFcn);
        c = onCleanup( @() ssaevent('done') );   
    else
        error('SSA:OptionsError',...
              'EventFcn option must be a valid function handle.'); 
    end
    
elseif ~isempty(options.OutputFcn)
    if ischar(options.OutputFcn)
        try options.OutputFcn = validatestring(options.OutputFcn, special_fcns);
        catch exc; end
        options.OutputFcn = str2func(options.OutputFcn);
    end      
    if isa(options.OutputFcn, 'function_handle')
        if isequal(options.OutputFcn, @ssaplot)
            ssaplot('init', tspan, x0(:), options.OutputSel);
            output_fcn = @(ti,xi) ssaplot('update',ti,xi);
            c = onCleanup( @() ssaplot('done') );
        elseif isequal(options.OutputFcn, @ssaphas2)
            ssaphas2('init', tspan, x0(:), options.OutputSel);
            output_fcn = @(ti,xi) ssaphas2('update',ti,xi);
            c = onCleanup( @() ssaphas2('done') );
        elseif isequal(options.OutputFcn, @ssaprogbar)
            ssaprogbar('init', tspan);
            output_fcn = @(ti,xi) ssaprogbar('update',ti,xi);
            c = onCleanup( @() ssaprogbar('done') );
        else
            output_fcn = options.OutputFcn;
        end
    else
        error('SSA:OptionsError',...
               'OutputFcn option must be a valid function handle.');
    end
end
      
% Run simulation algorithm
switch update_method
    case 'DIRECT'
        [ t, x ] = directMethod(stoich_matrix, propensity_fcn, tspan, x0,...
                                rate_params, output_fcn, max_num_events);
    case 'FIRST-REACTION'
        [ t, x ] = firstReactionMethod(stoich_matrix, propensity_fcn, tspan, x0,...
                                       rate_params, output_fcn, max_num_events);
end

% Output
if ~isempty(options.EventFcn) && nargout == 5
    [~,te,xe,ie] = ssaevent('get');
    varargout = {t,x,te,xe,ie};
else
    varargout = {t,x};
end

end

