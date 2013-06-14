function [ t, x ] = firstReactionMethod( stoich_matrix, propensity_fcn, tspan, x0,...
                                         params, output_fcn, MAX_OUTPUT_LENGTH)
%FIRSTREACTIONMETHOD Implementation of the First-Reaction Method variant of the Gillespie algorithm
%   Usage:
%       [t, x] = firstReactionMethod( stoich_matrix, propensity_fcn, tspan, x0 )
%       [t, x] = firstReactionMethod( stoich_matrix, propensity_fcn, tspan, x0, params )
%       [t, x] = firstReactionMethod( stoich_matrix, propensity_fcn, tspan, x0, params, output_fcn )
%
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)    
%
%   Required:
%       tspan:          Initial and final times, [t_init, t_final].
%
%       x0:             Initial species amounts, [S1_0, S2_0, ... ].
%
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%
%       prop_fcn:       Function that calculates reaction propensities.
%                       Function should be of the form
%                           ac = f( xc, params )
%                       where xc is the current state [S1, S2, ...], and
%                       params is the user-defined rate parameters.
%                       The function should return column vector ac of 
%                       propensities (Nreactions x 1) in the same order as
%                       the reactions given in stoich_matrix.
%
%   Optional:
%       params:         User-defined parameters, passed to prop_fun (e.g.,
%                       a struct of rate constants) <default=[]>
%
%       output_fcn:     Arbitrary function with signature
%                           status = f( tc, xc )
%                       The output_fcn is passed the current time and state
%                       after each step of the simulation. It can be used
%                       to locate events, monitor progress, write data,
%                       etc. If it returns 1, the simulation terminates.
%                       <default=none>
%
%   Reference: 
%       Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%       Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%
%   See also DIRECTMETHOD

%   Nezar Abdennur, 2012 <nabdennur@gmail.com>
%   Dynamical Systems Biology Laboratory, University of Ottawa
%   www.sysbiolab.uottawa.ca
%   Created: 2012-01-19

if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end
if ~exist('params', 'var')
    params = [];
end

% Initialize
num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;

% MAIN LOOP
while T(rxn_count) < tspan(2)        
    % Calculate reaction propensities
    a  = propensity_fcn(X(rxn_count,:), params);
    
    % Sample time-to-fire for each reaction channel
    % tau is the smallest time-to-fire and mu is the index of the channel
    r = rand(num_rxns,1);
    taus = -log(r)./a;
    [tau, mu] = min(taus);
    
    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count,:);
        warning('SSA:ExceededCapacity',...
                'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
        return;
    end
    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);   
    rxn_count = rxn_count + 1; 
    
    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count,:)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end
end  

% Return simulation time course
t = T(1:rxn_count);
x = X(1:rxn_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(rxn_count-1,:);
end    

end

