% gillespie
%
% Author:     Nezar Abdennur <nabdennur@gmail.com>
% Created:    2012-01-07
% Copyright:  (c) Nezar Abdennur 2012
% Release:    2.0.0
%
% Files
%   ssa        - Runs a chemical kinetics simulation using the Gillespie Stochastic Simulation Algorithm
%   ssaevent   - Used to determine the location of user-defined events in an ssa simulation
%   ssaphas2   - Used to update a 2D phase diagram at every time step during an ssa simulation
%   ssaplot    - Used to update a time course plot at every time step during an ssa simulation
%   ssaprogbar - Used to display a progress bar window during an ssa simulation
%
%   directMethod        - Implementation of the Direct Method variant of the Gillespie algorithm
%   firstReactionMethod - Implementation of the First-Reaction Method variant of the Gillespie algorithm
%   trapezoidalMethod   - Modified Gillespie algorithm to allow time-varying rate parameters
