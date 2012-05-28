function [ t, x ] = trapezoidalMethod( stoich_const, pfcn_const, stoich_var, pfcn_var,...
                                       tspan, x0, rate_params, timeseries, delta,...
                                       output_fcn, MAX_OUTPUT_LENGTH)
%TRAPEZOIDALMETHOD Modified Gillespie algorithm to allow time-varying rate parameters
%   Based on: V. Shahrezaei, J.F. Ollivier, P.S. Swain (2008) Colored
%   extrinsic fluctuations and stochastic gene expression. Mol Sys Biol, 4,
%   196.
%
%   Usage:
%       [t, x] = 
%
%   Author:     Nezar Abdennur
%   Copyright:  (c) Nezar Abdennur 2012
%
%   See also 

if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end

%% Initialize
num_const = size(stoich_const, 1);
num_var = size(stoich_var, 1);
stoich_matrix = [stoich_const; stoich_var];
num_species = size(stoich_matrix,2);

T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)   = tspan(1); %repmat(tspan(1), 2,1);
X(1,:) = x0; %repmat(x0, 2,1);
rxn_count = 2;

ti = timeseries.times;
k = timeseries.values;
if ti(1) > tspan(1) || ti(end) < tspan(2)
    error('The input time series [%g,...,%g] must contain the simulation time span [%g, %g] as a subinterval.', ti(1), ti(end),tspan(1), tspan(2));
end
start = 1;
time = T(1);
Tbarrier = time + delta;
[k1, start] = lookup(k,ti, T(1), start);
[k2, start] = lookup(k,ti, Tbarrier, start);

%% MAIN LOOP
while time <= tspan(2)        
    % Calculate propensities with constant rate parameters
    Aconst = pfcn_const( X(rxn_count,:), rate_params );
    r1 = rand(num_const,1);
    tau_const = -log(r1)./Aconst;
   
    % Estimate propensities with varying rate parameters
    Avar1  = pfcn_var( time,     X(rxn_count,:), rate_params, k1 );
    Avar2  = pfcn_var( Tbarrier, X(rxn_count,:), rate_params, k2 );
    Aprime = (Avar2-Avar1)/(Tbarrier-time);
    r2 = rand(num_var,1);
    tau_var = (Avar1./Aprime) .* ( sqrt( 1-(2*Aprime./Avar1.^2).*log(r2) ) - 1 );
 
    is_const_rate = (Aprime==0);
    tau_var(is_const_rate) = (1./Avar1(is_const_rate)).*(log(1./r2(is_const_rate))); 
    
    is_complex = logical(imag(tau_var));
    tau_var(is_complex) = inf;   

    % Next reaction channel and time
    [tau, mu] = min([tau_const;tau_var]);
    
    % Update
    if ( time + tau <= Tbarrier )

        if rxn_count + 1 > MAX_OUTPUT_LENGTH
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:ExceededCapacity',...
                    'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
            return;
        end
        
        % jump to firing time and fire reaction.
        time = time + tau;
        T(rxn_count+1) = time;
        X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);
        rxn_count = rxn_count + 1; 
        [k1, start] = lookup(k,ti, time, start);
    else
        % jump to barrier time and stay quiescent. No need to log anything.
        time = Tbarrier;
%         T(rxn_count+1) = t;
%         X(rxn_count_1,:) = X(rxn_count,:);
        Tbarrier = Tbarrier + delta;
        [k1, start] = lookup(k,ti, time, start);
        [k2, start] = lookup(k,ti, Tbarrier, start);
    end               
                
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

% Record output
t = T(1:rxn_count);
x = X(1:rxn_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(rxn_count-1,:);
elseif t(end) < time
    t(end+1) = tspan(2);
    x(end+1,:) = X(rxn_count,:);
end
    

end


function [ki, start] = lookup(k, t, t_interp, start)
% LOOKUP Estimate function value K(t) at t=t_interp by linear interpolation
%   from the time series (t, k). Search starts at index start.
%
%   NOTE: Assumes there are no repeated time values in ti.

    try
        while t(start+1) < t_interp
            start = start + 1;
        end
        tprev = t(start); kprev = k(start,:);
        tnext = t(start+1); knext = k(start+1,:);
        ki = kprev + ( (knext-kprev)/(tnext-tprev) ) * (t_interp-tprev); 
    catch %#ok
        ki = k(end,:);
    end       
end
