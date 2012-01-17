classdef Test_RandomModel < SSARunner
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stoich_matrix
        prop_fcn
        param
        x0
        tspan
        options
        output
    end
    
    properties
        Nspecies
        Nreactions
    end
    
    methods
        function self = TestSSA_RandomModel(Nspecies, Nreactions)
            %self.stoich_matrix =
            % for each column, sample gaussian distrib.
            % within 1 std, stoich:=0
            % between 1 and 2 std, stoich:=-1 or 1
            % greater than 2 stf, stoich:=-2 or 2
            % make sure at least one nonzero coeff per column
            
            % set propensities...
            % for each row, check for any negative stoich coefficients
            % negative coefficients should be accompanied by a propensity
            % which is at least first-order wrt to that species (to prevent
            % negative species quantities)
            
            %other possible approach:
            % make a list of feasible reaction types:
            % pseudo zero-order:
            % - birth         (0 -> A) *
            % - death         (A -> 0 if nA~=0) not very plausible
            %
            % unimolecular:
            % - isomerization/conversion (A -> B)
            % - decay/death/deg          (A -> 0)
            % - fragmentation            (A -> B + C) *
            % - catalytic synthesis      (A -> A + B) *
            % - autocatalytic synthesis  (A -> A + A) *
            %
            % bimolecular:
            % - complexation             (A + B -> C)
            % - catalytic degradation    (A + B -> A)
            % - annihilation             (A + B -> 0)
            % - dimerization             (A + A -> B)
            % - autocatalytic degrad     (A + A -> A)
            % - autoannihilation         (A + A -> 0)
            % - catalytic conversion     (A + B -> A + C)
            % - autocatalytic conversion (A + A -> A + B)
            % - mutual conversion        (A + B -> C + D) * 
            % - cat fragmentation        (A + B -> A + C + D) *
            % - autocat fragmentatn      (A + A -> A + B + C) *
            %        
            %
            % * possible to add an extra product to each rxn (zero-order
            % source)
        end
    end
    
end

