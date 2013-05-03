classdef SsaModel < handle
    %UNTITLED16 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        species   = struct()  %struct:name->amount
        reactions = {}        %cell:#->reaction , where reaction = k, {n1,n2}, {struct:name->S}
        time      = 0         %double
    end
    
    properties
        stoich_matrix
        propensity_fcn
        params
        x0
    end
    
    methods
        function model = SsaModel()
        end
        
        function addSpecies(self, name, amount)
            self.species.(name) = amount;
        end
        
        function addReaction(self, rate_constant, reactants, output)
            stoich = struct();
            for i = 1:2:length(output)
                stoich.(output{i}) = output{i+1};
            end
            num_rxns = length(self.reactions);
            self.reactions{num_rxns+1}.constant = rate_constant;
            self.reactions{num_rxns+1}.reactants = reactants;
            self.reactions{num_rxns+1}.stoichiometry = stoich;
        end
        
        function setTime(self, t)
            self.time = t;
        end
        
        function setAmount(self, name, amount)
            self.species{name} = amount;
        end
        
        function generate_model(self)
            % Make species initial amounts vector
            species_names = fieldnames(self.species);
            num_species = length(species_names);
            num_rxns = length(self.reactions);
            cSpecies = struct2cell(self.species);
            for i = 1:num_species
                self.x0(i) = cSpecies{i};
            end
            
            % Stamp the stoichiometry matrix
            S = zeros(num_rxns, num_species);
            for i = 1:num_rxns
                rxn = self.reactions{i};
                output_names = fieldnames(rxn.stoichiometry);
                for ii = 1:length(output_names)
                    output_name = output_names{ii};
                    j = strcmp(output_name, species_names);
                    S(i,j) = rxn.stoichiometry.(output_name);
                end
            end
            self.stoich_matrix = S;
                
            % Make propensity function
            % Reactions should be at most biomolecular
            for i = 1:num_rxns
                rxn = self.reactions{i};
                p.k(i) = rxn.constant;
                input_names = rxn.reactants;
                try
                    p.idx{i} = find(strcmp(input_names{1}, species_names));
                catch err
                    try
                        p.idx{i}(2) = find(strcmp(input_names{2}, species_names));
                        if p.idx{i}(1) == p.idx{i}(2)
                            p.idx{i}(2) = 0;
                        end
                    catch err
                    end
                end
            end          
            self.propensity_fcn = @pfun;
            self.params = p;
        end
        
    end
    
end

function a = pfun(x,p)
    nr = length(p.idx);
    a = zeros(nr,1);
    for i = 1:nr
        a(i) = p.k(i);
        for j = 1:length(p.idx{i})
            if j == 0
                a(i) = a(i)*(x(p.idx{i}(j)) - 1)/2;
            end
            a(i) = a(i)*x(p.idx{i}(j));
        end
    end
end
