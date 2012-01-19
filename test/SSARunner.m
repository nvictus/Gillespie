classdef SSARunner < handle
    %SSARUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        stoich_matrix
        prop_fcn
        param
        x0
        tspan
        options
        output
    end
    
    methods (Abstract)
        init(self)
    end
    
    
    
    properties (Hidden)
        is_init = false;
        hfigure = []
    end
    
    methods
        function [t,x] = run_ssa(self)
            if ~self.is_init
                error('Model not initialized! Call init method.');
            end
            
            self.output.t = [];
            self.output.x = [];
            self.output.te = [];
            self.output.xe = [];
            self.output.ie = [];
            if isempty(self.options)
                [t,x] = ssa(self.stoich_matrix, self.prop_fcn, self.tspan,...
                            self.x0, self.param);
            elseif isfield(self.options, 'EventFcn')
                [t,x,te,xe,ie] = ssa(self.stoich_matrix, self.prop_fcn, self.tspan,...
                                     self.x0, self.param, self.options);
                self.output.te = te;
                self.output.xe = xe;
                self.output.ie = ie;
            else
                [t,x] = ssa(self.stoich_matrix, self.prop_fcn, self.tspan,...
                            self.x0, self.param, self.options);
            end
            self.output.t = t;
            self.output.x = x;
        end
        
        function plot_output(self, name)
            if isempty(self.hfigure) || ~ishghandle(self.hfigure)
                self.hfigure = figure('Name', 'ssaplot');
            end
            figure(self.hfigure);
            hold on;
            num_species = size(self.output.x, 2);
            c = lines(num_species);
            for i = 1:num_species
                stairs(self.output.t,self.output.x(:,i),'Color', c(i,:));
            end
            hold off;
            if exist('name', 'var')
                title(name)
            end
        end   
        
        function plot_events(self)
            if ishghandle(self.hfigure) && ~isempty(self.output.ie)
                hold on;
                num_ie = max(self.output.ie);
                c = jet(num_ie);
                for i = 1:num_ie
                    idx = self.output.ie == i;
                    te = self.output.te(idx);
                    xe = self.output.xe(idx,:);
                    plot(te, xe, '*', 'color', c(i,:));
                end
                hold off;
            end
        end
    end    
end
