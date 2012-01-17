classdef Test_2StateModel < SSARunner
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
    
    methods        
        function init(self)
            self.tspan = [0 10000];
            self.x0 = [0, 0];
            self.prop_fcn = @prop_2state;
            self.stoich_matrix = [ 1  0  ;
                                  -1  0  ;
                                   0  1  ;
                                   0 -1 ];
            self.param.kR = 0.1;%0.01;      
            self.param.kP = 0.1;%1;                     
            self.param.gR = 0.1;                        
            self.param.gP = 0.002;
            self.options = [];
            self.output = [];
            
            self.is_init = true;
        end
    end
    
end

%----------------------------------------------
% Propensity function
function a = prop_2state(x, p)
R = x(1);
P = x(2);

a(1) = p.kR;
a(2) = p.gR*R;
a(3) = p.kP*R;
a(4) = p.gP*P;
end