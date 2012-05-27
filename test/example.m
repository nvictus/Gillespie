% Two-state model of gene expression
%
% Reactions:
%   0 -> mRNA
%   mRNA -> mRNA + protein
%   mRNA -> 0
%   protein -> 0

tspan = [0 10000];
x0 = [0, 0];
stoich_matrix = [ 1  0  ;
                 -1  0  ;
                  0  1  ;
                  0 -1 ];
p.kR = 0.1;%0.01;      
p.kP = 0.1;%1;                     
p.gR = 0.1;                        
p.gP = 0.002;

[t,x] = ssa(stoich_matrix, @prop_2state, tspan, x0, p);
figure(1)
stairs(t,x);
