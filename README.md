Gillespie
=========

Gillespie Stochastic Simulation Algorithm 

The two classic versions of the algorithm implemented in MATLAB:
- The _direct_ method
- The _first-reaction_ method

Example model
-------------
```
Reaction network:
    1. transcription:       0       --kR--> mRNA
    2. translation:         mRNA    --kP--> mRNA + protein
    3. mRNA decay:          mRNA    --gR--> 0
    4. protein decay:       protein --gP--> 0
```
```matlab
>> ssa_example
```
![Simulation output](ssa.png)


