implicit
========

Computational Physics; Solving Differential Equations using implicit methods for stability

4 methods used to solve problems:
Forward Euler
Midpoint
Backward Euler
Crank-Nicholson

Two problems solved:

1) Oscillator. Just a simple damped oscillator problem. Can see the differences in results from different schemes

2) Reaction rates. A system where elements decay into another, not exactly difficult. Interesting part is the use of 
variable time steps; was able to resolve these reactions fairly well with a low resolution due to adaptive step size! 
Can be improved--Paul resolved reaction in ~80 steps, which is half what I used...
