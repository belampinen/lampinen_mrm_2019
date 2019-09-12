# The stochastic Self Organizing Migrating Algorithm (SOMA)

The stochastic Self Organizing Migrating Algorithm (SOMA) generates a ‘population’ of ‘individuals’ (candidate solutions) that it proceeds to evolve in a number of ‘migrations’. Each individual is represented by a scalar vector (**v**) and is initialized through random selection from uniform distributions within boundaries specified by the problem. 

Evolution of individuals through migrations consists of the following steps:
1) Assign a ‘distance’ (*d*) to each individual using a metric (the objective function) and identify the individual with the lowest *d* (the ‘leader’; **v**L).
2) Mutate non-leader individuals by adjusting randomly chosen vector elements to become more similar to those of the leader.
3) Repeat from step 1 unless maximal number of iterations have been performed or the improvement in the leader’s d is below a threshold value. 

Step 2, in more detail, generates a set of mutated versions {**v**m} for each individual and selects the version with the lowest *d*. The set is generated according to

**v**m = **v** + *m* &sdot; *sl* &sdot; (**v**L – **v**) &cir; **x**, 
    
where *sl* is the ‘step length’ and &cir; denotes element-wise multiplication with a binary vector **x** created by randomly setting each element to unity with a fixed probability *p*. By default, *m* ranges from 0 (returning **v**) to 10 and *sl* = 0.21, creating 11 versions with elements mutated along a ‘path’ from the current position (**v**), crossing the leader position (**v**L), to a position on the near ‘opposite side’. Values created outside the allowed boundaries are replaced by randomized values. The random elements together with the ability of individuals to explore ‘jumps’ over the current best solution makes the SOMA algorithm robust to local minima. In our implementation, SOMA was executed using default settings but with 200 migrations and a step likelihood *p* = 0.2.


# Reference
Zelinka, I. (2004). SOMA—self-organizing migrating algorithm. New optimization techniques in engineering, Springer: 167-217.