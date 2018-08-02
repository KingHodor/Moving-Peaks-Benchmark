# Moving Peaks Benchmark

We carry out the following set of experiments on the MPB problems with different evolutionary algorithms. The corresponding algorithms are listed below:

- SGA (Simple Genetic Algorithm)
- HM (Hyper Mutation)
- RI (Rando Immigrants)
- MS (Memory Search)
- SOS (Self-Organizing Scouts)
- EIGA (Elitism-Based Immigrants)
- AMSO (Adaptive Multi-Swarm Optimizer)
- CPSOR (Clustering PSO with Restart)
- MIGA (GA with Memory-Based Immigrants)
- mQSO (The multi quantum swarm optimization)

We will explain these algorithms that we are using in our study in the following sections:

Simple Genetic Algorithm (SGA): One simple way to deal with DOPs is to regard the problem as a new one when a change occurs and restart GAs from scratch.

Hyper-mutation algorithm (HM): An enhanced version of the GA which consists in increasing the mutation rate whenever an environmental change is detected, in order to enhance the population diversity.

Random immigrants’ algorithm (RI): It is a modified GA trying to narrow down the problem of convergence, and this is done through a generational replacement of some individuals‟ randomly generated ones.

Memory Search algorithm (MS): The memory/search GA that combines the multi- population and memory schemes. In MS, in addition to the memory, MS maintains two populations P1 and P2 that evolve independently. When the memory is due to update, the best individual over P1 and P2 will replace the closest memory solution if it is fitter than the memory solution. The memory is reevaluated every generation.

Self-Organizing Scouts (SOS): In Self Organizing Scouts (SOS), the population splits a small fraction called the “scout population” watches over the peak, therest of the population called the “base population” spreads out and continues search for new peak.

Elitism-Based Immigrants (EIGA): The elitism-based immigrants‟ scheme combines the idea of elitism with the random immigrants‟ scheme. It differs from the aforementioned memory-based immigrants scheme in that the elite from the previous population instead of the best memory point is used to guide the immigrants toward the current environment.

Adaptive Multi-Swarm Optimizer(AMSO): The Adaptive Multi-Swarm Optimizerprovides a method to adaptively maintain the population diversity regarding tuning the number of populations for multi-population based algorithms without the assistance of change detection methods in dynamic environments.This algorithm uses the improved PSO with learning as a local search method for each population.

Clustering PSO with Restart (CPSOR): The CPSORapplies the random immigrants method without change detection based on a mechanism that can automatically reduce redundant individuals in the search space throughout the run.

GA with Memory-Based Immigrants (MIGA): The MIGA proposed a memory-based immigrants scheme for GAs in dynamic environments, where the memory is used to bias the immigrants toward the current environment.

The multiquantum swarm optimization (mQSO): In this paper, we explore new variants of particle swarm optimization (PSO) specifically designed to work well in dynamic environments. The main idea is to split the population of particles into a set of interacting swarms. These swarms interact locally by an exclusion parameter and globally through a new anti-convergence operator. In addition, each swarm maintains diversity either by using charged or quantum particles.

## Getting help

Please contact ahmet@evercoin.com

