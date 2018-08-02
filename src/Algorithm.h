/*
 * Algorithm.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */


#ifndef EVOLUTIONARYLIB_ALGORITHM_H
#define EVOLUTIONARYLIB_ALGORITHM_H


#include <glob.h>
#include "Population.h"
#include "ScoutPopulation.h"
#include "ParentPopulation.h"
#include "Statistics.h"
#include "PopPSO.h"
#include "Swarm.h"
#include "SwarmManager.h"


#define FITNESS_EVALUATION

class Algorithm {
public:
	static unsigned int TotalTime;
	static unsigned int ChangeTime;
	static int FitnessEvaluationSize;
    static unsigned int PopulationSize;
    static unsigned int GenerationCount;
    static unsigned int ChangePeriod;
    static double CrossoverProbability;
    static double MutationProbability;
    static unsigned int TournamentSize;
    static bool Elitism;
    static unsigned int RunCount;
    static Statistics statistics;
    static unsigned int NumberOfFitnessEvaluation;
    static unsigned int FinessEvaluationLimit;

    //EIGA
    static double RatioOfElitismBasedImmigrants;

    // Random Immigrants parameters
    static double RandomImmigrantsChangePercentage;

    // Memory Search parameters
    static unsigned int ExplicitMemorySize;
    static unsigned int MemoryUpdateFreq;

    // Hyper mutation parameters
    static double HyperMutationProbability;

    // Moving Peaks Parameters
    static unsigned int NumberOfPeaks;
    static double MinCoordinate;
    static double MaxCoordinate;
    static unsigned int DimensionSize;
    static double ShiftLength;


    // SOS parameters
    static double SOS_Lambda;
    static bool   SOS_UseBasisFunction;
    static unsigned int SOS_ForkingGenerationPeriod;
    static double SOS_MinScoutPopulationSizeRelative;
    static double SOS_MaxScoutPopulationSizeRelative;
    static double SOS_MinBasePopulationSizeRelative;
    static double SOS_MinDiameterRelative;
    static double SOS_MaxDiameterRelative;
    static double SOS_MinFitnessOfNewForkingPopulationsRelative;
    static double SOS_MinFitnessOfExistingForkingPopulationsRelative;
    static double SOS_DiameterReduceFactor;
    static double SOS_ALpha;

    //EIGA parametres
    static double EigaImmigrantsRatio;

    //CPSOR

    //MIGA parameters
    static double RandomImmigrantsRatioForMIGA;
    static double MutationRatioOfBestOrEliteIndividual;

    //MQSO parameters
	static unsigned int NumberOfSwarms;
	static unsigned int NumberOfNeutralParticlesInEachSwarm;
	static unsigned int NumberOfQuantumParticlesInEachSwarm;
	static double ControlAttractionOfGlobalBest;
	static double ControlAttractionOfPersonalBest;
	static double ConvergenceRadius;
	static double ExclusionRadius;
	static double ExclusionRadiusSquare;
	static double QuantumCloudRadius;

    static void InitMovPeaks();
    static void ReleaseMovPeaks();


    static void SelectParents( Population* population, int* momIndex, int* dadIndex, size_t tournamentSize);
    static std::vector<Individual> SimulatedBinaryCrossover( const Individual& mom, const Individual& dad );
    static std::vector<Individual> HillClimbingCrossover( Individual& mom, Individual& dad );
    static void GenerationalGA( Population* population );
    static void GenerationalGA( Population* population, size_t suggestedPopulationSize);

    static void ApplyElitism( Population* population, const Individual& best );

    // simple ga
    static void RunSimpleGA();

    // hyper mutation
    static void RunHyperMutation();

    // random immigranats
    static void RunRandomImmigrants();


    // memory search
    static void RunMemorySearch();
    static void UpdateMemoryPopulation(Population &memoryPopulation, vector<Individual> &explicitMemory);
    static void UpdateExplicitMemory( Population& memoryPopulation,
                                      Population& searchPopulation, std::vector<Individual>& explicitMemory );
    static Individual* MinDistIndividual( std::vector<Individual>& explicitMemory, Individual* bestIndividual );

    // self-organizing scouts
    static void RunSelfOrganizingScouts();
    static void SOSGenerationalGA( ParentPopulation* parentPopulation );
    static void SOSAdjustPopulationSizes(ParentPopulation *parentPopulation);
    static void SOSAdjustSearchSpace( ParentPopulation* parentPopulation );
    static void SOSDeleteNonFitScouts( ParentPopulation* parentPopulation );
    static double SOSGetRelativeDynamism(Population *population, double totalDynamism);
    static double SOSGetRelativeFitness(Population *population, double overallMinFitness, double totalFitness);
    static double SOSGetRelativeQuality(Population *population,
                                        double beta, double totalDynamism, double overallMinFitness,
                                        double totalFitness);
    static double SOSGetQuality(Population *population,
                                double beta, double totalDynamism, double overallMinFitness,
                                double totalFitness);
    static size_t SOSSuggestedSize(Population *population, double quality);
    static ScoutPopulation SOSFork(ParentPopulation *parentPopulation);
    static bool SOSMergeScoutPopulations( ParentPopulation *parentPopulation );
    static void MergeScoutPopulations( ParentPopulation* parentPopulation, ScoutPopulation* better, ScoutPopulation* worse );

    static void SOSReduceScoutDiameters( ParentPopulation* parentPopulation );

    //elitism-based immigrants GA
   static void RunEIGA();
   static void InitializeElitismImmigrants(vector<Individual> &elitismImmigrants, size_t elitismImmigrantsSize);
   static void MutateElite(Population* population, vector<Individual> &elitismImmigrants, Individual *base);
   static void ApplyEIGAElitism( Population* population, const Individual& best );

    // CPSOR
    static void RunCPSOR();
    static void ClusteringCPSOR(PopPSO *population, vector<PopPSO> *plst);
    static double** calculatemMatrix(std::vector<PopPSO>& clusters);
    static void deletemMatrix(double **matrix, int size);
    static double calculateDistance(IndPSO& i, IndPSO& j);
    static bool findNearestPairs(std::vector<PopPSO>& clusters, int *t, int *s, double **mMatrix);
    static void mergeClusters(PopPSO *t, PopPSO *s);
    static void removePLST(std::vector<PopPSO>& clusters);
    static double ratioOverlap(PopPSO *t, PopPSO *s);
    static double findBestFitnessforCPSOR(vector<PopPSO> &plst);

    //An Adaptive Multi-Swarm Optimizer
   static int getNextIndis(int preIndis, int curPops, int prePops, int *counter);
   static void calculateDistanceMatrix(vector<Swarm> &GList);
   static void mergeClusters(Swarm *tCluster, Swarm *sCluster);
   static bool findNearestPair(vector<Swarm> &GList, int *tIndex, int *sIndex);
   static void mergeAndRemoveOverlappedClusters(vector<Swarm> &plst, Swarm *clst);
  static double ratioOverlap(Swarm *t, Swarm *s);
   static void overlappingDetection(vector<Swarm> &, Swarm *clst);
   static void allocateDistanceMatrix(size_t gSize);
   static void deallocateDistanceMatrix( size_t gSize);
   static void clustering(Swarm *pop, vector<Swarm> *plst);
   static void RunAMSO();;
   static double findBestFitnessforAMSO(vector<Swarm> &plst);

   // Memory Based Immigrants GeneticAlgorithm
   static void RunMemoryBasedImmigrantsGeneticAlgorithm();

   // Multi Quantum Swarm Algorithm
   static void RunMultiQuantumSwarmAlgorithm();

};


#endif //EVOLUTIONARYLIB_ALGORITHM_H
