/*
 * Algorithm.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "Algorithm.h"
#include "movpeaks.h"
#include "Util.h"
#include "SwarmManager.h"

using namespace std;

#define subSize  7
#define E_ClusterRadiusThreshold 0.01
#define B_OverlabThreshold 0.1
#define A_DiversityThreshold 0.3

unsigned int Algorithm::TotalTime;
unsigned int Algorithm::ChangeTime;
unsigned int Algorithm::PopulationSize;
unsigned int Algorithm::GenerationCount;
unsigned int Algorithm::ChangePeriod;
double Algorithm::CrossoverProbability;
double Algorithm::MutationProbability;
unsigned int Algorithm::TournamentSize;
bool Algorithm::Elitism = true;
unsigned int Algorithm::RunCount;
Statistics Algorithm::statistics;
unsigned int Algorithm::NumberOfFitnessEvaluation;
unsigned int Algorithm::FinessEvaluationLimit;

// Random Immigrants parameters
double Algorithm::RandomImmigrantsChangePercentage;

// Memory Search parameters
unsigned int Algorithm::ExplicitMemorySize;
unsigned int Algorithm::MemoryUpdateFreq;

// Hyper mutation parameters
double Algorithm::HyperMutationProbability;

// Moving Peaks Parameters
unsigned int Algorithm::NumberOfPeaks;
double Algorithm::MinCoordinate;
double Algorithm::MaxCoordinate;
unsigned int Algorithm::DimensionSize;
double Algorithm::ShiftLength;

// SOS parameters
double Algorithm::SOS_Lambda;
bool   Algorithm::SOS_UseBasisFunction;
unsigned int Algorithm::SOS_ForkingGenerationPeriod;
double Algorithm::SOS_MinScoutPopulationSizeRelative;
double Algorithm::SOS_MaxScoutPopulationSizeRelative;
double Algorithm::SOS_MinBasePopulationSizeRelative;
double Algorithm::SOS_MinDiameterRelative;
double Algorithm::SOS_MaxDiameterRelative;
double Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative;
double Algorithm::SOS_MinFitnessOfExistingForkingPopulationsRelative;
double Algorithm::SOS_DiameterReduceFactor;
double Algorithm::SOS_ALpha;

//EIGA parameters
double Algorithm::EigaImmigrantsRatio;

int Algorithm::FitnessEvaluationSize;

//MIGA parameters
double Algorithm::RandomImmigrantsRatioForMIGA;
double Algorithm::MutationRatioOfBestOrEliteIndividual;

//MQSO parameters
unsigned int Algorithm::NumberOfSwarms;
unsigned int Algorithm::NumberOfNeutralParticlesInEachSwarm;
unsigned int Algorithm::NumberOfQuantumParticlesInEachSwarm;
double Algorithm::ControlAttractionOfGlobalBest;
double Algorithm::ControlAttractionOfPersonalBest;
double Algorithm::ConvergenceRadius;
double Algorithm::ExclusionRadius;
double Algorithm::ExclusionRadiusSquare;
double Algorithm::QuantumCloudRadius;


extern int number_of_peaks;
extern double vlength;
extern double lambda;
extern int use_basis_function;
extern double mincoordinate;
extern double maxcoordinate;

void Algorithm::InitMovPeaks() {
    number_of_peaks = Algorithm::NumberOfPeaks;
    vlength = Algorithm::ShiftLength;
    lambda = Algorithm::SOS_Lambda;
    use_basis_function = Algorithm::SOS_UseBasisFunction;
    mincoordinate = Algorithm::MinCoordinate;
    maxcoordinate = Algorithm::MaxCoordinate;

    init_peaks();
}


void Algorithm::SelectParents(Population* population, int *momIndex, int *dadIndex, size_t tournamentSize) {
    size_t size = std::min( population->getPopulationSize(), tournamentSize ); // population size may be less than tournament size

    vector<int> indices(size, -1);

    vector<Individual>& individuals = population->getIndividuals();

    unsigned int newGeneratedValue;

    // get indices
    for (int i = 0; i < size; ) {
        newGeneratedValue = (unsigned int) (rand() % population->getPopulationSize());

        vector<int>::iterator it = find(indices.begin(), indices.end(), newGeneratedValue);
        if (it != indices.end()) { // we already added this index to vector
            continue;
        }

        indices[i] = newGeneratedValue;
        i++;
    }

    // get momIndex
    double bestFitness = -DBL_MAX;
    for (int i = 0; i < size; i++) {
        double f = individuals[indices[i]].getFitness();

        if (f > bestFitness) {
            bestFitness = f;
            *momIndex = indices[i];
        }
    }

    // get dadIndex
    bestFitness = -DBL_MAX;
    for (int i = 0; i < size; i++) {
        if (indices[i] == *momIndex) { // this is mom index, so skip
            continue;
        }

        double f = individuals[indices[i]].getFitness();

        if (f > bestFitness) {
            bestFitness = f;
            *dadIndex = indices[i];
        }
    }
}



void Algorithm::ReleaseMovPeaks() {
    free_peaks();
}


void Algorithm::GenerationalGA(Population* population) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    assert( population->getPopulationSize() >= 4 );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < population->getPopulationSize() && (i+1) < population->getPopulationSize(); i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].EIGAMutate(population->getMutationStepSize());
            offspring[1].EIGAMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]) || !population->isValidIndividual(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( population->getPopulationSize() % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            //cout << ".->  " << offspring[0].getFitness() << endl;
            offspring[0].EIGAMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}


void Algorithm::RunSimpleGA() {
    srand((unsigned int)time(0));
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    int changeFreq = Algorithm::ChangePeriod;

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
        // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
    	if (total_t > changeTime) {
#endif

            statistics.addStat( population.best().getFitness() );
            statistics.addCalculateDiversity(population,Algorithm::DimensionSize);

    		statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();
            //printf(" --- Change Oldu --- Iterasyon:%d\n",i);

            change_peaks();
            population.generatePopulationRandom();

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreq= changeFreq + Algorithm::ChangePeriod;
        }

        GenerationalGA( &population );
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


void Algorithm::ApplyElitism(Population* population, const Individual &best) {
    if (Algorithm::Elitism) {
        int index;
        population->worst(&index);
        population->setIndividual( best, index );
    }
}

void Algorithm::RunRandomImmigrants() {
    srand((unsigned int)time(0));
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    int changeFreq = Algorithm::ChangePeriod;

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    size_t randomImmigrantsSize = (size_t) (Algorithm::PopulationSize * Algorithm::RandomImmigrantsChangePercentage);

    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
    	if (total_t > changeTime) {
#endif
			statistics.addStat( population.best().getFitness() );
			statistics.addCalculateDiversity(population,Algorithm::DimensionSize);

            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();
            //printf(" --- Change Oldu --- Iterasyon:%d\n",i);
            change_peaks();
            population.updateFitnesses();

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreq= changeFreq + Algorithm::ChangePeriod;
        }

        vector<Individual>& individuals = population.getIndividuals();

        vector<Individual*> worstIndividuals =
                Util::getWorstIndividuals(randomImmigrantsSize, individuals);

        for(vector<Individual*>::iterator it = worstIndividuals.begin(); it != worstIndividuals.end(); it++) {
            population.initIndividualRandom( *it );
        }

        GenerationalGA( &population );
		statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::RunMemorySearch() {
    srand((unsigned int)time(0));
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    int changeFreq = Algorithm::ChangePeriod;

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();
    int i;
    // generate population
    Population searchPopulation( 45 );
    Population memoryPopulation( 45 );

    vector<Individual> explicitMemory;

    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
    	if (total_t > changeTime) {
#endif
            statistics.addStat( memoryPopulation.best().getFitness() );
            statistics.addCalculateDiversity(memoryPopulation,Algorithm::DimensionSize);
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            // merge explicit mem and memory pop select best n
            change_peaks();

            // re-init search population
            searchPopulation.generatePopulationRandom();

            // refresh fitnesses of memory individuals
            memoryPopulation.updateFitnesses();

            statistics.addBestErrorAtChangeStat( memoryPopulation.best().getFitness() );

            // refresh fitnesses of explicit memory
            for (int j = 0; j < explicitMemory.size(); j++) {
                explicitMemory[j].updateFitness();
            }

            // update memory population
            Algorithm::UpdateMemoryPopulation( memoryPopulation, explicitMemory );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreq= changeFreq + Algorithm::ChangePeriod;
        }

        if (i > 0 && (i % Algorithm::MemoryUpdateFreq == 0)) {
            // update explicit memory
            Algorithm::UpdateExplicitMemory( memoryPopulation, searchPopulation, explicitMemory );
        }
        i++;
        // generational part of memory population
        GenerationalGA( &memoryPopulation );

        // generational part of search population
        GenerationalGA( &searchPopulation );
        statistics.addAbsoluteRecoveryRateStat( memoryPopulation.best().getFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::UpdateMemoryPopulation(Population &memoryPopulation, vector<Individual> &explicitMemory) {
    vector<Individual> mergedIndividuals = memoryPopulation.getIndividuals();
    mergedIndividuals.insert( mergedIndividuals.end(), explicitMemory.begin(), explicitMemory.end() );

    vector<Individual*> bestIndividuals = Util::getBestIndividuals( memoryPopulation.getPopulationSize(), mergedIndividuals );

    // copy best individuals to memory population
    for (int i = 0; i < memoryPopulation.getPopulationSize(); i++) {
        memoryPopulation.setIndividual( *(bestIndividuals[i]), i );
    }
}

void Algorithm::UpdateExplicitMemory(Population &memoryPopulation, Population &searchPopulation,
                                     std::vector<Individual> &explicitMemory) {
    Individual bestIndividualSearch = searchPopulation.best();
    Individual bestIndividualMemory = memoryPopulation.best();

    Individual* bestIndividual;

    if (bestIndividualMemory.getFitness() > bestIndividualSearch.getFitness())
        bestIndividual = &bestIndividualMemory;
    else
        bestIndividual = &bestIndividualSearch;


    if (explicitMemory.size() < Algorithm::ExplicitMemorySize) { // explicit memory with free slots
        explicitMemory.push_back( *bestIndividual );
    }
    else {
        // memory is full find mindist candidate to replace
        Individual* mindistIndividual = MinDistIndividual( explicitMemory, bestIndividual );

        // update explicit memory
        if (mindistIndividual->getFitness() < bestIndividual->getFitness()) {
            *mindistIndividual = *bestIndividual;
        }
    }
}

Individual *Algorithm::MinDistIndividual(std::vector<Individual> &explicitMemory, Individual *bestIndividual) {
    double mindist = DBL_MAX;
    Individual* mindistIndividual = nullptr;
    double distance;

    for (int i = 0; i < explicitMemory.size(); i++) {
        distance = bestIndividual->getEuclideanDistanceToIndividual(explicitMemory[i]);
        if ( distance < mindist ) {
            mindist = distance;
            mindistIndividual = &(explicitMemory[i]);
        }
    }

    return  mindistIndividual;
}

void Algorithm::RunSelfOrganizingScouts() {
    srand((unsigned int)time(0));
    int changeFreq = Algorithm::ChangePeriod;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    ParentPopulation parentPopulation( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif


        SOSGenerationalGA( &parentPopulation );

        SOSAdjustSearchSpace( &parentPopulation );

        SOSDeleteNonFitScouts( &parentPopulation );

        SOSAdjustPopulationSizes( &parentPopulation );

        ScoutPopulation sp = SOSFork( &parentPopulation );

        if (sp.getIndividuals().size() > 0) { // forked a child
            parentPopulation.addScoutPopulation( sp );
            parentPopulation.generatePopulationRandom();
        }

        while (SOSMergeScoutPopulations( &parentPopulation )) {
        }

        // any change?
        // any change?
#ifdef FITNESS_EVALUATION
        	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
        	if (total_t > changeTime) {
#endif
            statistics.addStat( parentPopulation.overallBestFitness() );
            statistics.addCalculateDiversitySOS(parentPopulation,Algorithm::DimensionSize);
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();
            //printf(" --- Change Oldu --- Iterasyon:%d\n",i);
            // merge explicit mem and memory pop select best n
            change_peaks();

            parentPopulation.updateAllFitnesses();

            SOSAdjustSearchSpace( &parentPopulation );

            Algorithm::statistics.addBestErrorAtChangeStat( parentPopulation.overallBestFitness() );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreq= changeFreq + Algorithm::ChangePeriod;
        }

        assert( parentPopulation.getPopulationSize() >= 8 );

        SOSReduceScoutDiameters( &parentPopulation );
//        vector<ScoutPopulation>& scouts = parentPopulation.getScoutPopulations();
//        for(int i = 0; i < scouts.size(); i++) {
//                assert( scouts[i].getPopulationSize() >= 4 &&  scouts[i].getPopulationSize() <= 10);
//        }
        statistics.addAbsoluteRecoveryRateStat( parentPopulation.overallBestFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


void Algorithm::SOSGenerationalGA(ParentPopulation *parentPopulation) {
    // generational ga steps for
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); ++it) {
        GenerationalGA( &(*it), (*it).getSuggestedSize());
    }

    GenerationalGA( parentPopulation, parentPopulation->getSuggestedSize());
}



void Algorithm::SOSAdjustPopulationSizes(ParentPopulation *parentPopulation) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);


    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();


    double obf = parentPopulation->overallBestFitness();

    // parent population (1) + number of child populations
    size_t numberOfPopulations = scouts.size() + 1;

    // get fit populations count and minimum fitness of all individuals (populations whose fitnesses are incrementing)
    size_t fitPopulationCount = 0;
    double minFitness = DBL_MAX;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {

        // get min fitness
        if ( minFitness > (*it).getFitness() ) {
            minFitness = (*it).getFitness();
        }

        // is incrementing?
        if ((*it).getFitness() > (*it).getPreviousBest().getFitness() ) {
            // will remove this population
            fitPopulationCount++;
        }
    }

    // get min fitness compare base population
    if ( minFitness > parentPopulation->getFitness() ) {
        minFitness = parentPopulation->getFitness();
    }

    // is incrementing for parent?
    // is incrementing?
    if (parentPopulation->getFitness() > parentPopulation->getPreviousBest().getFitness() ) {
        // will remove this population
        fitPopulationCount++;
    }

    double beta = (double) fitPopulationCount / numberOfPopulations;

    // calculate sumFitness
    double sumFitness = 0;
    double sumDynamism = 0;
    sumFitness += parentPopulation->getFitness() - minFitness;
    sumDynamism += parentPopulation->getDynamism();
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumFitness += (*it).getFitness() - minFitness;
        sumDynamism += (*it).getDynamism();
    }

    // arrange individuals between populations
    // first check out populations that is offered individual count > max individual count
    vector<ScoutPopulation*> minScoutSizeScouts;
    vector<ScoutPopulation*> maxScoutSizeScouts;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        double quality = Algorithm::SOSGetQuality(&(*it), beta, sumDynamism, minFitness, sumFitness);
        (*it).setQuality( quality );
    }

    double parentRelDynamism = SOSGetRelativeDynamism(parentPopulation, sumDynamism);
    double parentRelFitness = SOSGetRelativeFitness(parentPopulation, minFitness, sumFitness);
    double parentQuality = SOSGetRelativeQuality(parentPopulation, beta, sumDynamism, minFitness, sumFitness);

    double sumQuality = parentQuality;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumQuality += (*it).getQuality();
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setQuality( (*it).getQuality()/sumQuality );
    }

    parentQuality /= sumQuality;

    size_t totalScoutSizes = 0;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setSuggestedSize( SOSSuggestedSize( &(*it), (*it).getQuality() ) );
        totalScoutSizes += (*it).getSuggestedSize();
    }

    parentPopulation->setSuggestedSize( Algorithm::PopulationSize - totalScoutSizes );

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() > maxScoutSize ) {
            parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize() + (*it).getSuggestedSize()-maxScoutSize );
            (*it).setSuggestedSize( maxScoutSize );
        }
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() < minScoutSize ) {
            size_t maximumSizeScout = 0;
            ScoutPopulation* maximumScout = nullptr;
            for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
                if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                    maximumSizeScout = (*itMax).getSuggestedSize();
                    maximumScout = &(*itMax);
                }
            }

            if (maximumScout != nullptr && maximumScout->getSuggestedSize() > parentPopulation->getSuggestedSize()) {
                maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }
            else {
                parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }

            (*it).setSuggestedSize( minScoutSize );
        }
    }

    if (parentPopulation->getSuggestedSize() < minBaseSize) {
        size_t maximumSizeScout = 0;
        ScoutPopulation* maximumScout = nullptr;
        for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
            if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                maximumSizeScout = (*itMax).getSuggestedSize();
                maximumScout = &(*itMax);
            }
        }

        maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minBaseSize-parentPopulation->getSuggestedSize()) );
        parentPopulation->setSuggestedSize( minBaseSize );
    }


    size_t totalIndividualsOfScouts = parentPopulation->getSuggestedSize();
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        totalIndividualsOfScouts += (*it).getSuggestedSize();
    }

//    cout << "total: " << totalIndividualsOfScouts << endl;

    assert( totalIndividualsOfScouts == Algorithm::PopulationSize );
}

void Algorithm::SOSDeleteNonFitScouts(ParentPopulation *parentPopulation) {
    static double relativeMinFitness = 0.0;//AlgorithmInfo::instance()->getMinFitnessOfExistingForkingPopulationsRelative() * obf;

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

    vector<ScoutPopulation*> scoutsToRemove;

    // discard child populations whose fitness are less than minfitness
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        // check if child population's fitness is below minimum fitness, if so discard
        if ((*it).getFitness() < relativeMinFitness) {
            // will remove this population
            scoutsToRemove.push_back( &(*it) );
        }
    }

    // remove populations whose fitness is below min fitness
    for( vector<ScoutPopulation*>::iterator it = scoutsToRemove.begin(); it != scoutsToRemove.end(); it++) {
        parentPopulation->removeScoutPopulation( *it );
    }
}

void Algorithm::GenerationalGA(Population *population, size_t suggestedPopulationSize) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < suggestedPopulationSize && (i+1) < suggestedPopulationSize; i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());
            offspring[1].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]) || !population->isValidIndividual(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( suggestedPopulationSize % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}

void Algorithm::SOSAdjustSearchSpace(ParentPopulation *parentPopulation) {
    //FIXME parentten individuals da alınmış olunabilir.
    static double globalMinRadius = Algorithm::SOS_MinDiameterRelative
                                    * ( Algorithm::MaxCoordinate-Algorithm::MinCoordinate );

    static size_t minScoutSize = (size_t) (Algorithm::SOS_MinScoutPopulationSizeRelative * Algorithm::PopulationSize);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        vector<Individual> individualsToParent;

        // update center
        (*it).setCenter( (*it).best() );
//        // update radius
//        (*it).setDiameter( max((*it).getDiameter()*Algorithm::SOS_DiameterReduceFactor, globalMinRadius) );

        // get invalid individuals
        vector<Individual>& scoutIndividuals = (*it).getIndividuals();
        for( vector<Individual>::iterator ind = scoutIndividuals.begin(); ind != scoutIndividuals.end(); ind++ ) {
            if ( !(*it).isValidIndividual( *ind ) ) {
                individualsToParent.push_back( *ind );
            }
        }

        assert( (*it).getPopulationSize() >= 4 );
        size_t calculatedScoutSize = (*it).getPopulationSize() - individualsToParent.size();

        // add random individuals if number of individuals in scout dropped below min scout population threshold.
        if ( calculatedScoutSize < minScoutSize ) {
            size_t diff = minScoutSize - calculatedScoutSize;
            for (int i = 0; i < diff; i++) {
                (*it).initIndividualRandom( &(individualsToParent[individualsToParent.size()-i-1]) );
            }
            for (int i = 0; i < diff; i++) {
                individualsToParent.pop_back();
            }
        }

        // remove individuals from scout
        (*it).removeIndividuals( individualsToParent );

        // add these individuals to parent
        parentPopulation->addIndividuals( individualsToParent );
    }
}

double Algorithm::SOSGetRelativeDynamism(Population *population, double totalDynamism) {
    if ( totalDynamism <= 0 )
        return (double) 1/Algorithm::DimensionSize;

    return population->getDynamism()/totalDynamism;
}

double Algorithm::SOSGetRelativeFitness(Population *population, double overallMinFitness, double totalFitness) {
    if ( totalFitness <= 0 )
        return (double) 1/Algorithm::DimensionSize;

    return (population->getFitness()-overallMinFitness)/totalFitness;
}

double Algorithm::SOSGetRelativeQuality(Population *population,
                                        double beta, double totalDynamism, double overallMinFitness,
                                        double totalFitness) {
    return Algorithm::SOS_ALpha * beta * SOSGetRelativeDynamism(population, totalDynamism)
           + (1-Algorithm::SOS_ALpha*beta)* SOSGetRelativeFitness(population, overallMinFitness, totalFitness);
}

double Algorithm::SOSGetQuality(Population *population, double beta, double totalDynamism, double overallMinFitness,
                                double totalFitness) {
    return Algorithm::SOS_ALpha * beta *  population->getDynamism()
           + (1-Algorithm::SOS_ALpha*beta)* SOSGetRelativeFitness(population, overallMinFitness, totalFitness);
}

size_t Algorithm::SOSSuggestedSize(Population *population, double quality) {
    return (size_t)((quality*Algorithm::PopulationSize+population->getPopulationSize()) / 2);
}

/*
ScoutPopulation Algorithm::SOSForkDemet(ParentPopulation *parentPopulation) {
    ScoutPopulation scoutPopulation;

    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

//    // number of forking populations exceeded max. scout population size
//    while ( scouts.size() > (PopulationSize-minBaseSize)/minScoutSize) {
//        // delete worst scout
//        parentPopulation->deleteWorstScoutPopulation();
//    }

    // get most dense scout population

    // get elements whose fitness is better than min scout population fitness
    double overallBestIndividualFitness = parentPopulation->overallBestFitness();
    double minFitnessOfScout = Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative * overallBestIndividualFitness;
    if (overallBestIndividualFitness < 0)
        return scoutPopulation;


    struct FitIndividualInfo {
        int index; // index in m_individuals array
        double radius; // maximum radius calculated
        double maxRadius; // possible maximum radius (this value is used for next possible scout center)
        int maxRadiusIndex;
        vector<Individual> individuals;
    };

    // holds the index of fit individuals
    vector<FitIndividualInfo> fitIndividuals;
    vector<Individual>& individuals = parentPopulation->getIndividuals();
    for ( int i = 0; i < parentPopulation->getPopulationSize(); i++ ) {
        // cout << "curr fitness: " << m_individuals[i].getFitness() << " minfitness: " << minFitnessOfScout << endl;
        if ( individuals[i].getFitness() >= minFitnessOfScout ) {
            FitIndividualInfo fii;
            fii.index = i;
            fii.radius = 0;
            fii.maxRadius = globalMaxDiameter;
            fii.maxRadiusIndex = -1;
            fii.individuals.push_back( individuals[i] );

            fitIndividuals.push_back( fii );
        }
    }

    size_t fitSize = min( fitIndividuals.size(), maxScoutSize );

    // no suitable individual
    if (fitIndividuals.size() <= minScoutSize) {
        return scoutPopulation;
    }

    // fit individuals to other individuals distance array (distance: euclidean)
    double distanceArray[ fitIndividuals.size() ][ individuals.size() ];


    // fill distance array
    // get maximum possible radius
    for (size_t i = 0; i < fitSize; i++) {
        for (size_t j = 0; j < fitSize; j++) {
            if ( individuals[fitIndividuals[i].index].getFitness() >= individuals[fitIndividuals[j].index].getFitness() ) {
                distanceArray[i][j] = individuals[fitIndividuals[i].index].getEuclideanDistanceToIndividual(individuals[fitIndividuals[j].index]);
            }
            else {
                distanceArray[i][j] = -1;
            }

            if (distanceArray[i][j] >= 0 && distanceArray[i][j] < globalMaxDiameter) {
                fitIndividuals[i].individuals.push_back( individuals[ fitIndividuals[j].index ] );
                if (distanceArray[i][j] > fitIndividuals[i].radius) {
                    fitIndividuals[i].radius = distanceArray[i][j];
                    fitIndividuals[i].maxRadiusIndex = j;
                }
            }
        }

        if (fitIndividuals[i].maxRadiusIndex >= 0) {
            if (fitIndividuals[i].radius < globalMinDiameter) {
                fitIndividuals[i].radius = globalMinDiameter;
            }
        }
    }

    double maxDensity = 0;
    FitIndividualInfo mostDenseScout;
    for(int i = 0; i < fitSize; i++) {
        if (fitIndividuals[i].individuals.size() > (minScoutSize-1)) {
            double density = fitIndividuals[i].individuals.size() / fitIndividuals[i].radius;
            if (density > maxDensity) {
                maxDensity = density;
                mostDenseScout = fitIndividuals[i];
            }
        }
    }

    if (maxDensity <= 0)
        return scoutPopulation;


    scoutPopulation.setIndividuals( mostDenseScout.individuals );
    scoutPopulation.setSuggestedSize( mostDenseScout.individuals.size() );
    scoutPopulation.setDiameter( mostDenseScout.radius );


    vector<Individual> newIndividuals;
    // prepare individuals for parent population (scout now owns the individuals)
    size_t newIndividualsCounter = 0;

    for (int i = 0; i < individuals.size(); i++) {
        vector<Individual>::iterator found = find( mostDenseScout.individuals.begin(), mostDenseScout.individuals.end(), individuals[i] );

        if ( found != mostDenseScout.individuals.end() ) {
            continue;
        }

        newIndividuals.push_back( individuals[i] );
    }

    // update individuals
    parentPopulation->setIndividuals( newIndividuals );

    // check min parent pop size
    while ( parentPopulation->getPopulationSize() < minBaseSize
            || parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() < minBaseSize) {
        // delete worst scout
        parentPopulation->deleteWorstScoutPopulation();
    }

    parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() );

    return scoutPopulation;
}


ScoutPopulation Algorithm::SOSForkSupposedToBe(ParentPopulation *parentPopulation) {
    ScoutPopulation scoutPopulation;

    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

    // get most dense scout population

    // get elements whose fitness is better than min scout population fitness
    double overallBestIndividualFitness = parentPopulation->overallBestFitness();
    double minFitnessOfScout = Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative * overallBestIndividualFitness;
    if (overallBestIndividualFitness < 0) {
        minFitnessOfScout = (1-Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative) * overallBestIndividualFitness + overallBestIndividualFitness;
    }

    struct FitIndividualInfo {
        int index; // index in m_individuals array
        double radius; // maximum radius calculated
        double maxRadius; // possible maximum radius (this value is used for next possible scout center)
        vector<Individual> individuals;
    };

    // holds the index of fit individuals
    vector<FitIndividualInfo> fitIndividuals;
    vector<Individual>& individuals = parentPopulation->getIndividuals();
    for ( int i = 0; i < parentPopulation->getPopulationSize(); i++ ) {
        // cout << "curr fitness: " << m_individuals[i].getFitness() << " minfitness: " << minFitnessOfScout << endl;
        if ( individuals[i].getFitness() >= minFitnessOfScout ) {
            FitIndividualInfo fii;
            fii.index = i;
            fii.radius = 0;
            fii.maxRadius = globalMaxDiameter;
            fii.individuals.push_back( individuals[i] ); // add center to scout
            fitIndividuals.push_back( fii );
        }
    }

    // no suitable individual
    if (fitIndividuals.size() <= 0) {
        return scoutPopulation;
    }

    // fit individuals to other individuals distance array (distance: euclidean)
    double distanceArray[ fitIndividuals.size() ][ individuals.size() ];


    // fill distance array
    // get maximum possible radius
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < individuals.size(); j++) {
            if ( fitIndividuals[i].index == j ) {// distance to itself
                distanceArray[i][j] = 0;
            }
            else {
                distanceArray[i][j] = individuals[fitIndividuals[i].index].getEuclideanDistanceToIndividual(individuals[j]);

                // compared candidate center individual to other members, if compared member has better
                // fitness radius will not be incremented through this member
                if ( individuals[fitIndividuals[i].index].getFitness() < individuals[j].getFitness() ) {
                    if ( fitIndividuals[i].maxRadius > distanceArray[i][j] ) {
                        fitIndividuals[i].maxRadius = distanceArray[i][j];
                    }
                }
            }
        }
    }



    // get variables to calculate density
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < individuals.size(); j++) {
            if ( fitIndividuals[i].index == j ) {// distance to itself
                continue;
            }
            else {
                // FIXME burada ara elemanların eklenmesi density'i olumsuz etkileyebilir. Bu da dikkate alınmalı
                // check max scout size && radius constraints
                if ( distanceArray[i][j] < fitIndividuals[i].maxRadius
                     && fitIndividuals[i].individuals.size() < maxScoutSize) {
                    // update itemcount in scout (used in desity calculation)
                    fitIndividuals[i].individuals.push_back( individuals[j] );

                    // update radius
                    if ( fitIndividuals[i].radius < distanceArray[i][j] ) {
                        fitIndividuals[i].radius = distanceArray[i][j];
                    }
                }
            }
        }
    }


    // get density
    FitIndividualInfo mostDenseScout;
    double maxDensity = 0;
    bool foundValidScout = false;
    for (vector<FitIndividualInfo>::iterator it = fitIndividuals.begin(); it != fitIndividuals.end(); it++) {
        // check radius constraints
        if ((*it).radius <= 0 || (*it).radius < globalMinDiameter )
            continue;

        // check number of individuals constraints
        if ((*it).individuals.size() < minScoutSize)
            continue;

        double d = (double) (*it).individuals.size() / (*it).radius;
        if ( d > maxDensity ){
            foundValidScout = true;
            mostDenseScout = *it;
            maxDensity = d;
        }
    }

    // no suitable scout to fork
    if (!foundValidScout) {
        return scoutPopulation;
    }


    scoutPopulation.setIndividuals( mostDenseScout.individuals );
    scoutPopulation.setSuggestedSize( mostDenseScout.individuals.size() );
    scoutPopulation.setDiameter( mostDenseScout.radius );


    vector<Individual> newIndividuals;
    // prepare individuals for parent population (scout now owns the individuals)
    size_t newIndividualsCounter = 0;

    for (int i = 0; i < individuals.size(); i++) {
        vector<Individual>::iterator found = find( mostDenseScout.individuals.begin(), mostDenseScout.individuals.end(), individuals[i]);

        if (found != mostDenseScout.individuals.end()) {
            continue;
        }

        newIndividuals.push_back( individuals[i] );
    }

    // update individuals
    parentPopulation->setIndividuals( newIndividuals );

    // check min parent pop size
    while ( parentPopulation->getPopulationSize() < minBaseSize
            || parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() < minBaseSize) {
        // delete worst scout
        parentPopulation->deleteWorstScoutPopulation();
    }

    parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() );

    return scoutPopulation;
}

 */

ScoutPopulation Algorithm::SOSFork(ParentPopulation *parentPopulation) {
    ScoutPopulation scoutPopulation;

    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

//    // number of forking populations exceeded max. scout population size
//    while ( scouts.size() > (PopulationSize-minBaseSize)/minScoutSize) {
//        // delete worst scout
//        parentPopulation->deleteWorstScoutPopulation();
//    }

    // get most dense scout population

    // get elements whose fitness is better than min scout population fitness
    double overallBestIndividualFitness = parentPopulation->overallBestFitness();
    double minFitnessOfScout = Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative * overallBestIndividualFitness;
    if (overallBestIndividualFitness < 0) {
        minFitnessOfScout = (1-Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative) * overallBestIndividualFitness + overallBestIndividualFitness;
    }

    struct FitIndividualInfo {
        int index; // index in m_individuals array
        double radius; // maximum radius calculated
        double maxRadius; // possible maximum radius (this value is used for next possible scout center)
        vector<Individual*> individuals;
    };

    // holds the index of fit individuals
    vector<FitIndividualInfo> fitIndividuals;
    vector<Individual>& individuals = parentPopulation->getIndividuals();
    for ( int i = 0; i < parentPopulation->getPopulationSize(); i++ ) {
        // cout << "curr fitness: " << m_individuals[i].getFitness() << " minfitness: " << minFitnessOfScout << endl;
        if ( individuals[i].getFitness() >= minFitnessOfScout ) {
            FitIndividualInfo fii;
            fii.index = i;
            fii.radius = 0;
            fii.maxRadius = globalMaxDiameter;
            fii.individuals.push_back( &individuals[i] ); // add center to scout
            fitIndividuals.push_back( fii );
        }
    }

    // no suitable individual
    if (fitIndividuals.size() <= 0) {
        return scoutPopulation;
    }

    // fit individuals to other individuals distance array (distance: euclidean)
    double distanceArray[ fitIndividuals.size() ][ individuals.size() ];


    // fill distance array
    // get maximum possible radius
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < fitIndividuals.size(); j++) {
            if ( i == j ) {// distance to itself
                distanceArray[i][j] = 0;
            }
            else {
                distanceArray[i][j] = individuals[fitIndividuals[i].index].getEuclideanDistanceToIndividual(individuals[fitIndividuals[j].index]);

                // compared candidate center individual to other members, if compared member has better
                // fitness radius will not be incremented through this member
                if ( individuals[fitIndividuals[i].index].getFitness() < individuals[fitIndividuals[j].index].getFitness() ) {
                    if ( fitIndividuals[i].maxRadius > distanceArray[i][j] ) {
                        fitIndividuals[i].maxRadius = distanceArray[i][j];
                    }
                }
            }
        }
    }



    // get variables to calculate density
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < fitIndividuals.size(); j++) {
            if ( i == j ) {// distance to itself
                continue;
            }
            else {
                // FIXME burada ara elemanların eklenmesi density'i olumsuz etkileyebilir. Bu da dikkate alınmalı
                // check max scout size && radius constraints
                if ( distanceArray[i][j] < fitIndividuals[i].maxRadius
                     && fitIndividuals[i].individuals.size() < maxScoutSize) {
                    // update itemcount in scout (used in desity calculation)
                    fitIndividuals[i].individuals.push_back( &(individuals[fitIndividuals[j].index]) );

                    // update radius
                    if ( fitIndividuals[i].radius < distanceArray[i][j] ) {
                        fitIndividuals[i].radius = distanceArray[i][j];
                    }
                }
            }
        }
    }


    // get density
    FitIndividualInfo mostDenseScout;
    double maxDensity = 0;
    bool foundValidScout = false;
    for (vector<FitIndividualInfo>::iterator it = fitIndividuals.begin(); it != fitIndividuals.end(); it++) {
        // check radius constraints
        if ((*it).radius <= 0
            || (*it).radius < globalMinDiameter )
            continue;

        // check number of individuals constraints
        if ((*it).individuals.size() < minScoutSize)
            continue;

        double d = (double) (*it).individuals.size() / (*it).radius;
        if ( d > maxDensity ){
            foundValidScout = true;
            mostDenseScout = *it;
            maxDensity = d;
        }
    }

    // no suitable scout to fork
    if (!foundValidScout) {
        return scoutPopulation;
    }


    scoutPopulation.setIndividuals( mostDenseScout.individuals );
    scoutPopulation.setSuggestedSize( mostDenseScout.individuals.size() );
    scoutPopulation.setDiameter( mostDenseScout.radius );


    vector<Individual> newIndividuals;
    // prepare individuals for parent population (scout now owns the individuals)
    size_t newIndividualsCounter = 0;

    for (int i = 0; i < individuals.size(); i++) {
        bool found = false;
        for (vector<Individual*>::iterator it = mostDenseScout.individuals.begin();
             it != mostDenseScout.individuals.end(); it++) {
            if ( &(individuals[i]) == *it ) {
                found = true;
                break;
            }
        }

        if (found)
            continue;

        newIndividuals.push_back( individuals[i] );
    }

    // update individuals
    parentPopulation->setIndividuals( newIndividuals );

    // check min parent pop size
    while ( parentPopulation->getPopulationSize() < minBaseSize
            || parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() < minBaseSize) {
        // delete worst scout
        parentPopulation->deleteWorstScoutPopulation();
    }

    parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() );

    return scoutPopulation;
}

bool Algorithm::SOSMergeScoutPopulations(ParentPopulation *parentPopulation) {
    ScoutPopulation* scout1, * scout2;
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (size_t i = 0; i < scouts.size(); i++) {
        scout1 = &(scouts[i]);
        const Individual& scout1CenterIndividual = scout1->getCenter();
        for (size_t j = i+1; j < scouts.size(); j++) {
            scout2 = &(scouts[j]);
            const Individual& scout2CenterIndividual = scout2->getCenter();

            double scoutToScoutDistance =  scout1CenterIndividual.getEuclideanDistanceToIndividual( scout2CenterIndividual );

            // do we need merge?
            if ( scoutToScoutDistance < scout1->getDiameter() || scoutToScoutDistance < scout2->getDiameter()) {
                if ( scout1->getFitness() < scout2->getFitness() ) { // merge into scout2
                    MergeScoutPopulations( parentPopulation, scout2, scout1 );
                }
                else { // merge into scout1
                    MergeScoutPopulations( parentPopulation, scout1, scout2 );
                }

                return true;
            }
        }
    }

    return false;
}

void Algorithm::MergeScoutPopulations(ParentPopulation* parentPopulation, ScoutPopulation *better, ScoutPopulation *worse) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    // merge into better scout population

    double newRadius = pow(pow( better->getDiameter(), Algorithm::DimensionSize )
                           + pow( worse->getDiameter(), Algorithm::DimensionSize ), (double)1/Algorithm::DimensionSize);

    newRadius = max(min(newRadius, globalMaxDiameter), globalMinDiameter);

    double distance;
    vector<Individual>& betterScoutIndividuals = better->getIndividuals();
    vector<Individual>& worseScoutIndividuals = worse->getIndividuals();

    better->setDiameter( newRadius );
    better->setSuggestedSize( better->getSuggestedSize() + worse->getSuggestedSize() );

    vector<Individual> individualsToAdd;
    vector<Individual> individualsToParent;
    const Individual& centerInd = better->getCenter();

    // add individuals of better scout
    for (int i = 0; i < better->getPopulationSize(); i++) {
        distance = centerInd.getEuclideanDistanceToIndividual( betterScoutIndividuals[i] );
        if ( distance < better->getDiameter() ) { // add this individual to scout
            individualsToAdd.push_back( betterScoutIndividuals[i] );
        }
        else {
            individualsToParent.push_back( betterScoutIndividuals[i] );
        }
    }

    // add individuals of worse scout
    for (int i = 0; i < worse->getPopulationSize(); i++) {
        distance = centerInd.getEuclideanDistanceToIndividual( worseScoutIndividuals[i] );
        if ( distance < better->getDiameter() ) { // add this individual to scout
            individualsToAdd.push_back( worseScoutIndividuals[i] );
        }
        else {
            individualsToParent.push_back( worseScoutIndividuals[i] );
        }
    }

    better->setIndividuals( individualsToAdd );
    parentPopulation->addIndividuals( individualsToParent );

    // delete worst scout population
    for( vector<ScoutPopulation>::iterator it = parentPopulation->getScoutPopulations().begin();
         it != parentPopulation->getScoutPopulations().end(); it++) {
        if ( &(*it) == worse ) {
            parentPopulation->getScoutPopulations().erase( it );
            break;
        }
    }
}


void Algorithm::SOSReduceScoutDiameters(ParentPopulation *parentPopulation) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        // update radius
        (*it).setDiameter( max((*it).getDiameter()*Algorithm::SOS_DiameterReduceFactor, globalMinDiameter) );
    }
}

std::vector<Individual> Algorithm::SimulatedBinaryCrossover(const Individual &mom, const Individual &dad) {
    vector<Individual> offsprings;
    Individual sis = mom;
    Individual bro = dad;

    if (((double) rand() / RAND_MAX) < Algorithm::CrossoverProbability) {
        size_t lower, upper, rnd1, rnd2;

        double u, beta;

        rnd1 = rand() % Individual::getDimension();
        rnd2 = rand() % Individual::getDimension();

        if (rnd1 >= rnd2) {
            upper = rnd1;
            lower = rnd2;
        } else {
            upper = rnd2;
            lower = rnd1;
        }

        sis.setFitnessCalculated(false);
        bro.setFitnessCalculated(false);

        double* sisValues = sis.getValues();
        double* broValues = bro.getValues();
        double* momValues = mom.getValues();
        double* dadValues = dad.getValues();


        u = (double) rand() / RAND_MAX;
        if (u <= 0.5) {
            beta = pow(2 * u, (double) 1 / 3);
        } else {
            beta = (double) 1 / pow(2 * (1 - u), (double) 1 / 3);
        }
        for (size_t i = lower; i <= upper; i++) {
            sisValues[i] = 0.5 * ((1 + beta) * momValues[i] + (1 - beta) * dadValues[i]);
            broValues[i] = 0.5 * ((1 - beta) * momValues[i] + (1 + beta) * dadValues[i]);

            // fix min max coord values
            if (sisValues[i] < mincoordinate || sisValues[i] > maxcoordinate)
                sisValues[i] = (momValues[i] + dadValues[i]) / 2;

            if (broValues[i] < mincoordinate || broValues[i] > maxcoordinate)
                broValues[i] = (momValues[i] + dadValues[i]) / 2;
        }
    }

    offsprings.push_back( sis );
    offsprings.push_back( bro );

    return offsprings;
}



void Algorithm::RunHyperMutation() {
    srand((unsigned int)time(0));
    int changeFreq = Algorithm::ChangePeriod;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    clock_t start = clock();
    start_t = clock();
    // run generations
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
    	if (total_t > changeTime) {
#endif
            statistics.addStat( population.best().getFitness() );
            statistics.addCalculateDiversity(population,Algorithm::DimensionSize);
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            change_peaks();
            population.updateFitnesses();
            Individual::setHyperMutationEnabled( true );

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreq= changeFreq + Algorithm::ChangePeriod;
        }

        GenerationalGA( &population );

        Individual::setHyperMutationEnabled( false );
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


double Algorithm::findBestFitnessforCPSOR(vector<PopPSO> &plst)
{
	double bestFitness = -100;

	for (int k = 0; k < (int) plst.size(); ++k) {
		if (plst[k].gBest.getFitness() > bestFitness) {
			bestFitness = plst[k].gBest.getFitness() ;
		}
	}
	return bestFitness;
}

void Algorithm::RunEIGA() {
    srand((unsigned int)time(0));
    int changeFreg = Algorithm::ChangePeriod;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );
    int i = 0;
    Algorithm::Elitism = false;
    size_t eigaImmigrantsSize = (size_t) (Algorithm::PopulationSize * Algorithm::EigaImmigrantsRatio);

    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreg) {
#else
    	if (total_t > changeTime) {
#endif
            statistics.addCalculateDiversity(population,Algorithm::DimensionSize);
            statistics.addStat( population.best().getFitness() );
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();
            //printf(" --- Change Oldu --- Iterasyon:%d\n",i);
            change_peaks();
            population.generatePopulationRandom();

            statistics.addBestErrorAtChangeStat(population.best().getFitness() );
            changeFreg = changeFreg + Algorithm::ChangePeriod;
            changeTime= changeTime + Algorithm::ChangeTime;
        }

        Individual bestBefore = population.best();
        Individual worstNow = population.worst();

        if(i > 0){
        	bestBefore = population.getPreviousBest();

        	int index;
			population.worst(&index);
			population.setIndividual( bestBefore, index );
        }

        vector<Individual>& individuals = population.getIndividuals();

		vector<Individual*> worstIndividuals =
				Util::getWorstIndividuals(eigaImmigrantsSize, individuals);

		for(vector<Individual*>::iterator it = worstIndividuals.begin(); it != worstIndividuals.end(); it++) {
			Individual* individual = *it;
			Individual bestindividual = bestBefore;


			bestindividual.EIGAMutate(population.getMutationStepSize());
			bestindividual.updateFitness();

			*individual = bestindividual;
			individual->setPopulation( &population );

		}

        GenerationalGA( &population );

        population.setPreviousBest(population.best());


        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
        i++;
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
    Algorithm::Elitism = true;
}


void Algorithm::RunCPSOR() {
    srand((unsigned int)time(0));
    int changeFreg = Algorithm::ChangePeriod;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    Algorithm::InitMovPeaks();

    // generate population
    PopPSO population( Algorithm::PopulationSize );

    //MNA: Global sub_population
    vector<PopPSO> plst;

    //MNA:Particle Best ilk initial deger atamasi
    for(int i=0;i<(int)Algorithm::PopulationSize;i++){
    	IndPSO *tTemp =  population.getIndividualA(i);
    	for(int j=0;j<(int)Algorithm::DimensionSize;j++ ){
    		tTemp->particleBestLocation[j] = tTemp->location[j];
    		//printf(" %f, ",tTemp->particleBestLocation[j]);
    	}
    	//printf("\n");
    }

    clock_t start = clock();

    //clustering
    ClusteringCPSOR(&population, &plst);

    Algorithm::statistics.clear();
    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreg) {
#else
    	if (total_t > changeTime) {
#endif
            statistics.addStat( findBestFitnessforCPSOR(plst));
            statistics.addDiversityCPSOR(statistics.calcDiversityCPSOR(plst, Algorithm::DimensionSize));
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();
            //printf(" --- Change Oldu --- Iterasyon:%d\n",i);
            change_peaks();
			for (int j = 0; j < (int) plst.size(); j++) {
				plst[j].updateFitnesses();
				plst[j].gBest.updateFitness();
			}
            statistics.addBestErrorAtChangeStat(findBestFitnessforCPSOR(plst));
            changeFreg= changeFreg + Algorithm::ChangePeriod;
            changeTime= changeTime + Algorithm::ChangeTime;
        }

    	/**/
    	for(vector<PopPSO>::iterator it = plst.begin(); it != plst.end();it++){
    		PopPSO *cTemp = &(*it);
			//MNA: Burada herbir Cluster PSO ya gonderilecek
			cTemp->PSO();
		}

    	removePLST(plst);

    	int survived_ind=0;
    	for(vector<PopPSO>::iterator it = plst.begin(); it != plst.end();it++){
    		PopPSO cTemp = (*it);
    		survived_ind = survived_ind + cTemp.getPopulationSize();
    	}

    	if(survived_ind < Algorithm::PopulationSize * A_DiversityThreshold){
    		//cout<< "*** Random individuals Inserted randomly"<<endl;
    		PopPSO temp_population(Algorithm::PopulationSize - survived_ind );
    		ClusteringCPSOR(&temp_population, &plst);
    	}


        //printf("%d-BestFitness: %f PopSize:%d\t Sub:%d \n",i, gBest.getFitness(), survived_ind, plst.size());
        statistics.addAbsoluteRecoveryRateStat( findBestFitnessforCPSOR(plst));
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );


    Algorithm::ReleaseMovPeaks();
}

void Algorithm::ClusteringCPSOR(PopPSO* population, vector<PopPSO> *plst){
	vector<PopPSO> sub_population;

	for (int i = 0; i < population->getPopulationSize(); i++) {
		PopPSO cTemp;
		cTemp.addIndividual(population->getIndividual(i));
		sub_population.push_back(cTemp);
	}

	double **mMatrix = NULL;

	int mMatrix_Size = sub_population.size();
	mMatrix = calculatemMatrix(sub_population);

	int tIndex,sIndex;
	int k=0;
	while (true){
		if(!findNearestPairs(sub_population,&tIndex,&sIndex,mMatrix))
			break;

		//cout<<" Buldum t:"<< tIndex << " s:" << sIndex <<endl;

		//Merge islemi yapilmasi gerekiyor
		mergeClusters(&sub_population[tIndex],&sub_population[sIndex]);
		sub_population.erase(sub_population.begin()+sIndex);


		deletemMatrix(mMatrix,mMatrix_Size);
		mMatrix_Size = sub_population.size();
		mMatrix = calculatemMatrix(sub_population);


		//MNA: each cluster in sub_cluster has more than one individual
		bool moreThanOne = true;
		for (int i = 0; i < (int)sub_population.size() && moreThanOne; i++) {
			if(sub_population[i].getPopulationSize() <= 1){
				moreThanOne = false;
			}
		}
		if(moreThanOne)
			break;

	}
	/*
	for (int i = 0; i < sub_population.size(); i++) {
		cout<< " i:" << i<<endl;
		sub_population[i].display();
	}*/

	//MNA: sub_poputaion global sub_populationa aktarilmasi gerekiyor
	//delete population; HATA veriyor
	deletemMatrix(mMatrix,mMatrix_Size);

	for(vector<PopPSO>::iterator it = sub_population.begin(); it != sub_population.end();it++){
		PopPSO cTemp = (*it);
		cTemp.calculateRadius();
		//cout<< " Radius:"<<cTemp.getRadius()<<endl;
		plst->push_back(cTemp);
	}

}

bool Algorithm::findNearestPairs(std::vector<PopPSO>& clusters, int *t, int *s, double **mMatrix){
	bool found=false;
	double min_dist = DBL_MAX;
	int clusters_size = clusters.size();

	for(int i=0;i< clusters_size;i++){
		for(int j=i+1;j< clusters_size;j++){
			if((clusters[i].getPopulationSize()+clusters[j].getPopulationSize()) > subSize)
				continue;
			if(min_dist > mMatrix[i][j]){
				min_dist = mMatrix[i][j];
				*t=i;
				*s=j;
				found = true;
			}
		}
	}
	return found;
}

void Algorithm::mergeClusters(PopPSO *t, PopPSO *s){
	for(int i=0; i < s->getPopulationSize();i++ ){
		t->addIndividual(s->getIndividual(i));
	}
}

double** Algorithm::calculatemMatrix(std::vector<PopPSO>& clusters){
	int i,j;
	int numberofclusters = (int)clusters.size();
	double **Matrix = new double *[numberofclusters];

	vector<IndPSO> cluster_s;
	vector<IndPSO> cluster_t;

	for(i=0; i< numberofclusters; i++) {

		Matrix[i] = new double[numberofclusters];
		cluster_s = clusters[i].getIndividuals();

		for(j=i+1;j < numberofclusters; j++){

			cluster_t = clusters[j].getIndividuals();

			double temp_dist;
			double min_dist = DBL_MAX;
			for(int m=0;m< (int)cluster_s.size();m++){
				for(int n=0;n< (int)cluster_t.size();n++){
					temp_dist = calculateDistance(cluster_s[m], cluster_t[n]);
					if(temp_dist<min_dist)min_dist = temp_dist;
				}
			}
			Matrix[i][j] = min_dist;
			//cout << " ->" << " i:" << i << " j:" << j << " M:"<< Matrix[i][j] <<endl;
		}
		//cout << i << endl;
	}
	return Matrix;
}

void Algorithm::deletemMatrix(double **matrix,int size){
	for(int i=0;i<size;i++){
		delete [] matrix[i];
	}
	delete [] matrix;
	matrix = NULL;
}

double Algorithm::calculateDistance(IndPSO& i, IndPSO& j){
	//cout << " X " << i.getEuclideanDistanceToIndividual(j) <<endl;
	return i.getEuclideanDistanceToIndividualPSO(j);
}

void Algorithm::removePLST(std::vector<PopPSO>& clusters){

	/*
	for(int i=0;i< (int)clusters.size();i++){
		cout<<"Cluster "<<i<<":"<<endl;
		clusters[i].display();
	}*/

	//MNA: Merge and remove pair sets
	for(int i=0;i<(int)clusters.size();i++){
		for(int j=i+1;j<(int)clusters.size();j++){
			if( ratioOverlap(&clusters[i], &clusters[j]) > B_OverlabThreshold){
				mergeClusters(&clusters[i], &clusters[j]);
				clusters.erase(clusters.begin()+j);
				//cout<< " -- Merge And Remove Yapildi "<<endl;
			}
		}
	}

	//MNA: Control the cluster size with subSize
	for(int i=0;i< (int)clusters.size();i++){
		if(clusters[i].getPopulationSize() > subSize){
			int numberOfRemoveInd = clusters[i].getPopulationSize() - subSize;
			//cout<< " -C"<< i <<" Remove edilecek birey sayisi:"<< numberOfRemoveInd <<endl;
			vector<IndPSO>& tempInd = clusters[i].getIndividuals();
			vector<IndPSO*> removeWorsInd = Util::getWorstIndividualsPSO(numberOfRemoveInd,tempInd);
			clusters[i].removeIndividuals(removeWorsInd);
		}
	}


	//MNA: Chech cluster radius and delete
	for(int i=0;i< clusters.size();i++){

		clusters[i].calculateRadius();
		//printf("  -%d Radius:%f\n",i,clusters[i].getRadius());
		//clusters[i].display();
		if(clusters[i].getRadius() < E_ClusterRadiusThreshold){
			//remove Cluster members
			clusters.erase(clusters.begin() + i);
			//cout<< " -- Cluster Remove Edildi :"<< i <<endl;

		}
	}

}

double Algorithm::ratioOverlap(PopPSO *t, PopPSO *s){
	double Roverlap = 0.0;
	double clusters_distance;
	double total_cluster_radius;


	clusters_distance = Util::GetEuclideanDistanceByValues(t->getScenter(), s->getScenter());
	total_cluster_radius = t->getRadius() + s->getRadius();

	//MNA: Not overlaped
	if(total_cluster_radius < clusters_distance){
		//cout << "NOT OVERLAP"<<endl;
		return 0;
	}


	double t_ratio=0.0, s_ratio=0.0;
	for(vector<IndPSO>::iterator it = t->getIndividuals().begin(); it != t->getIndividuals().end();it++){
		double ind_cluster_distance;
		ind_cluster_distance = Util::GetEuclideanDistanceByValues(it->getLocation(), s->getScenter());
		//cout<<"  T>ind_cluster_distance:"<<ind_cluster_distance<< " Radius:"<<s->getRadius()<<endl;
		if (ind_cluster_distance < s->getRadius())
			t_ratio ++;
	}
	t_ratio = t_ratio / t->getPopulationSize();


	for(vector<IndPSO>::iterator it = s->getIndividuals().begin(); it != s->getIndividuals().end();it++){
		double ind_cluster_distance;
		ind_cluster_distance = Util::GetEuclideanDistanceByValues(it->getLocation(), t->getScenter());
		//cout<<"  S>ind_cluster_distance:"<<ind_cluster_distance<< " Radius:"<<t->getRadius()<<endl;
		if (ind_cluster_distance < t->getRadius())
			s_ratio ++;
	}
	s_ratio = s_ratio / s->getPopulationSize();


	if(t_ratio > s_ratio)
		Roverlap = s_ratio;
	else Roverlap = t_ratio;

	//cout <<" t_ratio:"<<t_ratio<<" s_ratio:"<<s_ratio<<endl;
	//cout<< " total_cluster_radius:"<<total_cluster_radius<<" clusters_distance:"<<clusters_distance<<" Roverlap:"<< Roverlap<<endl;

	return Roverlap;
}


void Algorithm::RunMemoryBasedImmigrantsGeneticAlgorithm() {
    srand((unsigned int)time(0));
    int changeFreg = Algorithm::ChangePeriod;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );
    Population memoryPopulation( Algorithm::ExplicitMemorySize );

    int timerForMemoryUpdate= 5 + rand()%6 ;
    size_t randomImmigrantsSize = (size_t) (Algorithm::PopulationSize * Algorithm::RandomImmigrantsRatioForMIGA);

    clock_t start = clock();

    int i;
    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif

    	int indexOfWorstIndividual;
    	population.worst(&indexOfWorstIndividual);
    	population.setIndividual(population.getPreviousBest(),indexOfWorstIndividual);
#ifdef FITNESS_EVALUATION
    	 if (i > 0 && ((Algorithm::FitnessEvaluationSize > changeFreg) || (timerForMemoryUpdate == i))) {

    	    if((Algorithm::FitnessEvaluationSize > changeFreg)){
#else
        if (i > 0 && ((total_t > changeTime) || (timerForMemoryUpdate == i))) {

        	if((total_t > changeTime)){
#endif
                statistics.addStat( population.best().getFitness() );
                statistics.addCalculateDiversity(population,Algorithm::DimensionSize);

		        //printeverythingforMIGA (population,memoryPopulation);
				statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
				statistics.clearAbsoluteRecoveryRateStat();
				change_peaks();
	            population.updateFitnesses();
	            memoryPopulation.updateFitnesses();
		        //printeverythingforMIGA (population,memoryPopulation);
            	statistics.addBestErrorAtChangeStat( population.best().getFitness() );
	            changeTime= changeTime + Algorithm::ChangeTime;
	            changeFreg= changeFreg + Algorithm::ChangePeriod;
        	}
        	Individual bestOfPopulation;
        	if(timerForMemoryUpdate == i)
        		bestOfPopulation = population.best();
        	if((i % Algorithm::ChangePeriod == 0))
        		bestOfPopulation = population.getPreviousBest();
        	if(memoryPopulation.hasRandomIndividual()){
        		memoryPopulation.replaceRandomIndividualOfMemoryWithBestOrElite(bestOfPopulation);
        	}else{
            	if(timerForMemoryUpdate == i){
            		Individual* nearestBestIndividual = memoryPopulation.getClosestIndivual(bestOfPopulation);
            		if(bestOfPopulation.getFitness()>nearestBestIndividual->getFitness()){
                    	memcpy(nearestBestIndividual->getValues(), bestOfPopulation.getValues(), Algorithm::DimensionSize*sizeof(double));
                    	nearestBestIndividual->setFitnessCalculated(false);
            		}
            	}
            	if((i % Algorithm::ChangePeriod == 0)){
            		Individual* nearestBestIndividual = memoryPopulation.getClosestIndivual(bestOfPopulation);
            		if(bestOfPopulation.getFitness()>nearestBestIndividual->getFitness()){
                    	memcpy(nearestBestIndividual->getValues(), bestOfPopulation.getValues(), Algorithm::DimensionSize*sizeof(double));
                    	nearestBestIndividual->setFitnessCalculated(false);
            		}
            	}
        	}

        	timerForMemoryUpdate= i + 5 + rand()%6 ;
        }

        vector<Individual>& individuals = population.getIndividuals();

        vector<Individual*> worstIndividuals = Util::getWorstIndividuals(randomImmigrantsSize, individuals);

        Individual bestIndividualInMemory = memoryPopulation.best();


        for(vector<Individual*>::iterator it = worstIndividuals.begin(); it != worstIndividuals.end(); ++it) {
            (*it)->mutateAccordingToIndividual( bestIndividualInMemory );
        }


        GenerationalGA( &population );
        i++;
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::RunMultiQuantumSwarmAlgorithm() {
    srand((unsigned int)time(0));
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    // Population population( Algorithm::PopulationSize );
    // generate Swarms
   int changeFreg = Algorithm::ChangePeriod;

    SwarmManager swarmManager( Algorithm::NumberOfSwarms, Algorithm::NumberOfNeutralParticlesInEachSwarm, Algorithm::NumberOfQuantumParticlesInEachSwarm);


    clock_t start = clock();

    // run generations
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
            // any change?
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreg) {
#else
    	if (total_t > changeTime) {
#endif
            statistics.addStat( swarmManager.getBestSwarmFitness()   );
            double max_diversity = 0;
    		double d_temp = 0;
    		for (int m = 0; m < (int) swarmManager.getSwarmManagerSize(); m++) {
    			d_temp = statistics.addCalculateDiversityInmQSO(swarmManager.getSwarm(m), Algorithm::DimensionSize);
    			if (m == 0)
    				max_diversity = d_temp;
    			else if (max_diversity < d_temp)
    				max_diversity = d_temp;
    		}
           	statistics.addDiversitymQSO(max_diversity);
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            change_peaks();
            //Add environment change code from here

            swarmManager.setReInitsForEnvironmentChange();
            swarmManager.triggerFitnessCalculation();
            swarmManager.updateFitnesses();

            //population.generatePopulationRandom();
            //Up to here

            statistics.addBestErrorAtChangeStat( swarmManager.getBestSwarmFitness() );
            changeTime= changeTime + Algorithm::ChangeTime;
            changeFreg= changeFreg + Algorithm::ChangePeriod;
        }else{
            swarmManager.updateFitnesses();
        }

        //Add environment change code from here

        swarmManager.testForConvergance();
        swarmManager.testForExlusion();

        swarmManager.move();
#ifdef MQSO_DEBUG_EXIT
        printf("mqSQ Evaulation size = %d - gbest fitness= %f\n\r", Algorithm::FitnessEvaluationSize, swarmManager.getBestSwarmFitness());
#endif
        statistics.addAbsoluteRecoveryRateStat( swarmManager.getBestSwarmFitness()   );
        end_t = clock();
    	total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();

}
