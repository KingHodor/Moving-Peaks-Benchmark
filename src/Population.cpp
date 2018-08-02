

#include <float.h>
#include <assert.h>
#include "Population.h"
#include "Algorithm.h"
#include "Util.h"

#include <vector>
#include <algorithm>
#include <assert.h>
#include <cstring>

Population::Population()
        : mutationStepSize(3.3) {
}

Population::Population(size_t size)
        : mutationStepSize(3.3) {
    // allocate individuals array

    for (int i = 0; i < size; i++) {
        Individual individual;
        initIndividualRandom( &individual );
        individual.setPopulation( this );
        individuals.push_back( individual );
    }

    clock_gettime( CLOCK_MONOTONIC, &populationCreationTime );
}

Population::~Population() {

}

void Population::generatePopulationRandom() {
    for (vector<Individual>::iterator i = individuals.begin(); i != individuals.end(); i++) {
        initIndividualRandom( &(*i) );
    }
}

vector<Individual>& Population::getIndividuals() {
    return individuals;
}

Individual& Population::getIndividual(size_t index) {
    assert( index < individuals.size() && index >= 0 );
    return individuals[index];
}

void Population::setIndividuals(const vector<Individual> &inds) {
    individuals = inds;
}

void Population::setIndividuals(const vector<Individual *> &inds) {
    individuals.clear();
    for(vector<Individual*>::const_iterator it = inds.begin(); it != inds.end(); it++) {
        Individual ind;
        ind = *(*it);
        ind.setPopulation( this );
        individuals.push_back( ind );
    }
}

void Population::setIndividual(const Individual &ind, size_t index) {
    assert( index < individuals.size() && index >= 0 );
    individuals[index] = ind;
}

Individual Population::best(int *index) {
    int bestIndex = 0;
    double bestFitness = -DBL_MAX;
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].getFitness() > bestFitness) {
            bestFitness = individuals[i].getFitness();
            bestIndex = i;
        }
    }

    if (index)
        *index = bestIndex;

    return individuals[bestIndex];
}

Individual Population::worst(int *index) {
    int worstIndex = 0;
    double worstFitness = DBL_MAX;
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].getFitness() < worstFitness) {
            worstFitness = individuals[i].getFitness();
            worstIndex = i;
        }
    }

    if (index)
        *index = worstIndex;

    return individuals[worstIndex];
}

size_t Population::getPopulationSize() const {
    return individuals.size();
}

void Population::display() const {
    for (vector<Individual>::const_iterator it = individuals.begin(); it != individuals.end(); it++) {
        (*it).display();
    }
}

double Population::getMutationStepSize() const {
    return mutationStepSize;
}

void Population::setMutationStepSize(double size) {
    mutationStepSize = size;
}

bool Population::isValidIndividual( const Individual& individual ) const {
    return true;
}

bool Population::isValidIndividualMovingCenter(const Individual &individual) {
    return true;
}

void Population::updateFitnesses() {
    for (vector<Individual>::iterator i = individuals.begin(); i != individuals.end(); ++i) {
        (*i).updateFitness();
    }
}

void Population::setPreviousBest(const Individual &prevBest) {
    previousBest = prevBest;
}

double Population::getFitness() {
    return best().getFitness();
}

double Population::getDynamism() {
    return std::max(0.0, (getFitness()-previousBest.getFitness())/previousBest.getFitness());
}

void Population::initIndividualRandom(Individual *individual) {
    double* values = individual->getValues();
    for (int i = 0; i < Algorithm::DimensionSize; i++) {
        values[i] = Util::RandomDouble(Algorithm::MinCoordinate, Algorithm::MaxCoordinate);
    }

    individual->updateFitness();
}

void Population::replaceIndividualElitist(Individual *individual, Individual *mutateElite)
{
	double* values = individual->getValues();
	double* eliteValues = mutateElite->getValues();
	for (int i = 0; i < Algorithm::DimensionSize; i++) {
		values[i] = eliteValues[i];
	}

	individual->updateFitness();
}

void Population::initIndividualRandomScoutAware(Individual *individual) {
    initIndividualRandom( individual );
}

void Population::removeIndividuals(const vector<Individual> &individualsToRemove) {
    for( vector<Individual>::const_iterator it = individualsToRemove.begin(); it != individualsToRemove.end(); it++ ) {
        vector<Individual>::iterator found = find( individuals.begin(), individuals.end(), *it );
        assert( found != individuals.end() );
        individuals.erase( found );
    }
}

void Population::addIndividuals(const vector<Individual *> &individualsToAdd) {
    for( vector<Individual*>::const_iterator it = individualsToAdd.begin(); it != individualsToAdd.end(); it++ ) {
        Individual ind = *(*it);
        ind.setPopulation( this );
        individuals.push_back( ind );
    }
}

void Population::addIndividuals(const vector<Individual> &individualsToAdd) {
    for( vector<Individual>::const_iterator it = individualsToAdd.begin(); it != individualsToAdd.end(); it++ ) {
        Individual ind = *it;
        ind.setPopulation( this );
        individuals.push_back( ind );
    }
}

Individual &Population::getPreviousBest() {
    return previousBest;
}


void Population::addIndividual(const Individual &ind) {
    Individual i = ind;
    i.setPopulation( this );
    individuals.push_back( i );
}

void Population::removeIndividuals(const vector<Individual *> &individualsToRemove) {
    for( vector<Individual*>::const_iterator it = individualsToRemove.begin(); it != individualsToRemove.end(); it++ ) {
//        for(int i = 0; i < individuals.size(); i++) {
//            if ( individuals[i] == *(*it) ) {
//                std::cout << "FOUND" << std::endl;
//                break;
//            }
//        }
        vector<Individual>::iterator found = find( individuals.begin(), individuals.end(), *(*it) );
//        if (found == individuals.end()) {
//            (*it)->display();
//            std::cout << "-------------" << std::endl;
//
//            display();
//        }

        //assert( found != individuals.end() );
        individuals.erase( found );
    }
}


vector<double> Population::getTotalLinRank() {
    double sum = 0;
    vector<double> linRanks;
    for(int i=0;i<(int)individuals.size();i++){
        sum+=individuals[i].getLinRank();
        linRanks.push_back(sum);
    }

    return linRanks;
}


bool Population::hasRandomIndividual() {
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].isRandomlyCreated()) {
        	return true;
        }
    }
    return false;
}

void Population::replaceRandomIndividualOfMemoryWithBestOrElite(const Individual &ind) {
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].isRandomlyCreated()) {
        	std::memcpy(individuals[i].getValues(), ind.getValues(), Algorithm::DimensionSize*sizeof(double));
        	individuals[i].setRandomlyCreated(false);
        	individuals[i].setFitnessCalculated(false);
        	return;
        }
    }
    return;
}

Individual* Population::getClosestIndivual(const Individual &ind) {
	double closestDistance = DBL_MAX;
	int closestDistanceIndex;

    for (int i = 0; i < individuals.size(); i++) {
    	double distance = 0;
    	for(int dimension=0; dimension < Algorithm::DimensionSize; dimension++){
    		double currentDiff = abs(individuals[i].getValues()[dimension] - ind.getValues()[dimension]);
    		distance += currentDiff * currentDiff;
        }
    	if(distance<closestDistance){
    		closestDistanceIndex = i;
    		closestDistance = distance;
    	}
    }
    return &individuals[closestDistanceIndex];
}
