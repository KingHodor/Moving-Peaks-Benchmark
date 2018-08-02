
#include "ScoutPopulation.h"
#include "Util.h"
#include "Algorithm.h"

#include <iostream>
#include <cmath>


using namespace std;

ScoutPopulation::ScoutPopulation()
{

}

ScoutPopulation::ScoutPopulation(size_t size) : Population(size), suggestedSize(size) {

}

bool ScoutPopulation::isValidIndividual(const Individual &individual) const {
    return Util::GetEuclideanDistance( center, individual ) <= diameter;
}


void ScoutPopulation::initIndividualRandom(Individual *individual) {
    static size_t numDimension = Algorithm::DimensionSize;

    double* coords = center.getValues();
    double* values = individual->getValues();

    for (int i = 0; i < numDimension; i++) {
        values[i] = Util::RandomDouble(coords[i]-diameter/sqrt(numDimension), coords[i]+diameter/sqrt(numDimension));
    }

    individual->updateFitness();
}


void ScoutPopulation::setIndividuals(const vector<Individual *> &inds) {
    Population::setIndividuals(inds);
    center = best();
}

void ScoutPopulation::setIndividuals(const vector<Individual> &inds) {
    Population::setIndividuals(inds);
    center = best();
}

bool ScoutPopulation::isValidIndividualMovingCenter(const Individual &individual) {
    return Util::GetEuclideanDistance( best(), individual ) <= diameter;
}
