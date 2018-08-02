/*
 * Individual.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <assert.h>
#include "Population.h"
#include "Algorithm.h"
#include <stdlib.h>
#include "Util.h"
#include "movpeaks.h"

size_t Individual::Dimension = 5;

bool Individual::HyperMutationEnabled = false;


Individual::Individual()
        : fitnessCalculated(false),
          fitness(-100), population(nullptr), linRank(0), rank(-1) {
    values = new double[Individual::Dimension];
}

Individual::Individual(const Individual &ind) {
    values = new double[Individual::Dimension];
    fitness = ind.fitness;
    fitnessCalculated = ind.fitnessCalculated;
    population = ind.population;
    mutation = ind.mutation;
    rank = ind.rank;
    linRank = ind.linRank;

    memcpy(values, ind.values, Individual::Dimension*sizeof(double));
}

Individual::~Individual() {
    if (values)
        delete [] values;
}


Individual &Individual::operator=(const Individual &ind) {
    if (this == &ind)
        return *this;

    if (values)
        delete [] values;

    values = new double[Individual::Dimension];
    fitness = ind.fitness;
    fitnessCalculated = ind.fitnessCalculated;
    population = ind.population;
    mutation = ind.mutation;
    rank = ind.rank;
    linRank = ind.linRank;

    memcpy(values, ind.values, Individual::Dimension*sizeof(double));

    return *this;
}

bool Individual::operator==(const Individual &ind) const {
    assert( values != nullptr && ind.values != nullptr );

    if (fitness != ind.fitness
        || fitnessCalculated != ind.fitnessCalculated
        || population != ind.population) {
        return false;
    }

    if ( values != nullptr && ind.values != nullptr ) {
        for (int i = 0; i < Individual::Dimension; ++i) {
            if ( values[i] != ind.values[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }

    return true;
}

void Individual::setDimension(size_t dimension) {
    Individual::Dimension = dimension;
}

size_t Individual::getDimension() {
    return Individual::Dimension;
}

double *Individual::getValues() const {
    return values;
}

void Individual::display() const {
    std::cout << "Individual [";
    for (int i = 0; i < Individual::Dimension-1; i++)
        std::cout << values[i] << ",";
    std::cout << values[Individual::Dimension-1] << "] fitness: " << fitness << std::endl;
}

double Individual::getFitness() {
    if (!fitnessCalculated)
        updateFitness();

    return fitness;
}

void Individual::updateFitness() {
    fitness = dummy_eval_fitness( values );
    fitnessCalculated = true;
    Algorithm::NumberOfFitnessEvaluation ++;
}

void Individual::gaussianMutate(double std) {
    double temp;

    double mutationProb = Algorithm::MutationProbability;
    if (Individual::HyperMutationEnabled)
        mutationProb = Algorithm::HyperMutationProbability;

    for (int i = 0; i < Individual::Dimension; i++) {
        if ( ((double)rand() / RAND_MAX) < mutationProb ) {
            fitnessCalculated = false;

            temp = Util::gaussRND(values[i], std);

            if (temp > Algorithm::MaxCoordinate)
                values[i] = Algorithm::MaxCoordinate;
            else if (temp < Algorithm::MinCoordinate)
                values[i] = Algorithm::MinCoordinate;
            else
                values[i] = temp;
        }
    }
}

void Individual::EIGAMutate(double std) {
	double temp;

	double mutationProb = Algorithm::MutationProbability;

	for (int i = 0; i < Individual::Dimension; i++) {
		if ( ((double)rand() / RAND_MAX) < mutationProb ) {

			fitnessCalculated = false;

			temp = Util::gaussRND(values[i], std);

			if (temp > Algorithm::MaxCoordinate)
				values[i] = Algorithm::MaxCoordinate;
			else if (temp < Algorithm::MinCoordinate)
				values[i] = Algorithm::MinCoordinate;
			else
				values[i] = temp;
		}
	}
}


bool Individual::isFitnessCalculated() const {
    return fitnessCalculated;
}

void Individual::setFitnessCalculated(bool fitnessCalculated) {
    Individual::fitnessCalculated = fitnessCalculated;
}

double Individual::getEuclideanDistanceToIndividual(const Individual &individual) const {
    return Util::GetEuclideanDistance( *this, individual );
}

void Individual::setHyperMutationEnabled(bool enabled) {
    Individual::HyperMutationEnabled = enabled;
}

bool Individual::isHyperMutationEnabled() {
    return Individual::HyperMutationEnabled;
}

bool Individual::isRandomlyCreated(){
	return this->randomness;
}

void Individual::setRandomlyCreated(bool randomlyCreated) {
    this->randomness = randomlyCreated;
}

void Individual::mutateAccordingToIndividual(const Individual &reference) {

    double mutationProb = Algorithm::MutationRatioOfBestOrEliteIndividual;

    for (int i = 0; i < Individual::Dimension; i++) {
        if ( ((double)rand() / RAND_MAX) < mutationProb ) {
            this->fitnessCalculated = false;

            double temp = Util::gaussRND(reference.values[i], 3.3);
            //cout << "dimension:" << i << " oldVal: " << this->values[i] << "refVal: " << reference.values[i] <<  " newVal: " << temp << "\n";

            if (temp > Algorithm::MaxCoordinate)
                this->values[i] = Algorithm::MaxCoordinate;
            else if (temp < Algorithm::MinCoordinate)
            	this->values[i] = Algorithm::MinCoordinate;
            else
            	this->values[i] = temp;
        }
    }
}
