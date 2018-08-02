/*
 * IndPSO.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <assert.h>
#include "Algorithm.h"
#include <stdlib.h>

#include "IndPSO.h"
#include "Util.h"
#include "movpeaks.h"

size_t IndPSO::Dimension = 5;

IndPSO::IndPSO()
        : fitnessCalculated(false),
          fitness(-100), populationPSO(nullptr), linRank(0), rank(-1) {
    location = new double[IndPSO::Dimension];
	velocity = new double[IndPSO::Dimension];
	particleBestLocation = new double[IndPSO::Dimension];
}

IndPSO::IndPSO(const IndPSO &ind) {
    location = new double[IndPSO::Dimension];
	velocity = new double[IndPSO::Dimension];
	particleBestLocation = new double[IndPSO::Dimension];
	
    fitness = ind.fitness;
    fitnessCalculated = ind.fitnessCalculated;
    populationPSO = ind.populationPSO;
    rank = ind.rank;
    linRank = ind.linRank;

    memcpy(location, ind.location, IndPSO::Dimension*sizeof(double));
	memcpy(velocity, ind.velocity, IndPSO::Dimension*sizeof(double));
	memcpy(particleBestLocation, ind.particleBestLocation, IndPSO::Dimension*sizeof(double));
}

IndPSO::~IndPSO() {
    if (location)
        delete [] location;
	if (velocity)
        delete [] velocity;
	if (particleBestLocation)
	        delete [] particleBestLocation;
}


IndPSO &IndPSO::operator=(const IndPSO &ind) {
    if (this == &ind)
        return *this;

    if (location)
        delete [] location;
	if (velocity)
        delete [] velocity;
	if (particleBestLocation)
	        delete [] particleBestLocation;

    location = new double[IndPSO::Dimension];
	velocity = new double[IndPSO::Dimension];
	particleBestLocation = new double[IndPSO::Dimension];
	
    fitness = ind.fitness;
    fitnessCalculated = ind.fitnessCalculated;
    populationPSO = ind.populationPSO;
    rank = ind.rank;
    linRank = ind.linRank;

    memcpy(location, ind.location, IndPSO::Dimension*sizeof(double));
	memcpy(velocity, ind.velocity, IndPSO::Dimension*sizeof(double));
	memcpy(particleBestLocation, ind.particleBestLocation, IndPSO::Dimension*sizeof(double));

    return *this;
}

bool IndPSO::operator==(const IndPSO &ind) const {
    assert( location != nullptr && ind.location != nullptr );
	assert( velocity != nullptr && ind.velocity != nullptr );
	assert( particleBestLocation != nullptr && ind.particleBestLocation != nullptr );

    if (fitness != ind.fitness
        || fitnessCalculated != ind.fitnessCalculated
        || populationPSO != ind.populationPSO) {
        return false;
    }

    if ( location != nullptr && ind.location != nullptr ) {
        for (int i = 0; i < IndPSO::Dimension; ++i) {
            if ( location[i] != ind.location[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }
	
	
	if ( velocity != nullptr && ind.velocity != nullptr ) {
        for (int i = 0; i < IndPSO::Dimension; ++i) {
            if ( velocity[i] != ind.velocity[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }
	
	if ( particleBestLocation != nullptr && ind.particleBestLocation != nullptr ) {
        for (int i = 0; i < IndPSO::Dimension; ++i) {
            if ( particleBestLocation[i] != ind.particleBestLocation[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }


    return true;
}

void IndPSO::setDimension(size_t dimension) {
    IndPSO::Dimension = dimension;
}

size_t IndPSO::getDimension() {
    return IndPSO::Dimension;
}

double *IndPSO::getLocation() const {
    return location;
}

double *IndPSO::getVelocity() const {
    return velocity;
}

void IndPSO::display() const {
    std::cout << "IndPSO [";
    for (int i = 0; i < IndPSO::Dimension-1; i++)
        std::cout << location[i] << ",";
    std::cout << location[IndPSO::Dimension-1] << "] fitness: " << fitness << std::endl;
}

double IndPSO::getFitness() {
    if (!fitnessCalculated)
        updateFitness();

    return fitness;
}

void IndPSO::updateFitness() {
    fitness = dummy_eval_fitness( location );
    fitnessCalculated = true;
    Algorithm::NumberOfFitnessEvaluation ++;
}


bool IndPSO::isFitnessCalculated() const {
    return fitnessCalculated;
}

void IndPSO::setFitnessCalculated(bool fitnessCalculated) {
    IndPSO::fitnessCalculated = fitnessCalculated;
}

double IndPSO::getEuclideanDistanceToIndividualPSO(const IndPSO &individual) const {
    return Util::GetEuclideanDistancePSO( *this, individual);
}

double IndPSO::getParticleBestFitness(){
	Algorithm::NumberOfFitnessEvaluation ++;
	return dummy_eval_fitness( particleBestLocation );
}

