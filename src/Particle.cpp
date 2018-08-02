/*
 * Particle.cpp
 *
 *  Created on: 31 Mar 2016
 *      Author: alptekin
 */

#include "Particle.h"
#include <string.h>
#include <assert.h>
#include "Population.h"
#include "Algorithm.h"
#include <stdlib.h>
#include "Util.h"
#include "movpeaks.h"

size_t Particle::Dimension = 5;

Particle::Particle()
        : fitnessCalculated(false),fitnessPbestCalculated(false),
          fitness(-100), pBestFitness(-100), swarm(nullptr) {
    location = new double[Particle::Dimension];
    velocity = new double[Particle::Dimension];
    pBestLocation = new double[Particle::Dimension];
}

Particle::Particle(const Particle &particle) {
	location = new double[Particle::Dimension];
	velocity = new double[Particle::Dimension];
	pBestLocation = new double[Particle::Dimension];

	fitness = particle.fitness;
	pBestFitness = particle.pBestFitness;
	fitnessPbestCalculated = particle.fitnessPbestCalculated;
	fitnessCalculated = particle.fitnessCalculated;
	swarm = particle.swarm;

	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		location[i] = particle.location[i];
	}
	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		velocity[i] = particle.velocity[i];
	}
	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		pBestLocation[i] = particle.pBestLocation[i];
	}
}

Particle::~Particle() {
    if (location)
        delete [] location;
	if (velocity)
        delete [] velocity;
	if (pBestLocation)
	        delete [] pBestLocation;
}


Particle &Particle::operator=(const Particle &ind) {
    if (this == &ind)
        return *this;

    if (location)
        delete [] location;
    if (velocity)
           delete [] velocity;
    if (pBestLocation)
                    delete [] pBestLocation;

    location = new double[Particle::Dimension];
    velocity = new double[Particle::Dimension];
	pBestLocation = new double[Particle::Dimension];

    fitness = ind.fitness;
    fitnessCalculated = ind.fitnessCalculated;
	pBestFitness = ind.pBestFitness;
	fitnessPbestCalculated = ind.fitnessPbestCalculated;
    swarm = ind.swarm;

	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		location[i] = ind.location[i];
	}
	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		velocity[i] = ind.velocity[i];
	}
	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		pBestLocation[i] = ind.pBestLocation[i];
	}
    return *this;
}

bool Particle::operator==(const Particle &ind) const {
    assert( location != nullptr && ind.location != nullptr );
    assert( velocity != nullptr && ind.velocity != nullptr );
    assert( pBestLocation != nullptr && ind.pBestLocation != nullptr );

    if (fitness != ind.fitness
        || fitnessCalculated != ind.fitnessCalculated
        || pBestFitness != ind.pBestFitness
        || fitnessPbestCalculated != ind.fitnessPbestCalculated
        || swarm != ind.swarm) {
        return false;
    }

    if ( location != nullptr && ind.location != nullptr ) {
        for (int i = 0; i < (int)Particle::Dimension; ++i) {
            if ( location[i] != ind.location[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }

    if ( velocity != nullptr && ind.velocity != nullptr ) {
        for (int i = 0; i < (int)Particle::Dimension; ++i) {
            if ( velocity[i] != ind.velocity[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }

    if ( pBestLocation != nullptr && ind.pBestLocation != nullptr ) {
        for (int i = 0; i < (int)Particle::Dimension; ++i) {
            if ( pBestLocation[i] != ind.pBestLocation[i] ) {
                return false;
            }
        }
    }
    else {
        return false;
    }
    return true;
}

void Particle::setDimension(size_t dimension) {
    Particle::Dimension = dimension;
}

size_t Particle::getDimension() {
    return Particle::Dimension;
}

double Particle::getFitness() {
	if (!fitnessCalculated)
	{
		updateFitness();
	}
    return fitness;
}

void Particle::updateFitness() {
    fitness = dummy_eval_fitness( location );
    fitnessCalculated = true;
}

bool Particle::isFitnessCalculated() const {
    return fitnessCalculated;
}

void Particle::setFitnessCalculated(bool fitnessCalculated) {
	Particle::fitnessCalculated = fitnessCalculated;
}

double Particle::getPbestFitness() {
	if (!fitnessPbestCalculated)
	{
		updatePbestFitness();
	}
    return pBestFitness;
}

void Particle::updatePbestFitness() {
	pBestFitness = dummy_eval_fitness( pBestLocation );
	fitnessPbestCalculated = true;
}

bool Particle::isFitnessPbestCalculated() const {
    return fitnessPbestCalculated;
}

void Particle::setFitnessPbestCalculated(bool fitnessPbestCalculated) {
	Particle::fitnessPbestCalculated = fitnessPbestCalculated;
}

double Particle::getEuclideanDistanceToParticle(const Particle &particle) const {
    return Util::GetEuclideanParticleDistance( *this, particle );
}

void Particle::clear()
{
	fitness = -100;
	pBestFitness = -100;
}
