/*
 * Swarm.cpp
 *
 *  Created on: 31 Mar 2016
 *      Author: alptekin
 */

#include <float.h>
#include <assert.h>
#include "Swarm.h"
#include "Algorithm.h"
#include "Util.h"
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "movpeaks.h"
#include "Particle.h"

#define N1 1.7
#define N2 1.7
#define W_UP 1.0
#define W_LO 0.4

#define W_INIT 0.6
#define MOVE_LIMIT_INIT 0.8

Swarm::Swarm() {
	radius = 0.0;
	sCenter = new double[Algorithm::DimensionSize];
	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		sCenter[i] = 0.0;
	}
	removeFlag = false;
}

Swarm::Swarm(size_t size){
    // allocate particles array
	radius = 0.0;
    for (int i = 0; i < (int)size; i++) {
        Particle Particle;
        initParticleRandom( &Particle );
        Particle.setSwarm( this );
        particles.push_back( Particle );
    }

	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		gBest.location[i] = bestParticle().location[i];
	}
	gBest.updateFitness();

	sCenter = new double[Algorithm::DimensionSize];
	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		sCenter[i] = 0.0;
	}
	radius = 0.0;
	removeFlag = false;
}

Swarm::~Swarm() {

}


void Swarm::generatePopulationRandom() {
    for (vector<Particle>::iterator i = particles.begin(); i != particles.end(); i++) {
        initParticleRandom( &(*i) );
    }
}

size_t Swarm::getParticleSize() {
	return particles.size();
}

void Swarm::initParticleRandom(Particle *particle) {

	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		particle->location[i] = Util::RandomDouble(Algorithm::MinCoordinate,
				Algorithm::MaxCoordinate);
	}
	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		particle->velocity[i] = Util::RandomDouble(Algorithm::MinCoordinate,
				Algorithm::MaxCoordinate);
	}

	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		particle->pBestLocation[i] = Util::RandomDouble(Algorithm::MinCoordinate,
				Algorithm::MaxCoordinate);
	}
}

vector<Particle>& Swarm::getParticles() {
    return particles;
}

Particle& Swarm::getParticle(size_t index) {
    assert( index < particles.size() && index >= 0 );
    return particles[index];
}
Particle* Swarm::getParticleAdress(size_t index) {
    assert( index < particles.size() && index >= 0 );
    return &particles[index];
}

void Swarm::setRemoveFlag(bool value)
{
	removeFlag = value;
}

bool Swarm::isRemoveFlag()
{
	return removeFlag;
}

void Swarm::addParticle(const Particle &particle) {
	Particle i = particle;
	i.setSwarm(this);
	particles.push_back(i);
}

void Swarm::removeParticles(const vector<Particle*> &particlessToRemove) {
	  for( vector<Particle*>::const_iterator it = particlessToRemove.begin(); it != particlessToRemove.end(); it++ ) {
	        vector<Particle>::iterator found = find( particles.begin(), particles.end(), *(*it) );
	        particles.erase( found );
	   }
}

void Swarm::updateFitnesses() {

    for (int i = 0; i < (int)particles.size(); ++i) {
    	Particle * particle = &particles[i];
    	particle->updateFitness();
    	particle->updatePbestFitness();
	}
}

void Swarm::calculateSCenter() {
	for (int i = 0; i < (int) Algorithm::DimensionSize; ++i) {
		sCenter[i] = 0;
	}

	for (vector<Particle>::const_iterator it = particles.begin();
			it != particles.end(); it++) {

		Particle particle = (*it);
		for (int i = 0; i < (int) Algorithm::DimensionSize; ++i) {
			sCenter[i] = sCenter[i] + particle.location[i];
		}
	}

	for (int j = 0; j < (int) Algorithm::DimensionSize; ++j) {
		sCenter[j] = sCenter[j] / (double)getParticleSize();
	}

}
void Swarm::calculateRadius(){

	double totalDistance = 0.0;

	for (int i = 0; i < (int)Algorithm::DimensionSize; ++i) {
		sCenter[i] = 0;
	}

	for (int i = 0; i < (int) Algorithm::DimensionSize; ++i) {
		for (int j = 0; j < (int)particles.size(); ++j) {
			sCenter[i] = sCenter[i] + particles[j].location[i];
		}

	}

	for (int k = 0; k < (int)Algorithm::DimensionSize; ++k) {
		sCenter[k] = sCenter[k] / (double)particles.size();
	}

	for (int l = 0; l < (int)particles.size(); ++l) {
		totalDistance = totalDistance + Util::GetEuclideanDistanceByValues(particles[l].location, sCenter);
	}
	radius = totalDistance / (double)particles.size();
}

Particle Swarm::bestParticle(int *index) {
    int bestIndex = 0;
    double bestFitness = -DBL_MAX;
    for (int i = 0; i < (int)particles.size(); i++) {
        if (particles[i].getFitness() > bestFitness) {
            bestFitness = particles[i].getFitness();
            bestIndex = i;
        }
    }

    if (index)
        *index = bestIndex;

    return particles[bestIndex];
}

void Swarm::clear()
{
	particles.clear();
}

void Swarm::gBestLearn(Particle *particle) {
	double *probabilty = new double[Algorithm::DimensionSize];
	double total = 0.0;
	for (int i = 0; i < (int) Algorithm::DimensionSize; ++i) {
		double * gbestLocation = gBest.location;
		probabilty[i] = fabs(particle->location[i] - gbestLocation[i]);
		total = total + probabilty[i];
	}
	for (int j = 0; j < (int) Algorithm::DimensionSize; ++j) {
		probabilty[j] = probabilty[j] / total;
		if (((double) rand() / RAND_MAX) < probabilty[j]) {
			Particle  temporaryGBest = gBest;

			temporaryGBest.location[j] = particle->location[j];
			temporaryGBest.setFitnessCalculated(false);
			if (temporaryGBest.getFitness() > gBest.getFitness()) {

				gBest.location[j] = temporaryGBest.location[j];
				gBest.setFitnessCalculated(false);
			}
		}
	}
}

void Swarm::PSO(){
		int i, k, j;
		double w = W_INIT;
		double rightmoveLimit = MOVE_LIMIT_INIT * Algorithm::MaxCoordinate;
		double leftmoveLimit = MOVE_LIMIT_INIT * Algorithm::MinCoordinate;

		for (i = 0; i < (int)getParticleSize(); ++i) {
			Particle *particle = getParticleAdress(i);
			Particle tempParticle = *particle;

			/*Update particle*/
			double r1 = Util::RandomDouble(0, 1);
			double r2 = Util::RandomDouble(0, 1);

			//update velocity
			double newVel[Algorithm::DimensionSize];
			for (k = 0; k < (int) Algorithm::DimensionSize; ++k) {
				newVel[k] = (w * particle->velocity[k])
						+ (r1 * N1) * (particle->pBestLocation[k] - particle->location[k])
						+ (r2 * N2) * (gBest.location[k] - particle->location[k]);

				particle->velocity[k] = newVel[k];
			}

			//update location
			double newLoc[Algorithm::DimensionSize];
			for (j = 0; j < (int) Algorithm::DimensionSize; ++j) {
				newLoc[j] = particle->location[j] + newVel[j];
				if (newLoc[j] > Algorithm::MaxCoordinate) {
					newLoc[j] = Algorithm::MaxCoordinate;
				}
				if (newLoc[j] < Algorithm::MinCoordinate) {
					newLoc[j] = newLoc[j] < Algorithm::MinCoordinate;
				}

				particle->location[j] = newLoc[j];
			}
			particle->setFitnessCalculated(false);

			if (particle->getFitness() > particle->getPbestFitness()) {

				for (int i = 0; i < Algorithm::DimensionSize; ++i) {
					particle->pBestLocation[i] = particle->location[i];
				}
				particle->setFitnessPbestCalculated(false);

				if (particle->getFitness() > gBest.getFitness()) {
					gBest = *particle;
					gBest.setFitnessCalculated(false);
				}

				if (particle->getFitness() > tempParticle.getFitness() ) {
					gBestLearn(particle);
				}
			}
		}

}

