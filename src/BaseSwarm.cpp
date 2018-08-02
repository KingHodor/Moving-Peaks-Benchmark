/*
 * BaseSwarm.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */



#include <float.h>
#include <assert.h>
#include "Algorithm.h"
#include "Util.h"

#include <vector>
#include <algorithm>
#include <assert.h>


BaseSwarm::BaseSwarm():reinit(false),fitnessCalculated(false),convergance(false),bestOfSwarm(nullptr)
{
	this->baseParticles = new std::vector<BaseParticle*>;
}

BaseSwarm::BaseSwarm(size_t numberOfNeutralParticles,size_t numberOfQuantumParticles):
	reinit(false),fitnessCalculated(false),convergance(false),bestOfSwarm(nullptr),countOfNeutralParticles(numberOfNeutralParticles),countOfQuantumParticles(numberOfQuantumParticles){
    // allocate Particles array
	this->baseParticles = new std::vector<BaseParticle*>;
    for (int i = 0; i < numberOfNeutralParticles; i++) {
    	createRandomNeutralParticle();
    }
    for (int i = 0; i < numberOfQuantumParticles; i++) {
    	createRandomQuantumParticle();
    }

//	int test;
//	std::cin >> test;

    clock_gettime( CLOCK_MONOTONIC, &SwarmCreationTime );
}

BaseSwarm::~BaseSwarm() {
	for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
		delete *i;
	}
	baseParticles->clear();
	delete baseParticles;
	//baseParticles->clear();
}

template<typename Base, typename T>
bool BaseSwarm::instanceof(const T *ptr) {
    return dynamic_cast<const Base*>(ptr) != nullptr;
}

BaseParticle* BaseSwarm::createRandomNeutralParticle(){
	NeutralParticle* neutralParticle = new NeutralParticle();
	BaseParticle* particle = neutralParticle;
	particle->setBaseSwarm(this);
	this->baseParticles->push_back(particle);
	return particle;
}

BaseParticle* BaseSwarm::createRandomQuantumParticle(){
	QuantumParticle* neutralParticle = new QuantumParticle();
	BaseParticle *particle = neutralParticle;
	particle->setBaseSwarm(this);
	this->baseParticles->push_back(particle);
	return particle;
}

vector<BaseParticle*>* BaseSwarm::getBaseParticles() {
    return baseParticles;
}

BaseParticle* BaseSwarm::getBaseParticle(size_t index) {
    assert( index < baseParticles->size() && index >= 0 );
    return baseParticles->at(index);
}

size_t BaseSwarm::getSwarmSize() const {
    return baseParticles->size();
}

void BaseSwarm::addBaseParticle(const BaseParticle &particle) {
	BaseParticle i = particle;

	baseParticles->push_back(&i);
}

void BaseSwarm::updateFitnesses() {
	for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
		if((*i)->isReinitSet()){
			(*i)->resetVelocity();
			(*i)->setReinit(false);
		}
        (*i)->updateFitness();
    }
	this->fitnessCalculated = true;
}

double BaseSwarm::getFitness() {
    return bestOfSwarm->getBestPersonalFitness();
}

bool BaseSwarm::isFitnessCalculated() const {
    return fitnessCalculated;
}

void BaseSwarm::setFitnessCalculated(bool fitnessCalculated) {
	this->fitnessCalculated = fitnessCalculated;
    for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
    	(*i)->setFitnessCalculated(fitnessCalculated);
    }
}

void BaseSwarm::resetBestFitnessMemory() {
	this->bestOfSwarm = nullptr;
	this->fitnessCalculated = false;
	this->convergance = false;

    for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
    	(*i)->setBestPersonalFitnessCalculated(false);
    	(*i)->setFitnessCalculated(false);
    }
}


void BaseSwarm::cancelAllReInits(){
	this->setReinit(false);
    for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
    	(*i)->setReinit(true);
    }
}

bool BaseSwarm::isReinitSet() const{
	return reinit;
}
void BaseSwarm::setReinit(bool reinit){
	this->reinit = reinit;
}

bool BaseSwarm::isConverged(){
	bool atLeast2NeutralParticles = false; // If we dont have at least 2 Neutral Particles, then convergance returns always false which means it is skipped.
	for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
		// Check if this particle is Neutral or not. Since we calculate convergance with Neutral particles, others to be skipped
		if(!dynamic_cast<NeutralParticle*>((BaseParticle*)(*i))){
			continue;
		}
    	for (vector<BaseParticle*>::iterator j = i+1; j!=baseParticles->end(); ++j){
    		// Check if i is the last element. If so, there would be no element in iteration j.
    		if(j==baseParticles->end()){
    			break;
    		}
    		// Check if this particle is Neutral or not. Since we calculate convergance with Neutral particles, others to be skipped
    		if(!dynamic_cast<NeutralParticle*>((BaseParticle*)(*j))){
    			continue;
    		}
    		atLeast2NeutralParticles = true;
    		// If any distance in any dimension between any Neutral particles passes ConvergenceRadius, this means it is not converged yet.
    		for(int k=0;k<Algorithm::DimensionSize;k++){
    			double distanceForCurrentDimension = abs((*i)->getCoordinates()[k]-(*j)->getCoordinates()[k]);
    			if(distanceForCurrentDimension>Algorithm::ConvergenceRadius)
    				return false;
    		}
    	}
	}
	if(!atLeast2NeutralParticles)
		return false;
	else
		return true;
}


void BaseSwarm::setBestOfSwarm(BaseParticle* bestOfSwarm ){
	this->bestOfSwarm = bestOfSwarm;
}


BaseParticle* BaseSwarm::getBestOfSwarm(){
	return this->bestOfSwarm;
}


void BaseSwarm::move(){
	for (vector<BaseParticle*>::iterator i = baseParticles->begin(); i != baseParticles->end(); ++i) {
		(*i)->move();
	}
}
