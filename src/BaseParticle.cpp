/*
 * BaseParticle.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "BaseParticle.h"
#include "Algorithm.h"

#include "Util.h"
#include "movpeaks.h"


BaseParticle::BaseParticle():
	fitnessCalculated(false),
	bestPersonalFitnessCalculated(false),
	bestPersonalFitness(-101),
	fitness(-100), swarm(nullptr),reinit(false),velocity(nullptr){
		coordinates = new double[Algorithm::DimensionSize];
		bestPersonalFitnessCoordinates = new double[Algorithm::DimensionSize];
		for(int i=0;i<Algorithm::DimensionSize;i++){
			this->coordinates[i] = (Algorithm::MaxCoordinate - Algorithm::MinCoordinate) * ((double)rand() / (double)(RAND_MAX)) + Algorithm::MinCoordinate;
		}
}


BaseParticle::~BaseParticle(){
	if(coordinates)
		delete [] coordinates;
	if(bestPersonalFitnessCoordinates)
		delete [] bestPersonalFitnessCoordinates;
//	if(velocity)
}


void BaseParticle::resetVelocity(){
	free(this->velocity);
	this->velocity=nullptr;
}


void BaseParticle::setBaseSwarm(BaseSwarm *swarm){
	this->swarm = swarm;
}

double *BaseParticle::getCoordinates() const {
    return coordinates;
}

double *BaseParticle::getBestPersonalFitnessCoordinates() const {
    return bestPersonalFitnessCoordinates;
}

bool BaseParticle::isFitnessCalculated() const {
    return fitnessCalculated;
}

void BaseParticle::setFitnessCalculated(bool fitnessCalculated) {
	this->fitnessCalculated = fitnessCalculated;
}

bool BaseParticle::isBestPersonalFitnessCalculated() const {
    return bestPersonalFitnessCalculated;
}

void BaseParticle::setBestPersonalFitnessCalculated(bool bestPersonalFitnessCalculated) {
	this->bestPersonalFitnessCalculated = bestPersonalFitnessCalculated;
}

bool BaseParticle::isReinitSet() const{
	return reinit;
}
void BaseParticle::setReinit(bool reinit){
	this->reinit = reinit;
}

BaseSwarm* BaseParticle::getBaseSwarm(){
	return this->swarm;
}

void BaseParticle::setBestPersonalFitness(double fitness, double* coordinates){
	this->bestPersonalFitness = fitness;
	memcpy(this->bestPersonalFitnessCoordinates, coordinates, Algorithm::DimensionSize*sizeof(double));
	checkAndSetBestOfSwarm();
	this->bestPersonalFitnessCalculated=true;
}

void BaseParticle::checkAndSetBestOfSwarm(){
	if(this->getBaseSwarm()->getBestOfSwarm()==0){
		this->getBaseSwarm()->setBestOfSwarm(this);
	}else{
		if(this->bestPersonalFitness>this->getBaseSwarm()->getBestOfSwarm()->bestPersonalFitness)
			this->getBaseSwarm()->setBestOfSwarm(this);
	}
}



double BaseParticle::getFitness() {
    if (!fitnessCalculated)
        updateFitness();

    return fitness;
}

double BaseParticle::getBestPersonalFitness() {
    if (!bestPersonalFitnessCalculated){
    	if(!fitnessCalculated){
    		this->fitness = dummy_eval_fitness( this->coordinates );
    		this->fitnessCalculated=true;
    	}
    	setBestPersonalFitness(this->fitness,this->coordinates);
    }
    return this->bestPersonalFitness;
}

void BaseParticle::updateFitness(){
	this->fitness = dummy_eval_fitness( this->coordinates );
	//cout << "coordinates : " << this->coordinates[0] << " " << this->coordinates[1] << " " << this->coordinates[2] << " " << this->coordinates[3] << " " << this->coordinates[4] << "\n";
	//cout << "fitness : " << this->fitness << "BestPersonelfitness : " << this->getBestPersonalFitness() << "\n";
	this->fitnessCalculated = true;
	if(this->fitness>this->getBestPersonalFitness()){
    	setBestPersonalFitness(this->fitness,this->coordinates);
    	if(this->swarm!=nullptr){
			if(this->fitness>this->swarm->getFitness()){
				this->swarm->setBestOfSwarm(this);
			}
    	}else{
    		this->swarm->setBestOfSwarm(this);
    	}
	}
}

void BaseParticle::move(){
	cout << "Wrong usage!!!!!!!!!!!\n This instance probably created from BaseParticle class whish is abstract!!! \n\n";
}
