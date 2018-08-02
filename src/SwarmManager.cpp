/*
 * SwarmManager.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <float.h>
#include <assert.h>
#include "SwarmManager.h"
#include "Algorithm.h"
#include "Util.h"

#include <vector>
#include <algorithm>
#include <assert.h>


SwarmManager::SwarmManager(size_t numberOfSwarms,size_t numberOfNeutralParticlesInEachSwarm,size_t numberOfQuantumParticlesInEachSwarm){

    swarms = new std::vector<BaseSwarm*>;

    // allocate Particles array
    for (int i = 0; i < numberOfSwarms; i++) {
    	createRandomSwarm(numberOfNeutralParticlesInEachSwarm,numberOfQuantumParticlesInEachSwarm);
    }

    clock_gettime( CLOCK_MONOTONIC, &SwarmManagerCreationTime );
}

SwarmManager::~SwarmManager() {

	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		delete *i;
	}
	swarms->clear();
	free(swarms);
}


BaseSwarm* SwarmManager::createRandomSwarm(size_t numberOfNeutralParticlesInEachSwarm,size_t numberOfQuantumParticlesInEachSwarm){
	BaseSwarm* swarm = new BaseSwarm(numberOfNeutralParticlesInEachSwarm,numberOfQuantumParticlesInEachSwarm);
	swarms->push_back(swarm);
	return swarm;
}

vector<BaseSwarm*>* SwarmManager::getSwarms() {
    return swarms;
}

BaseSwarm* SwarmManager::getSwarm(size_t index) {
    assert( index < swarms->size() && index >= 0 );
    return swarms->at(index);
}

BaseSwarm* SwarmManager::getBestSwarm(){
	int index=0;
	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		double fitness = -DBL_MAX;
		if((*i)->getFitness()>fitness){
			fitness = (*i)->getFitness();
			index = std::distance(swarms->begin(),i);
		}
	}
	return swarms->at(index);
}

double SwarmManager::getBestSwarmFitness(){
	double fitness = -DBL_MAX;
	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		if((*i)->getFitness()>fitness){
			fitness = (*i)->getFitness();
		}
	}
	return fitness;
}

size_t SwarmManager::getSwarmManagerSize() const {
    return swarms->size();
}

void SwarmManager::updateFitnesses() {
	int countOfErasedSwarms=0;
	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end();) {
		if((*i)->isReinitSet()){
			delete *i;
			i=swarms->erase(i);
			countOfErasedSwarms++;
			if(i==swarms->end())
				break;
		}
		else
			++i;
	}
	while(countOfErasedSwarms-->0){
		this->createRandomSwarm(Algorithm::NumberOfNeutralParticlesInEachSwarm,Algorithm::NumberOfQuantumParticlesInEachSwarm);
	}
	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		(*i)->updateFitnesses();
	}
}

void SwarmManager::setReInitsForEnvironmentChange(){
    for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
    	(*i)->cancelAllReInits();
    }
}


void SwarmManager::triggerFitnessCalculation(){
    for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
    	//i->setFitnessCalculated(false);
    	(*i)->resetBestFitnessMemory();
    }
}

void SwarmManager::testForConvergance(){
	if(swarms->size() < 2)
		return;
	for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		// Check if all Swarms are converged. If at least 1 BaseSwarm is not converged then return.
    	if(!(*i)->isConverged())
    		return;
    }
	double lowestFitness = swarms->at(0)->getFitness();
	int lowestFitnessIndex = 0;
	for (vector<BaseSwarm*>::iterator i = ++swarms->begin(); i != swarms->end(); ++i) {
		if((*i)->getFitness() < swarms->at(lowestFitnessIndex)->getFitness()){
			lowestFitness = (*i)->getFitness();
			lowestFitnessIndex = std::distance(swarms->begin(),i);
		}
	}
	//cout << "lowestFitnessIndex : " << lowestFitnessIndex << " lowestFitness : " << lowestFitness << "\n";
	swarms->at(lowestFitnessIndex)->setReinit(true);
}

void SwarmManager::testForExlusion(){
    for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
		if((*i)->isReinitSet())
			continue;
    	for (vector<BaseSwarm*>::iterator j = i+1; j!=swarms->end(); ++j){
    		if(j==swarms->end())
    			break;
    		if((*j)->isReinitSet())
    			continue;
    		double radius=0;
    		for(int k=0;k<Algorithm::DimensionSize;k++){
    			double distanceForCurrentDimension = abs((*i)->getBestOfSwarm()->getBestPersonalFitnessCoordinates()[k] - (*j)->getBestOfSwarm()->getBestPersonalFitnessCoordinates()[k]);
    			radius+=distanceForCurrentDimension*distanceForCurrentDimension;
    		}

    		if(radius < Algorithm::ExclusionRadiusSquare){
        		//cout << "radius between " << std::distance(swarms->begin(),i)+1 << " and " << std::distance(swarms->begin(),j)+1 << " for exclusion is : " << radius << " and Exclusion radius is: " << Algorithm::ExclusionRadiusSquare << "\n";
    			if((*i)->getBestOfSwarm()->getBestPersonalFitness()>(*j)->getBestOfSwarm()->getBestPersonalFitness()){
    				//cout << std::distance(swarms->begin(),j)+1 << " th swarm will be reinitialized!\n";
    				(*j)->setReinit(true);
    			}
    			else{
    				//cout << std::distance(swarms->begin(),i)+1 << " th swarm will be reinitialized!\n\n";
    				(*i)->setReinit(true);
    			}
    		}
    	}
    }
}

void SwarmManager::move(){
    for (vector<BaseSwarm*>::iterator i = swarms->begin(); i != swarms->end(); ++i) {
    	(*i)->move();
    }
}





