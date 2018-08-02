/*
 * NeutralParticle.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <assert.h>
#include "Algorithm.h"
#include <stdlib.h>


#include "Util.h"
#include "movpeaks.h"
#include "BaseParticle.h"
#include "NeutralParticle.h"

NeutralParticle::~NeutralParticle(){
	free(this->velocity);
}


NeutralParticle::NeutralParticle(){

}

void NeutralParticle::move(){
	unsigned int* EpsilonForGlobal = (unsigned int *) calloc(Algorithm::DimensionSize, sizeof(unsigned int));
	unsigned int* EpsilonForPersonal = (unsigned int *) calloc(Algorithm::DimensionSize, sizeof(unsigned int));
	for(int i=0; i<Algorithm::DimensionSize; i++){
		EpsilonForGlobal[i] = rand()%2;
		EpsilonForPersonal[i] = rand()%2;
	}
	if(this->velocity==nullptr){
		this->velocity = (double *) calloc(Algorithm::DimensionSize, sizeof(double));
//		this->velocity = new double[Algorithm::DimensionSize];
		for(int i=0; i<Algorithm::DimensionSize; i++)
			this->velocity[i] = (Algorithm::MaxCoordinate - Algorithm::MinCoordinate) * ((double)rand() / (double)(RAND_MAX)) + Algorithm::MinCoordinate;
	}
	for(int i=0; i<Algorithm::DimensionSize; i++){
		this->velocity[i] = ((double)rand() / (double)(RAND_MAX)) * (this->velocity[i] + Algorithm::ControlAttractionOfGlobalBest*EpsilonForGlobal[i]*(this->swarm->getBestOfSwarm()->getBestPersonalFitnessCoordinates()[i] - this->coordinates[i]) + Algorithm::ControlAttractionOfPersonalBest*EpsilonForPersonal[i]*(this->getBestPersonalFitnessCoordinates()[i] - this->coordinates[i]));
		this->coordinates[i] += this->velocity[i];
	}
	this->setFitnessCalculated(false);
	free(EpsilonForGlobal);
	free(EpsilonForPersonal);
}
