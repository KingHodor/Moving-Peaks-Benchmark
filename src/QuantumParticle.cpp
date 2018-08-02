
#include <string.h>
#include <assert.h>
#include "Algorithm.h"
#include <stdlib.h>

#include "Util.h"
#include "movpeaks.h"

#include "QuantumParticle.h"
#include "BaseParticle.h"

QuantumParticle::~QuantumParticle(){

}

QuantumParticle::QuantumParticle(){

}

void QuantumParticle::move(){
	for(int i=0; i<Algorithm::DimensionSize; i++){
		this->coordinates[i] = this->swarm->getBestOfSwarm()->getBestPersonalFitnessCoordinates()[i]+((2*((double)rand() / (double)(RAND_MAX)) - 1) * Algorithm::QuantumCloudRadius);
	}
	this->setFitnessCalculated(false);
}
