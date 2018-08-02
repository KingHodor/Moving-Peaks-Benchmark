/*
 * PopPSO.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */


#include "PopPSO.h"

#include <float.h>
#include <assert.h>
#include "Algorithm.h"
#include "Util.h"

#include <vector>
#include <algorithm>
#include <assert.h>

PopPSO::PopPSO()
        : mutationStepSize(3.3) {
	radius = 0.0;
	numberOfIterator = 0;

	Scenter = new double[Algorithm::DimensionSize];
	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		Scenter[i] = 0.0;
	}
}

PopPSO::PopPSO(size_t size)
        : mutationStepSize(3.3) {
    // allocate individuals array

    for (int i = 0; i < size; i++) {
        IndPSO individual;
        initIndividualRandom( &individual );
        individual.setPopulationPSO( this );
        individuals.push_back( individual );
    }

	for (int i = 0; i < (int) Algorithm::DimensionSize; i++) {
		gBest.location[i] = best().location[i];
	}
	gBest.updateFitness();

    clock_gettime( CLOCK_MONOTONIC, &populationCreationTime );
}

PopPSO::~PopPSO() {

}

void PopPSO::generatePopulationRandom() {
    for (vector<IndPSO>::iterator i = individuals.begin(); i != individuals.end(); i++) {
        initIndividualRandom( &(*i) );
    }
}

vector<IndPSO>& PopPSO::getIndividuals() {
    return individuals;
}

IndPSO& PopPSO::getIndividual(size_t index) {
    assert( index < individuals.size() && index >= 0 );
    return individuals[index];
}
IndPSO* PopPSO::getIndividualA(size_t index) {
    assert( index < individuals.size() && index >= 0 );
    return &individuals[index];
}

void PopPSO::setIndividuals(const vector<IndPSO> &inds) {
    individuals = inds;
}

void PopPSO::setIndividuals(const vector<IndPSO *> &inds) {
    individuals.clear();
    for(vector<IndPSO*>::const_iterator it = inds.begin(); it != inds.end(); it++) {
        IndPSO ind;
        ind = *(*it);
        ind.setPopulationPSO( this );
        individuals.push_back( ind );
    }
}

void PopPSO::setIndividual(const IndPSO &ind, size_t index) {
    assert( index < individuals.size() && index >= 0 );
    individuals[index] = ind;
}

IndPSO PopPSO::best(int *index) {
    int bestIndex = 0;
    double bestFitness = -DBL_MAX;
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].getFitness() > bestFitness) {
            bestFitness = individuals[i].getFitness();
            bestIndex = i;
        }
    }

    if (index)
        *index = bestIndex;

    return individuals[bestIndex];
}

IndPSO PopPSO::worst(int *index) {
    int worstIndex = 0;
    double worstFitness = DBL_MAX;
    for (int i = 0; i < individuals.size(); i++) {
        if (individuals[i].getFitness() < worstFitness) {
            worstFitness = individuals[i].getFitness();
            worstIndex = i;
        }
    }

    if (index)
        *index = worstIndex;

    return individuals[worstIndex];
}

size_t PopPSO::getPopulationSize() const {
    return individuals.size();
}

void PopPSO::display() const {
    for (vector<IndPSO>::const_iterator it = individuals.begin(); it != individuals.end(); it++) {
        (*it).display();
    }
}

double PopPSO::getMutationStepSize() const {
    return mutationStepSize;
}

void PopPSO::setMutationStepSize(double size) {
    mutationStepSize = size;
}

bool PopPSO::isValidIndividual( const IndPSO& individual ) const {
    return true;
}

bool PopPSO::isValidIndividualMovingCenter(const IndPSO &individual) {
    return true;
}

void PopPSO::updateFitnesses() {
    for (vector<IndPSO>::iterator i = individuals.begin(); i != individuals.end(); ++i) {
        (*i).updateFitness();
    }
}

void PopPSO::setPreviousBest(const IndPSO &prevBest) {
    previousBest = prevBest;
}

double PopPSO::getFitness() {
    return best().getFitness();
}

double PopPSO::getDynamism() {
    return std::max(0.0, (getFitness()-previousBest.getFitness())/previousBest.getFitness());
}

void PopPSO::initIndividualRandom(IndPSO *individual) {
    double* location = individual->getLocation();
    double* velocity = individual->getVelocity();
    for (int i = 0; i < Algorithm::DimensionSize; i++) {
    	location[i] = Util::RandomDouble(Algorithm::MinCoordinate, Algorithm::MaxCoordinate);
    	//velocity[i] = Util::RandomDouble(Algorithm::MinCoordinate, Algorithm::MaxCoordinate);
    	velocity[i] = 0;
    }

    individual->updateFitness();
}


void PopPSO::initIndividualRandomScoutAware(IndPSO *individual) {
    initIndividualRandom( individual );
}

void PopPSO::removeIndividuals(const vector<IndPSO> &individualsToRemove) {
    for( vector<IndPSO>::const_iterator it = individualsToRemove.begin(); it != individualsToRemove.end(); it++ ) {
        vector<IndPSO>::iterator found = find( individuals.begin(), individuals.end(), *it );
        assert( found != individuals.end() );
        individuals.erase( found );
    }
}

void PopPSO::addIndividuals(const vector<IndPSO *> &individualsToAdd) {
    for( vector<IndPSO*>::const_iterator it = individualsToAdd.begin(); it != individualsToAdd.end(); it++ ) {
        IndPSO ind = *(*it);
        ind.setPopulationPSO( this );
        individuals.push_back( ind );
    }
}

void PopPSO::addIndividuals(const vector<IndPSO> &individualsToAdd) {
    for( vector<IndPSO>::const_iterator it = individualsToAdd.begin(); it != individualsToAdd.end(); it++ ) {
        IndPSO ind = *it;
        ind.setPopulationPSO( this );
        individuals.push_back( ind );
    }
}

IndPSO &PopPSO::getPreviousBest() {
    return previousBest;
}


void PopPSO::addIndividual(const IndPSO &ind) {
    IndPSO i = ind;
    i.setPopulationPSO( this );
    individuals.push_back( i );
}

void PopPSO::removeIndividuals(const vector<IndPSO *> &individualsToRemove) {
    for( vector<IndPSO*>::const_iterator it = individualsToRemove.begin(); it != individualsToRemove.end(); it++ ) {
//        for(int i = 0; i < individuals.size(); i++) {
//            if ( individuals[i] == *(*it) ) {
//                std::cout << "FOUND" << std::endl;
//                break;
//            }
//        }
        vector<IndPSO>::iterator found = find( individuals.begin(), individuals.end(), *(*it) );
//        if (found == individuals.end()) {
//            (*it)->display();
//            std::cout << "-------------" << std::endl;
//
//            display();
//        }

        //assert( found != individuals.end() );
        individuals.erase( found );
    }
}


vector<double> PopPSO::getTotalLinRank() {
    double sum = 0;
    vector<double> linRanks;
    for(int i=0;i<individuals.size();i++){
        sum+=individuals[i].getLinRank();
        linRanks.push_back(sum);
    }

    return linRanks;
}



double PopPSO::calculateRadius(){

	double TotalDistance = 0.0;
	//double center[Algorithm::DimensionSize];

	for (int i = 0; i < Algorithm::DimensionSize; ++i) {
		Scenter[i] = 0;
	}

	for (vector<IndPSO>::iterator it = individuals.begin();it != individuals.end(); it++) {
		IndPSO temp = (*it);

		for (int i = 0; i < (int)Algorithm::DimensionSize; ++i) {
			Scenter[i] = Scenter[i] + temp.getLocation()[i];
		}
	}

	for (int j = 0; j < Algorithm::DimensionSize; ++j) {
		Scenter[j] = Scenter[j] / this->getPopulationSize();
	}

	for (vector<IndPSO>::iterator it = individuals.begin();it != individuals.end(); it++) {
		IndPSO temp = (*it);

		TotalDistance = TotalDistance + Util::GetEuclideanDistanceByValues(temp.getLocation(), Scenter);
	}
	radius = TotalDistance / this->getPopulationSize();

	return radius;
}

double PopPSO::getRadius(){
	return this->radius;
}
double* PopPSO::getScenter(){
	return this->Scenter;
}




void PopPSO::PSO(){
	double w = 0.6, n1=1.7, n2=1.7;


	for(int i = 0;i<(int)this->getPopulationSize();i++){
		IndPSO xTemp = getIndividual(i);
		IndPSO *currentInd = getIndividualA(i);
		IndPSO *pBest = getIndividualA(i);

		//pBest->display();

		double r1=Util::RandomDouble(0,1);
		double r2=Util::RandomDouble(0,1);

		//MNA: Update the particle based on (4) and (5)
		double updatedVelocity[Algorithm::DimensionSize];
		double updatedLocation[Algorithm::DimensionSize];


		for(int j=0;j<(int)Algorithm::DimensionSize;j++){
			//printf("%d[%f-%f ,%f]",j,currentInd->location[j],currentInd->velocity[j],currentInd->getFitness());

			r1=Util::RandomDouble(0,1);
			r2=Util::RandomDouble(0,1);

			updatedVelocity[j] = (w)*currentInd->velocity[j] + (n1)*(r1)*(pBest->particleBestLocation[j] - currentInd->location[j]) + (n2)*(r2)*(gBest.location[j] - currentInd->location[j]);
			currentInd->velocity[j] = updatedVelocity[j];

			//printf("\n [%f-%f ,%f]\n",j,currentInd->location[j],currentInd->velocity[j],currentInd->getFitness());
			//adresi gostermediginden boyle yapiyor
		}
		//printf("\n");


		for(int j=0;j<(int)Algorithm::DimensionSize;j++){
			//printf("L%d[%f-%f ,%f]",j,currentInd->location[j],currentInd->velocity[j],currentInd->getFitness());

			updatedLocation[j]  = currentInd->location[j] + currentInd->velocity[j];

			if(updatedLocation[j] > Algorithm::MaxCoordinate)
				updatedLocation[j] = Algorithm::MaxCoordinate;
			if(updatedLocation[j] < Algorithm::MinCoordinate)
				updatedLocation[j] = Algorithm::MinCoordinate;

			currentInd->location[j] = updatedLocation[j];

			//printf("\nL [%f-%f ,%f]\n",j,currentInd->location[j],currentInd->velocity[j],currentInd->getFitness());

		}
		//printf("\n");

		//printf("\n F:%f, PbestF:%f \n",currentInd->getFitness(),pBest->getParticleBestFitness() );
		currentInd->updateFitness();
		if(currentInd->getFitness() > pBest->getParticleBestFitness()){
			//particleBest update olmasi gerek

			for(int m=0;m<Algorithm::DimensionSize;m++){
				pBest->particleBestLocation[m] = currentInd->location[m];
			}

			if(currentInd->getFitness() > gBest.getFitness()){
				//globalBest update olmasi gerekiyor
				gBest = *currentInd;
				gBest.setFitnessCalculated(false);
				//printf(" >Global Best Degisti< ");
			}
			//printf(" >Particle Best Degisti< ");

			//printf("**** Current F:%f  xTempF: %f * \n",currentInd->getFitness(),xTemp.getFitness());
			if(currentInd->getFitness() > xTemp.getFitness()){
				//gBestLearn fonksiyonuna gitmesi gerekiyor
				//printf("**** Current Daha iyi Oldu ****** \n");
				//printf(" >* gBestF:%f  BEFORE * \n",gBest->getFitness());
				for(int j=0;j<Algorithm::DimensionSize;j++){
					IndPSO t_gBest = gBest;

					t_gBest.location[j] = currentInd->location[j];
					t_gBest.updateFitness();

					if(t_gBest.getFitness() > gBest.getFitness()){
						gBest.location[j] = t_gBest.location[j];
						gBest.setFitnessCalculated(false);
					}

				}
				//printf(" >* gBestF:%f  AFTER * \n",gBest->getFitness());

			}


		}/**/


	}

}
int PopPSO::getNumberOfIterator(){
	return numberOfIterator;
}
void PopPSO::setNumberOfIterator(int val){
	numberOfIterator = val;
}


