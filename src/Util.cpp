/*
 * Util.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */
#include "Util.h"
#include "Algorithm.h"

#include <cmath>
#include <algorithm>
#include <libconfig.h++>
#include <assert.h>
#include <float.h>
#include <stdlib.h>

using namespace libconfig;

using namespace std;

double Util::RandomDouble(double min, double max) {
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}


double Util::gaussRND(double mean, double std) {
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return std * X + mean;
}


bool Util::Load_RC_File(const char *filename) {
    Config cfg;

    cfg.readFile(filename);

    const Setting& root = cfg.getRoot();

    Algorithm::GenerationCount = root["GenerationCount"];
	Algorithm::TotalTime = root["TotalTime"];
	Algorithm::ChangeTime = root["ChangeTime"];
    Algorithm::ChangePeriod = root["ChangePeriod"];
    Algorithm::CrossoverProbability = root["CrossoverProbability"];
    Algorithm::MutationProbability = root["MutationProbability"];
    Algorithm::RandomImmigrantsChangePercentage = root["RandomImmigrantsChangePercentage"];
    Algorithm::PopulationSize = root["PopulationSize"];
    Algorithm::TournamentSize = root["TournamentSize"];
    Algorithm::ExplicitMemorySize = root["ExplicitMemorySize"];
    Algorithm::MemoryUpdateFreq = root["MemoryUpdateFreq"];
    Algorithm::NumberOfPeaks = root["NumberOfPeaks"];
    Algorithm::ShiftLength = root["ShiftLength"];
    Algorithm::SOS_Lambda = root["Lambda"];
    Algorithm::SOS_UseBasisFunction = root["UseBasisFunction"];
    Algorithm::HyperMutationProbability = root["HyperMutationProbability"];
    Algorithm::SOS_ForkingGenerationPeriod = root["ForkingGenerationPeriod"];
    Algorithm::SOS_MinScoutPopulationSizeRelative = root["MinSubPopulationSizeRelative"];
    Algorithm::SOS_MaxScoutPopulationSizeRelative = root["MaxSubPopulationSizeRelative"];
    Algorithm::SOS_MinDiameterRelative = root["MinDiameterRelative"];
    Algorithm::SOS_MaxDiameterRelative = root["MaxDiameterRelative"];
    Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative = root["MinFitnessOfNewForkingPopulationsRelative"];
    Algorithm::SOS_MinFitnessOfExistingForkingPopulationsRelative = root["MinFitnessOfExistingForkingPopulationsRelative"];
    Algorithm::MinCoordinate = root["MinCoord"];
    Algorithm::MaxCoordinate = root["MaxCoord"];
    Algorithm::SOS_MinBasePopulationSizeRelative = root["MinBasePopulationSizeRelative"];
    Algorithm::DimensionSize = root["Dimension"];
    Algorithm::SOS_DiameterReduceFactor = root["SOSDiameterReduceFactor"];
    Algorithm::SOS_ALpha = root["SOSAlpha"];
    Algorithm::RunCount = root["RunCount"];
    Algorithm::EigaImmigrantsRatio = root["EigaImmigrantsRatio"];
	Algorithm::RandomImmigrantsRatioForMIGA = root["RandomImmigrantsRatioForMIGA"];
	Algorithm::MutationRatioOfBestOrEliteIndividual = root["MutationRatioOfBestOrEliteIndividual"];
    Algorithm::NumberOfSwarms = root["NumberOfSwarms"];
	Algorithm::NumberOfNeutralParticlesInEachSwarm = root["NumberOfNeutralParticlesInEachSwarm"];
	Algorithm::NumberOfQuantumParticlesInEachSwarm = root["NumberOfQuantumParticlesInEachSwarm"];
	Algorithm::ControlAttractionOfGlobalBest = root["ControlAttractionOfGlobalBest"];
	Algorithm::ControlAttractionOfPersonalBest = root["ControlAttractionOfPersonalBest"];
	Algorithm::ConvergenceRadius = root["ConvergenceRadius"];
	Algorithm::ExclusionRadius = root["ExclusionRadius"];
	Algorithm::ExclusionRadiusSquare = Algorithm::ExclusionRadius * Algorithm::ExclusionRadius;
	Algorithm::QuantumCloudRadius = root["QuantumCloudRadius"];
    return true;
}

std::vector< Individual* > Util::getWorstIndividuals(size_t worstK, std::vector<Individual> &individuals) {
    vector<Individual*> worstIndividuals;

    double worstFitness;
    double fitness;
    Individual* worstIndividualAtCurrentIteration;

    // get worst k individuals of population
    for (int k = 0; k < worstK; k++) {
        worstFitness = DBL_MAX;

        for (vector<Individual>::iterator i = individuals.begin(); i != individuals.end(); i++) {

            vector<Individual*>::iterator found = find( worstIndividuals.begin(), worstIndividuals.end(), &(*i) );

            // already added to worst individuals array so, skip this element
            if (found != worstIndividuals.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness < worstFitness) {
                worstFitness = fitness;
                worstIndividualAtCurrentIteration = &(*i);
            }
        }

        worstIndividuals.push_back( worstIndividualAtCurrentIteration );
    }

    return worstIndividuals;
}


std::vector< IndPSO* > Util::getWorstIndividualsPSO(size_t worstK, std::vector<IndPSO> &individuals) {
    vector<IndPSO*> worstIndividuals;

    double worstFitness;
    double fitness;
    IndPSO* worstIndividualAtCurrentIteration;

    // get worst k individuals of population
    for (int k = 0; k < worstK; k++) {
        worstFitness = DBL_MAX;

        for (vector<IndPSO>::iterator i = individuals.begin(); i != individuals.end(); i++) {

            vector<IndPSO*>::iterator found = find( worstIndividuals.begin(), worstIndividuals.end(), &(*i) );

            // already added to worst individuals array so, skip this element
            if (found != worstIndividuals.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness < worstFitness) {
                worstFitness = fitness;
                worstIndividualAtCurrentIteration = &(*i);
            }
        }

        worstIndividuals.push_back( worstIndividualAtCurrentIteration );
    }

    return worstIndividuals;
}

std::vector< Particle* > Util::getWorstParticles(size_t worstK, std::vector<Particle> &particles) {
    vector<Particle*> worstParticles;

    double worstFitness;
    double fitness;
    Particle* worstParticleAtCurrentIteration;

    // get worst k Particles of population
    for (int k = 0; k < worstK; k++) {
        worstFitness = DBL_MAX;

        for (vector<Particle>::iterator i = particles.begin(); i != particles.end(); i++) {

            vector<Particle*>::iterator found = find( worstParticles.begin(), worstParticles.end(), &(*i) );

            // already added to worst Particles array so, skip this element
            if (found != worstParticles.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness < worstFitness) {
                worstFitness = fitness;
                worstParticleAtCurrentIteration = &(*i);
            }
        }

        worstParticles.push_back( worstParticleAtCurrentIteration );
    }

    return worstParticles;
}


std::vector< Individual* > Util::getBestIndividuals(size_t bestK, std::vector<Individual> &individuals) {
    vector<Individual*> bestIndividuals;

    double bestFitness;
    size_t bestIndex;
    double fitness;
    Individual* tempBest;

    // get worst k individuals of population
    for (int k = 0; k < bestK; k++) {
        bestFitness = -DBL_MAX;

        for (vector<Individual>::iterator i = individuals.begin();
             i != individuals.end(); i++) {

            vector<Individual*>::iterator found = find( bestIndividuals.begin(), bestIndividuals.end(), &(*i) );

            // already added to worst individuals array so, skip this element
            if (found != bestIndividuals.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness > bestFitness) {
                bestFitness = fitness;
                tempBest = &(*i);
            }
        }

        bestIndividuals.push_back(tempBest);
    }

    return bestIndividuals;
}

double Util::GetEuclideanDistance(const Individual &from, const Individual &to) {
    double* fromValues = from.getValues();
    double* toValues = to.getValues();

    assert( fromValues != nullptr && toValues != nullptr );

    double d = 0;
    for (int i = 0; i < Individual::getDimension(); i++) {
        d += pow(fromValues[i]-toValues[i], 2.0);
    }

    return sqrt( d );
}


double Util::GetEuclideanDistancePSO(const IndPSO &from, const IndPSO &to) {
    double* fromValues = from.getLocation();
    double* toValues = to.getLocation();

    assert( fromValues != nullptr && toValues != nullptr );

    double d = 0;
    for (int i = 0; i < (int)IndPSO::getDimension(); i++) {
        d += pow(fromValues[i]-toValues[i], 2.0);
    }

    return sqrt( d );
}

double Util::GetEuclideanParticleDistance(const Particle &from, const Particle &to) {
    double* fromValues = from.location;
    double* toValues = to.location;

    assert( fromValues != nullptr && toValues != nullptr );

    double d = 0;
    for (int i = 0; i < (int)Particle::getDimension(); i++) {
        d += pow(fromValues[i]-toValues[i], 2.0);
    }

    return sqrt( d );
}

double Util::GetEuclideanDistanceByValues(const double* fromValues , const double* toValues) {

    assert( fromValues != nullptr && toValues != nullptr );

    double d = 0;
    for (int i = 0; i < Individual::getDimension(); i++) {
        d += pow(fromValues[i]-toValues[i], 2.0);
    }

    return sqrt( d );
}

void Util::DisplayPercentage(size_t step, size_t total, PercentageType pt) {
    if (pt == PT_Arrow) {
        cout << "[" << flush;
        for (int i = 0; i < total; i++) {
            if (i == step)
                cout << ">" << flush;
            else if (i < step)
                cout << "=" << flush;
            else
                cout << "." << flush;
        }
        cout << "]" << endl;
    }
    else if (pt == PT_Dot) {
        cout << "." << flush;

        if ( step == total-1 ) {
            cout << endl;
        }
    }
}

void Util::BubbleSort(double *arr, int *sortedIndices, size_t size) {
    int pass,j,hold,in;
    for(pass=1;pass<=size-1;pass++){
        for(j=0;j<=size-2;j++){
            if(arr[j]<arr[j+1]){
                hold = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = hold;
                in = sortedIndices[j];
                sortedIndices[j]=sortedIndices[j+1];
                sortedIndices[j+1]=in;
            }
        }
    }
}

void Util::initIndividual(Individual *individual) {
    double* values = individual->getValues();
    for (int i = 0; i < Algorithm::DimensionSize; i++) {
        values[i] = Util::RandomDouble(Algorithm::MinCoordinate, Algorithm::MaxCoordinate);
    }

    individual->updateFitness();
}
