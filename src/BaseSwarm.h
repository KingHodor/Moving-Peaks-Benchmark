/*
 * BaseSwarm.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_SWARM_H
#define EVOLUTIONARYLIB_SWARM_H

#include <vector>
#include <time.h>
#include "BaseParticle.h"
#include "NeutralParticle.h"
#include "QuantumParticle.h"
#include <typeinfo>

using std::vector;


class BaseSwarm {
protected:
    vector<BaseParticle*>* baseParticles;

    // Self-organizing scouts
    timespec SwarmCreationTime;
    BaseParticle* bestOfSwarm;
    bool reinit;
    bool fitnessCalculated;
    bool convergance;

    template<typename Base, typename T>
    bool instanceof(const T *ptr);


public:
    BaseSwarm();
    BaseSwarm(size_t numberOfNeutralParticlesInEachSwarm,size_t numberOfQuantumParticlesInEachSwarm);
    virtual ~BaseSwarm();

    //virtual void generateSwarmRandom();

    size_t countOfNeutralParticles,countOfQuantumParticles;

    vector<BaseParticle*>* getBaseParticles();
    BaseParticle* getBaseParticle(size_t index) ;
    size_t getSwarmSize() const;
    void addBaseParticle(const BaseParticle &particle);

    void setBestOfSwarm(BaseParticle* bestOfSwarm);
    BaseParticle* getBestOfSwarm();
    BaseParticle* createRandomNeutralParticle();
    BaseParticle* createRandomQuantumParticle();

    bool isFitnessCalculated() const;
    void setFitnessCalculated(bool fitnessCalculated);

    void resetBestFitnessMemory();

    virtual void updateFitnesses();

    double getFitness();
    void cancelAllReInits();

    bool isReinitSet() const;
    void setReinit(bool reinit);

    void move();

    bool isConverged();
    //void setConvergance(bool convergance);

};

#endif //EVOLUTIONARYLIB_SWARM_H
