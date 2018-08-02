/*
 * SwarmManager.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_SWARMMANAGER_H
#define EVOLUTIONARYLIB_SWARMMANAGER_H

#include <vector>
#include <time.h>

#include "BaseSwarm.h"

using std::vector;

class SwarmManager {
protected:
    vector<BaseSwarm*>* swarms;

    // Self-organizing scouts
    timespec SwarmManagerCreationTime;



public:
    SwarmManager(size_t numberOfSwarms,size_t numberOfNeutralParticlesInEachSwarm,size_t numberOfQuantumParticlesInEachSwarm);
    virtual ~SwarmManager();

    //virtual void generateSwarmsRandom();


    vector<BaseSwarm*>* getSwarms();
    BaseSwarm* getSwarm(size_t index) ;
    BaseSwarm* getBestSwarm();
    double getBestSwarmFitness();
    size_t getSwarmManagerSize() const;

    BaseSwarm* createRandomSwarm(size_t numberOfNeutralParticlesInEachSwarm,size_t numberOfQuantumParticlesInEachSwarm);

    virtual void updateFitnesses();
    void setReInitsForEnvironmentChange();
    void triggerFitnessCalculation();
    void testForConvergance();
    void testForExlusion();
    void move();
};

#endif //EVOLUTIONARYLIB_SWARMMANAGER_H
