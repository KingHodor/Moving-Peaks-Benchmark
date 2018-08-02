/*
 * BaseParticle.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_PARTICLE_H
#define EVOLUTIONARYLIB_PARTICLE_H
#include <stddef.h>

class BaseSwarm;

class BaseParticle {
public:
	BaseParticle();
	virtual ~BaseParticle();

    static void setDimension(size_t dimension);
    static size_t getDimension();
    bool isFitnessCalculated() const;
    void setFitnessCalculated(bool fitnessCalculated);
    bool isBestPersonalFitnessCalculated() const;
    void setBestPersonalFitnessCalculated(bool bestPersonalFitnessCalculated);
    bool isReinitSet() const;
    void setReinit(bool reinit);
    BaseSwarm *getBaseSwarm();
    void setBaseSwarm(BaseSwarm *swarm);
    void resetVelocity();

    double* getCoordinates() const;
    double* getBestPersonalFitnessCoordinates() const;

    double getFitness();
    double getBestPersonalFitness();

    void updateFitness();
    //    virtual void move() = 0;
    virtual void move();




protected:
    double* coordinates;
    double* bestPersonalFitnessCoordinates;
    double fitness;
    double bestPersonalFitness;
    bool fitnessCalculated;
    bool bestPersonalFitnessCalculated;
    double* velocity;
    bool reinit;
    BaseSwarm *swarm;
    void setBestPersonalFitness(double fitness, double* coordinates);
    void checkAndSetBestOfSwarm();
};


#endif /* EVOLUTIONARYLIB_PARTICLE_H */
