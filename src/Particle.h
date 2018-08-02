/*
 * Particle.h
 *
 *  Created on: 31 Mar 2016
 *      Author: alptekin
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <stddef.h>

class Swarm;

class Particle {
public:
	double fitness;
	double pBestFitness;
    double* location;
    double* velocity;
    double* pBestLocation;

    Particle();
    Particle(const Particle& ind);
    virtual ~Particle();

    Particle& operator=(const Particle& ind);
    bool operator==(const Particle& ind) const;
    inline bool operator!=(const Particle& ind) const {
        return !(*this == ind);
    }

    static void setDimension(size_t dimension);
    static size_t getDimension();

    bool isFitnessCalculated() const;
    void setFitnessCalculated(bool fitnessCalculated);
    double getFitness();
  	void updateFitness();

    bool isFitnessPbestCalculated() const;
    void setFitnessPbestCalculated(bool fitnessPbestCalculated);
    double getPbestFitness();
  	void updatePbestFitness();

  	double getEuclideanDistanceToParticle(const Particle &particle) const;
  	void clear();

	Swarm *getSwarm() const {
		return swarm;
	}

	void setSwarm(Swarm *swarm) {
		Particle::swarm = swarm;
	}

private:
    static size_t Dimension;
    bool fitnessPbestCalculated;
    bool fitnessCalculated;
    Swarm* swarm;
};

#endif /* PARTICLE_H_ */
