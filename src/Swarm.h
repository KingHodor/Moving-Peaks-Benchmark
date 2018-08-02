/*
 * Swarm.h
 *
 *  Created on: 31 Mar 2016
 *      Author: alptekin
 */

#ifndef SWARM_H_
#define SWARM_H_

#include <vector>
#include <time.h>
#include "Particle.h"

using std::vector;

class Swarm {
protected:
    vector<Particle> particles;

public:
	Particle gBest;
	Swarm();
	Swarm(size_t size);
	virtual ~Swarm();

	virtual void initParticleRandom(Particle* particle);
	virtual void generatePopulationRandom();
	Particle&  getParticle(size_t index);
    Particle* getParticleAdress(size_t index) ;
	vector<Particle>& getParticles();
	void addParticle(const Particle& particle);
	void removeParticles(const vector<Particle*> &particlessToRemove);
	Particle bestParticle(int* index = nullptr);
	void updateFitnesses();
	void gBestLearn(Particle *particle);
	void PSO();
	void clear();
	void calculateRadius();
	void calculateSCenter();
	void setRemoveFlag(bool value);
	bool isRemoveFlag();
	size_t getParticleSize();
	double getRadius() const {
		return radius;
	}
	double* getSCenter() const {
		return sCenter;
	}

private:
	double radius;
	double* sCenter;
	bool removeFlag;

};


#endif /* SWARM_H_ */
