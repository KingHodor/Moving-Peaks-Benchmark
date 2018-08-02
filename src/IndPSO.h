/*
 * IndPSO.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_IndPSO_H
#define EVOLUTIONARYLIB_IndPSO_H

#include <stddef.h>

class PopPSO;

class IndPSO {
public:

    double* location;
    double* velocity;
    double* particleBestLocation;

    IndPSO();
    IndPSO(const IndPSO& ind);
    virtual ~IndPSO();

    IndPSO& operator=(const IndPSO& ind);
    bool operator==(const IndPSO& ind) const;
    inline bool operator!=(const IndPSO& ind) const {
        return !(*this == ind);
    }

    static void setDimension(size_t dimension);
    static size_t getDimension();

    double* getLocation() const;
	double* getVelocity() const;

    void display() const;

    double getFitness();

    void updateFitness();

    bool isFitnessCalculated() const;
    void setFitnessCalculated(bool fitnessCalculated);

    double getEuclideanDistanceToIndividualPSO(const IndPSO &individualPSO) const;



    PopPSO *getPopulationPSO() const {
        return populationPSO;
    }

    void setPopulationPSO(PopPSO *populationPSO) {
        IndPSO::populationPSO = populationPSO;
    }

    int getRank() const {
        return rank;
    }

    void setRank(int rank) {
        IndPSO::rank = rank;
    }

    double getLinRank() const {
        return linRank;
    }

    void setLinRank(double linRank) {
        IndPSO::linRank = linRank;
    }

    double getParticleBestFitness();


private:
    static size_t Dimension;

    double fitness;
    bool fitnessCalculated;
    PopPSO* populationPSO;
    int rank;
    double linRank;

};


#endif //EVOLUTIONARYLIB_IndividualPSO_H
