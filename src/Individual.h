/*
 * Individual.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */
#ifndef EVOLUTIONARYLIB_INDIVIDUAL_H
#define EVOLUTIONARYLIB_INDIVIDUAL_H

#include <stddef.h>

class Population;

class Individual {
public:
    Individual();
    Individual(const Individual& ind);
    virtual ~Individual();

    Individual& operator=(const Individual& ind);
    bool operator==(const Individual& ind) const;
    inline bool operator!=(const Individual& ind) const {
        return !(*this == ind);
    }

    static void setDimension(size_t dimension);
    static size_t getDimension();

    static void setHyperMutationEnabled(bool enabled);
    static bool isHyperMutationEnabled();

    double* getValues() const;

    void display() const;

    double getFitness();

    void updateFitness();

    void gaussianMutate(double std = 3.3);

    void EIGAMutate(double std = 3.3);

    bool isFitnessCalculated() const;
    void setFitnessCalculated(bool fitnessCalculated);

    double getEuclideanDistanceToIndividual(const Individual &individual) const;

    double getMutation() const {
        return mutation;
    }

    void setMutation(double mutation) {
        Individual::mutation = mutation;
    }

    Population *getPopulation() const {
        return population;
    }

    void setPopulation(Population *population) {
        Individual::population = population;
    }

    int getRank() const {
        return rank;
    }

    void setRank(int rank) {
        Individual::rank = rank;
    }

    double getLinRank() const {
        return linRank;
    }

    void setLinRank(double linRank) {
        Individual::linRank = linRank;
    }
    bool isRandomlyCreated();

    void setRandomlyCreated(bool randomlyCreated);

    void mutateAccordingToIndividual(const Individual &reference);

private:
    static size_t Dimension;
    static bool HyperMutationEnabled;

    double* values;
    double fitness;
    bool fitnessCalculated;
    double mutation;
    Population* population;
    int rank;
    double linRank;
    bool randomness=true;
};


#endif //EVOLUTIONARYLIB_INDIVIDUAL_H
