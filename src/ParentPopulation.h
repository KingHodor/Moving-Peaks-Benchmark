/*
 * ParentPopulation.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_PARENTPOPULATION_H
#define EVOLUTIONARYLIB_PARENTPOPULATION_H


#include "ScoutPopulation.h"

#include <vector>


class ParentPopulation : public Population {
protected:
    std::vector<ScoutPopulation> scoutPopulations;
    size_t suggestedSize;

public:
    ParentPopulation(size_t size);

    std::vector<ScoutPopulation>& getScoutPopulations();

    virtual bool isValidIndividual( const Individual& individual ) const;

    virtual void display() const;

    size_t getTotalSize() const;

    double overallBestFitness();

    void removeScoutPopulation( ScoutPopulation* scout );

    void deleteWorstScoutPopulation();

    void updateAllFitnesses();

    size_t getSuggestedSize() const {
        return suggestedSize;
    }

    void setSuggestedSize(size_t suggestedSize) {
        ParentPopulation::suggestedSize = suggestedSize;
    }

    virtual void generatePopulationRandom();

    void addScoutPopulation( const ScoutPopulation& scoutPopulation );

    void fixIndividuals();

    virtual bool isValidIndividualMovingCenter(const Individual &individual);


};


#endif //EVOLUTIONARYLIB_PARENTPOPULATION_H
