
#ifndef EVOLUTIONARYLIB_SCOUTPOPULATION_H
#define EVOLUTIONARYLIB_SCOUTPOPULATION_H


#include "Population.h"

class ScoutPopulation : public Population {
protected:
    Individual center;
    double diameter;
    Population* parentPopulation;
    size_t suggestedSize;
    double quality;

public:
    ScoutPopulation();
    ScoutPopulation(size_t size);

    const Individual &getCenter() const {
        return center;
    }

    void setCenter(const Individual &center) {
        ScoutPopulation::center = center;
    }

    double getDiameter() const {
        return diameter;
    }

    void setDiameter(double diameter) {
        ScoutPopulation::diameter = diameter;
    }

    Population *getParentPopulation() const {
        return parentPopulation;
    }

    void setParentPopulation(Population *parentPopulation) {
        ScoutPopulation::parentPopulation = parentPopulation;
    }

    virtual double getMutationStepSize() const {
        return mutationStepSize *
                (3 * diameter / 100);
    }

    virtual bool isValidIndividual( const Individual& individual ) const;
    virtual bool isValidIndividualMovingCenter(const Individual &individual);


    double getQuality() const {
        return quality;
    }

    void setQuality(double quality) {
        ScoutPopulation::quality = quality;
    }

    size_t getSuggestedSize() const {
        return suggestedSize;
    }

    void setSuggestedSize(size_t suggestedSize) {
        ScoutPopulation::suggestedSize = suggestedSize;
    }

    virtual void initIndividualRandom( Individual* individual );

    virtual void setIndividuals( const vector<Individual*>& inds );
    virtual void setIndividuals( const vector<Individual>& inds );

};


#endif //EVOLUTIONARYLIB_SCOUTPOPULATION_H
