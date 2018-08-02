
#ifndef EVOLUTIONARYLIB_POPULATION_H
#define EVOLUTIONARYLIB_POPULATION_H

#include <vector>
#include <time.h>
#include "Individual.h"

using std::vector;

class Population {
protected:
    vector<Individual> individuals;

    // Self-organizing scouts
    timespec populationCreationTime;
    Individual previousBest;
    double mutationStepSize;


public:
    Population();
    Population(size_t size);
    virtual ~Population();


    virtual void generatePopulationRandom();

    vector<Individual>& getIndividuals();
    Individual& getIndividual(size_t index) ;
    size_t getPopulationSize() const;
    virtual void setIndividuals( const vector<Individual>& inds );
    virtual void setIndividuals( const vector<Individual*>& inds );
    void setIndividual(const Individual& ind, size_t index);
    void setPreviousBest( const Individual& prevBest );
    Individual& getPreviousBest();


    Individual best(int* index = nullptr);
    Individual worst(int* index = nullptr);

    virtual void display() const;

    virtual double getMutationStepSize() const;
    virtual void setMutationStepSize( double size );

    virtual bool isValidIndividual( const Individual& individual ) const;
    virtual bool isValidIndividualMovingCenter( const Individual& individual );

    virtual void updateFitnesses();

    double getDynamism();

    double getFitness();

    virtual void initIndividualRandom( Individual* individual );
    virtual void initIndividualRandomScoutAware( Individual* individual );

    void removeIndividuals( const vector<Individual>& individualsToRemove );
    void removeIndividuals( const vector<Individual*>& individualsToRemove );



    void addIndividuals( const vector<Individual*>& individualsToAdd );
    void addIndividuals( const vector<Individual>& individualsToAdd );
    void addIndividual( const Individual& ind );

    vector<double> getTotalLinRank();

    virtual void replaceIndividualElitist(Individual *individual, Individual *base);

    bool hasRandomIndividual();
    void replaceRandomIndividualOfMemoryWithBestOrElite(const Individual &ind);
    Individual* getClosestIndivual(const Individual &ind);
};


#endif //EVOLUTIONARYLIB_POPULATION_H
