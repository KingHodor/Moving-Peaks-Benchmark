/*
 * PopPSO.h
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#ifndef EVOLUTIONARYLIB_POP_H
#define EVOLUTIONARYLIB_POP_H

#include <vector>
#include <time.h>

#include "IndPSO.h"

using std::vector;

class PopPSO {
protected:
    vector<IndPSO> individuals;

    // Self-organizing scouts
    timespec populationCreationTime;
    IndPSO previousBest;
    double mutationStepSize;


public:
    IndPSO gBest;
	double radius;
	double* Scenter;
	vector<IndPSO> pBestIndividuals;
	int numberOfIterator;


    PopPSO();
    PopPSO(size_t size);
    virtual ~PopPSO();


    virtual void generatePopulationRandom();

    vector<IndPSO>& getIndividuals();
    IndPSO& getIndividual(size_t index) ;
    IndPSO* getIndividualA(size_t index) ;
    size_t getPopulationSize() const;
    virtual void setIndividuals( const vector<IndPSO>& inds );
    virtual void setIndividuals( const vector<IndPSO*>& inds );
    void setIndividual(const IndPSO& ind, size_t index);
    void setPreviousBest( const IndPSO& prevBest );
    IndPSO& getPreviousBest();


    IndPSO best(int* index = nullptr);
    IndPSO worst(int* index = nullptr);

    virtual void display() const;

    virtual double getMutationStepSize() const;
    virtual void setMutationStepSize( double size );

    virtual bool isValidIndividual( const IndPSO& individual ) const;
    virtual bool isValidIndividualMovingCenter( const IndPSO& individual );

    virtual void updateFitnesses();

    double getDynamism();

    double getFitness();

    virtual void initIndividualRandom( IndPSO* individual );
    virtual void initIndividualRandomScoutAware( IndPSO* individual );

    void removeIndividuals( const vector<IndPSO>& individualsToRemove );
    void removeIndividuals( const vector<IndPSO*>& individualsToRemove );


    void addIndividuals( const vector<IndPSO*>& individualsToAdd );
    void addIndividuals( const vector<IndPSO>& individualsToAdd );
    void addIndividual( const IndPSO& ind );

    vector<double> getTotalLinRank();
	
	
	double calculateRadius();
	double getRadius();
	double* getScenter();
	void PSO();
	int getNumberOfIterator();
	void setNumberOfIterator(int val);
};


#endif //EVOLUTIONARYLIB_POPULATION_H
