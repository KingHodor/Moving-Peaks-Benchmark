/*
 * ParentPopulation.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */


#include "ParentPopulation.h"

#include <vector>
#include <iostream>
#include <float.h>

using namespace std;

ParentPopulation::ParentPopulation(size_t size) : Population(size), suggestedSize(size) {

}

vector<ScoutPopulation> &ParentPopulation::getScoutPopulations() {
    return scoutPopulations;
}

bool ParentPopulation::isValidIndividual(const Individual &individual) const {
    for(vector<ScoutPopulation>::const_iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        if ( (*it).getCenter().getEuclideanDistanceToIndividual(individual) < (*it).getDiameter() )
            return false;
    }

    return true;
}

void ParentPopulation::display() const {
    cout << "Parent individuals: " << endl;
    Population::display();
    cout << "Child individuals: " << endl;
    for(vector<ScoutPopulation>::const_iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        (*it).display();
    }
}

size_t ParentPopulation::getTotalSize() const {
    size_t total = getPopulationSize();
    for(vector<ScoutPopulation>::const_iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        total += (*it).getPopulationSize();
    }

    return total;
}

double ParentPopulation::overallBestFitness() {
    double bf;

    // population fitness is the best individual's fitness
    bf = best().getFitness();

    for ( std::vector<ScoutPopulation>::iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        if ( (*it).getFitness() > bf )
            bf = (*it).getFitness();

    }

    return bf;
}

void ParentPopulation::removeScoutPopulation(ScoutPopulation *scout) {
    for( vector<ScoutPopulation>::iterator it = scoutPopulations.begin(); it != scoutPopulations.end();  ) {
        if ( scout == &(*it) ) {
            vector<Individual> individualsToAdd = scout->getIndividuals();
            addIndividuals( individualsToAdd );
            suggestedSize += scout->getSuggestedSize();

            scoutPopulations.erase(it);
        }
        else {
            it++;
        }
    }
}

void ParentPopulation::deleteWorstScoutPopulation() {
    double fworst = DBL_MAX;

    ScoutPopulation* worstScout = nullptr;
    for( vector<ScoutPopulation>::iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        if ( (*it).getFitness() < fworst ) {
            fworst =(*it).getFitness();
            worstScout = &(*it);
        }
    }

    // release individuals to parent
    removeScoutPopulation( worstScout );
}

void ParentPopulation::updateAllFitnesses() {
    for (int j = 0; j < individuals.size(); j++) {
        individuals[j].updateFitness();
    }

    // update child populations' individuals (fitness)
    for ( std::vector<ScoutPopulation>::iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        for (int j = 0; j < (*it).getPopulationSize(); j++) {
            (*it).getIndividual(j).updateFitness();
        }
    }
}

void ParentPopulation::generatePopulationRandom() {
    Population::generatePopulationRandom();
}



void ParentPopulation::addScoutPopulation(const ScoutPopulation &scoutPopulation) {
    scoutPopulations.push_back( scoutPopulation );
}

void ParentPopulation::fixIndividuals() {
    // parent population fixation
    for (int i = 0; i < individuals.size(); ++i) {
        while ( ! isValidIndividual( individuals[i] ) ) {
            initIndividualRandom(&(individuals[i]));
        }
    }

    // child populations fixation
    for(vector<ScoutPopulation>::iterator i = scoutPopulations.begin(); i != scoutPopulations.end(); i++) {
        vector<Individual>& childIndividuals = (*i).getIndividuals();
        for(vector<Individual>::iterator ind = childIndividuals.begin(); ind != childIndividuals.end(); ind++) {
            while ( ! (*i).isValidIndividual(*ind) ) {
                (*i).initIndividualRandom( &(*ind) );
            }
        }
    }
}

bool ParentPopulation::isValidIndividualMovingCenter(const Individual &individual) {
    for(vector<ScoutPopulation>::iterator it = scoutPopulations.begin(); it != scoutPopulations.end(); it++) {
        if ( (*it).best().getEuclideanDistanceToIndividual(individual) < (*it).getDiameter() )
            return false;
    }

    return true;
}
