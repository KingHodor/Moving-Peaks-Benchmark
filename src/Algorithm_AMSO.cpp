/*
 * Algorithm_AMSO.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "Algorithm.h"
#include "movpeaks.h"
#include "Particle.h"
#include "Swarm.h"
#include "Util.h"
#include "SwarmQueue.h"

using namespace std;

//#define AMSO_DEBUG_CHANGE
//#define AMSO_DEBUG_EXIT

double **distanceMatrix = NULL;

#define MAX_INDIS 300    /*maximum number of individuals*/
#define MIN_INDIS 70     /*minimum number of individuals*/
#define SUB_SIZE 7       /*maximum individuals in a sub-pop*/
#define B 0.5 			 /*overlapping ratio*/
#define EPSILON 0.0001   /*convergence threshold*/
#define DELTA 1500	     /*trace gap*/
#define ALFA 3           /*population decrease threshold*/
#define STEP 10			 /*population adjustment step size*/
#define MIN_POP_SIZE 5



/* Note: counter is the number of successive increasing points that have the same number of individuals;
 * nextIndis and preIndis are the number of individuals for the next and previous check points, respectively;
 * curPops and prePops are the number of populations in the current and previous increasing points,
 * respectively; step and α are constants with values of 10 and 3, respectively.
 */
int Algorithm::getNextIndis(int preIndis, int curPops, int prePops,
		int *counter) {

	int nextIndis;

	if (*counter == 1) {
		nextIndis = preIndis;
	} else {
		if (curPops - prePops > 0) {
			nextIndis = preIndis + (STEP * (curPops - prePops));
		} else if (prePops - curPops > ALFA) {
			nextIndis = preIndis - (STEP * (prePops - curPops));
		} else {
			nextIndis = preIndis;
		}
	}

	if (nextIndis == preIndis) {
		(*counter) = (*counter) + 1;
		if (curPops > prePops) {
			prePops = curPops;
		}
	} else {
		*counter = 1;
		prePops = curPops;
	}
	return nextIndis;
}


void Algorithm::clustering(Swarm* pop, vector<Swarm> *plst) {

	vector<Swarm> GList; //Create a temporary cluster list G of size |pop|;
	size_t size = pop->getParticleSize();
	size_t gSize, preGSize;

	for (int i = 0; i < (int) size; ++i) {
		Swarm G(1);
		Particle *particle = G.getParticleAdress(0);

		for (int j = 0; j < (int) Algorithm::DimensionSize; j++) {
			particle->location[j] = pop->getParticle(i).location[j];
		}
		for (int k = 0; k < (int) Algorithm::DimensionSize; k++) {
			particle->velocity[k] = pop->getParticle(i).velocity[k];
		}

		GList.push_back(G);
	}

	gSize = GList.size();

	/*allocate distance matrix*/
	allocateDistanceMatrix(gSize);

	/*Calculate the distance between all clusters (i.e., individuals) in G and construct a distance matrix M of size |G| × |G|;*/
	calculateDistanceMatrix(GList);

	while (true) {
		int tIndex;
		int sIndex;
		preGSize = gSize;

		if (!findNearestPair(GList, &tIndex, &sIndex)) {
			break;
		}

		mergeClusters(&GList[tIndex], &GList[sIndex]); // t := t + s; {i.e., merge clusters t and s into t}
		GList.erase(GList.begin() + (sIndex)); // Delete the cluster s from G;

		gSize = GList.size();

		/*clean distance matrix*/
		deallocateDistanceMatrix(preGSize);

		/*allocate distance matrix*/
		allocateDistanceMatrix(gSize);

		/*Recalculate all distances in M which have been affected by the merge of t and s;*/
		calculateDistanceMatrix(GList);

		bool loopContinue = true;
		for (int j = 0; j < (int) gSize; ++j) {
			if (GList[j].getParticleSize() < 2) {
				loopContinue = false;
				break;
			}
		}

		/*if each cluster in G has more than one individual then break*/
		if (loopContinue) {
			break;
		}
	}
	/*clean distance matrix*/
	deallocateDistanceMatrix(gSize);

	/*plst := plst + G;*/
	for (vector<Swarm>::const_iterator it = GList.begin(); it != GList.end();
			it++) {
		Swarm swarm = (*it);
		swarm.calculateRadius();
		plst->push_back(swarm);
	}
}

bool Algorithm::findNearestPair(vector<Swarm> &GList, int *tIndex,
		int *sIndex) {
	size_t size = GList.size();
	bool found = false;
	double mindist = DBL_MAX;
	for (int i = 0; i < (int) (size - 1); ++i) {
		for (int j = i + 1; j < (int) size; ++j) {
			Swarm tCluster = GList[i];
			Swarm sCluster = GList[j];
			size_t totalClusterSize = tCluster.getParticleSize()
					+ sCluster.getParticleSize();
			if (totalClusterSize > SUB_SIZE) {
				continue;
			}
			if (mindist > distanceMatrix[i][j]) {
				mindist = distanceMatrix[i][j];
				*tIndex = i;
				*sIndex = j;
				found = true;
			}

		}
	}
	return found;
}

void Algorithm::mergeClusters(Swarm *tCluster, Swarm *sCluster){
	for(int i=0; i <(int) sCluster->getParticleSize();i++ ){
		tCluster->addParticle(sCluster->getParticle(i));
	}
}


/* Calculate the distance between all clusters (i.e., individuals) in G and
 * construct a distance matrix M of size |G| × |G|;
 */
void Algorithm::calculateDistanceMatrix(vector<Swarm> &GList) {
	double mindist = DBL_MAX;
	double distance;

	size_t size = GList.size();

	for (int l = 0; l < (int) size; ++l) {
		for (int m = 0; m < (int) size; ++m) {
			mindist = DBL_MAX;
			Swarm tCluster = GList[l];
			Swarm sCluster = GList[m];

			int sizeOftCluster = tCluster.getParticleSize();
			int sizeOfsCluster = sCluster.getParticleSize();

			for (int n = 0; n < sizeOftCluster; ++n) {
				for (int o = 0; o < sizeOfsCluster; ++o) {
					distance =
							tCluster.getParticle(n).getEuclideanDistanceToParticle(
									sCluster.getParticle(o));
					if (distance < mindist) {
						mindist = distance;
					}
				}
			}
			distanceMatrix[l][m] = mindist;
		}
	}
}

void Algorithm::mergeAndRemoveOverlappedClusters(vector<Swarm> &plst, Swarm *clst) {

	double rOverlap_t;
	double rOverlap_s;
	double rOverlap;
	for (int i = 0; i < (int) (plst.size() - 1); ++i) {
		for (int j = i + 1; j < (int) plst.size(); ++j) {

			size_t tClusterOverlappingSize = 0;
			size_t sClusterOverlappingSize = 0;

			Swarm *tCluster = &plst[i];
			Swarm *sCluster = &plst[j];

			sCluster->calculateRadius();
			tCluster->calculateRadius();

//			double overLappingDistance = Util::GetEuclideanDistanceByValues(tCluster->getSCenter(), sCluster->getSCenter());
//			double overLappingRadius = tCluster->getRadius() + sCluster->getRadius();

//			if (overLappingRadius > overLappingDistance)
			{
				/*Overlapping occurs*/
				for (int k = 0; k < (int) tCluster->getParticleSize(); ++k) {
					Particle* part = tCluster->getParticleAdress(k);
					double tClusterDistance =
							Util::GetEuclideanDistanceByValues(part->location,
									sCluster->getSCenter());
					if (tClusterDistance < sCluster->getRadius()) {
						tClusterOverlappingSize++;
					}
				}

				for (int l = 0; l < (int) sCluster->getParticleSize(); ++l) {
					Particle* part = sCluster->getParticleAdress(l);
					double tClusterDistance =
							Util::GetEuclideanDistanceByValues(part->location,
									tCluster->getSCenter());
					if (tClusterDistance < tCluster->getRadius()) {
						sClusterOverlappingSize++;
					}
				}

				rOverlap_t = (double)tClusterOverlappingSize
						/ (double)tCluster->getParticleSize();
				rOverlap_s = (double)sClusterOverlappingSize
						/ (double)sCluster->getParticleSize();

				if (rOverlap_t > rOverlap_s) {
					rOverlap = rOverlap_s;
				} else {
					rOverlap = rOverlap_t;
				}

				/*rOverlap (t, s) > β then*/
				if (rOverlap > B) {
					/*Merge t and s into t;*/
					mergeClusters(tCluster, sCluster);
					/*Remove s from plst;*/
					//plst.erase(plst.begin() + j);
					sCluster->setRemoveFlag(true);
				};
			}
		}
	}
}

void Algorithm::overlappingDetection(vector<Swarm> &plst, Swarm *clst){
	if (plst.size() < 2)
			return;


	mergeAndRemoveOverlappedClusters(plst, clst);

	/*Remove worst (|t| − subSize) individuals from t;*//*Over-crowded*/
	for (int i = 0; i < (int) plst.size(); ++i) {
		Swarm *tCluster = &plst[i];
		size_t clusterSize = tCluster->getParticleSize();
		if (clusterSize > SUB_SIZE) {
			size_t worstK = clusterSize - SUB_SIZE;
			vector<Particle>& particles = tCluster->getParticles();

			vector<Particle*> worstParticles = Util::getWorstParticles(worstK,
					particles);
			tCluster->removeParticles(worstParticles);
		}
	}

	for (int j = 0; j < (int) plst.size(); ++j) {

		Swarm *sCluster = &plst[j];
		sCluster->calculateRadius();
		if ((sCluster->getRadius()) < EPSILON) {
			Particle elite = sCluster->bestParticle();
			clst->addParticle(elite);
			//plst.erase(plst.begin() + j);
			sCluster->setRemoveFlag(true);
		}
	}

	for (int k = 0; k < (int) plst.size(); ++k) {
		Swarm *sCluster = &plst[k];

		if ((sCluster)->isRemoveFlag()) {
			plst.erase(plst.begin() + k);
		}

	}
}

void Algorithm::allocateDistanceMatrix(size_t gSize) {
	distanceMatrix = new double *[gSize];

	for (int i = 0; i < (int) gSize; i++)
		distanceMatrix[i] = new double[gSize];
}

void Algorithm::deallocateDistanceMatrix( size_t gSize) {
	for (int i = 0; i < (int) gSize; ++i)
		delete[] distanceMatrix[i];
	delete[] distanceMatrix;
}

double Algorithm::findBestFitnessforAMSO(vector<Swarm> &plst)
{
	double bestFitness = -100;

	for (int k = 0; k < (int) plst.size(); ++k) {
		if (plst[k].gBest.getFitness() > bestFitness) {
			bestFitness = plst[k].gBest.getFitness() ;
		}
	}
	return bestFitness;
}
/*
 * AMSO algorithm employs a single-linkage hierarchical clustering method to
 * generate populations. An overlapping detection scheme is introduced to remove
 * redundant populations during the running process.
 */
void Algorithm::RunAMSO() {
    srand((unsigned int)time(0));
    int changeFreq= Algorithm::ChangePeriod ;
    clock_t start_t, end_t;
    double total_t = 0.0;
    int changeTime= Algorithm::ChangeTime ;
    Algorithm::InitMovPeaks();
    Swarm clst;
    Swarm population( Algorithm::PopulationSize );
    vector<Swarm> plst;

    int counter = 1;
    int pre_indis = 0;
    int s_indis =0;
    int prePops = 0;
    int curPops = 0;
    SwarmQueue Q;

    clock_t start = clock();

    clustering(&population, &plst);
    pre_indis = 0;//Algorithm::PopulationSize;
    //prePops = plst.size();

    Algorithm::statistics.clear();
    start_t = clock();
#ifdef FITNESS_EVALUATION
      while(Algorithm::GenerationCount  > Algorithm::FitnessEvaluationSize){
#else
      while (total_t < Algorithm::TotalTime){
#endif
    	// any change?
#ifdef AMSO_DEBUG_CHANGE
    	/*Population size*/
    	printf("Evaulation size = %d - Population Size = %d\n\r", Algorithm::FitnessEvaluationSize, plst.size());

    	/*Radius*/
    	double average_radius = 0.0;
		for (int i = 0; i < (int) plst.size(); ++i) {
			plst[i].calculateRadius();
			average_radius = average_radius + plst[i].getRadius();
		}
		average_radius = average_radius / double(plst.size());
		printf("Evaulation size = %d - Average Radius = %f\n\r", Algorithm::FitnessEvaluationSize, average_radius);
		/*Q size*/
#endif
#ifdef FITNESS_EVALUATION
    	if (Algorithm::FitnessEvaluationSize > changeFreq) {
#else
    	if (total_t > changeTime) {
#endif
			statistics.addDiversityAMSO(statistics.calcDiversityAMSO(plst, Algorithm::DimensionSize));
			statistics.addStat(findBestFitnessforAMSO(plst));

    		statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            change_peaks();

            for (int j = 0; j < (int)plst.size(); j++) {
            	plst[j].updateFitnesses();
            	plst[j].gBest.updateFitness();
            }

            statistics.addBestErrorAtChangeStat(findBestFitnessforAMSO(plst));

            changeFreq= changeFreq + Algorithm::ChangePeriod;
            changeTime= changeTime + Algorithm::ChangeTime;
    	}

    	for (int i = 0; i < (int)plst.size(); ++i) {
        	plst[i].PSO();
        }

    	overlappingDetection(plst, &clst);

    	s_indis=0;
    	for(vector<Swarm>::iterator it = plst.begin(); it != plst.end();it++){
    		Swarm swarm = (*it);
    		s_indis = s_indis + swarm.getParticleSize();
    	}

    	Q.push( Algorithm::FitnessEvaluationSize, (int) plst.size());

		int period = Q.back->t -  Q.front->t;
		if ((period > DELTA)
				&& ((double) (Q.front->popSize - Q.back->popSize)
						/ (double) period) < 0.002) {

			curPops = plst.size();

			/*Count the number of individuals s_indis;*/
			s_indis = 0;
			for (int i = 0; i < (int)plst.size(); ++i) {
					s_indis = s_indis + plst[i].getParticleSize();
			}
			int e_indis = getNextIndis(pre_indis, curPops, prePops, &counter);


			if (e_indis > MAX_INDIS) {
				e_indis = MAX_INDIS;
			}
			if (e_indis < MIN_INDIS) {
				e_indis = MIN_INDIS;
			}

			int t_indis = e_indis - s_indis - clst.getParticleSize();

			if (t_indis > 0) {

				/*Create a temporal population t pop with t indis random individuals;*/
				Swarm t_pop(t_indis);

				/*Add the individuals of clst into t pop*/
				for (int i = 0; i < (int) clst.getParticleSize(); ++i) {
					Particle particle = clst.getParticle(i);
					t_pop.addParticle(particle);
				}
				clst.clear();
				clustering(&t_pop, &plst);
				Q.clear();

			}
			pre_indis =s_indis;
			prePops = curPops;
		}
		if (period > DELTA) {
			Q.pop();
		}
		statistics.addAbsoluteRecoveryRateStat(findBestFitnessforAMSO(plst));

//#ifdef AMSO_DEBUG_EXIT
		//printf("Evaulation size = %d - gbest fitness= %f - Population Size = %d\n\r", Algorithm::FitnessEvaluationSize, findBestFitnessforAMSO(plst), plst.size());
//#endif
	    end_t = clock();
	    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC * 1000;
    }
    statistics.addTimeSpan( start );
    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );
    Algorithm::ReleaseMovPeaks();
}


