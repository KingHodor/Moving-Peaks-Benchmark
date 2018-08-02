
#include "Statistics.h"
#include "time.h"
#include "math.h"
#include "Algorithm.h"
#include <string.h>

extern double global_max;

using namespace std;

Statistics::Statistics() {
}

Statistics::~Statistics() {
}


double Statistics::getOfflineError() const {
    double d = 0.;
    for (vector<double>::const_iterator it = offlineError.begin();
         it != offlineError.end(); it++) {
        d += *it;
    }

    return d / (offlineError.size());
}

double Statistics::getOfflinePerformance() const {
    double d = 0.;

    for (vector<double>::const_iterator it = offlinePerformance.begin();
         it != offlinePerformance.end(); it++) {
        d += *it;
    }

    return d / (offlinePerformance.size());
}

void Statistics::addSwarmDiversity(vector<Swarm> &plst) {
	for (int m = 0; m < plst.size(); ++m) {
		Swarm * swarm = &plst[m];
		int particleSize = swarm->getParticleSize();
		double centroid[particleSize];
		double locations[particleSize][Algorithm::DimensionSize];
		memset(centroid, 0, sizeof(float) * particleSize);
		int i, j, k;
		for (i = 0; i < particleSize; i++) {
			for (j = 0; j < (int) Algorithm::DimensionSize; j++) {
				locations[i][j] = swarm->getParticle(i).location[j];
				centroid[j] += locations[i][j];

			}
		}
		for (j = 0; j < Algorithm::DimensionSize; j++) {
			centroid[j] /= particleSize;
		}

		//then calculate moment of inertia
		float inertia = 0;
		for (k = 0; k < particleSize; k++) {
			for (j = 0; j < (int) Algorithm::DimensionSize; j++) {
				double difference = locations[k][j] - centroid[j];
				inertia += (difference) * (difference);
			}
		}
		diversity.push_back(inertia);
	}
}

void Statistics::addBestOfGeneration( double val )
{
	runbestOfGeneration.push_back(val);
}

void Statistics::addStat(double bestFitness) {

	if (global_max - bestFitness < 0) {
		int g = 5;
	}
    offlineError.push_back( global_max - bestFitness );
    offlinePerformance.push_back( bestFitness );
    accuracy.push_back( bestFitness/global_max );
}

void Statistics::clear() {
    offlineError.clear();
    offlinePerformance.clear();
    accuracy.clear();
    bestErrorAtChange.clear();
    calculatedAbsoluteRecoveryRate.clear();
    absoluteRecoveryRate.clear();


    runTimeSpan.clear();
    diversity.clear();
    Algorithm::FitnessEvaluationSize= 0;
}

void Statistics::addRunOfflinePerformance(double performance) {
    runOfflinePerformance.push_back( performance );
}


void Statistics::addRunOfflineError(double error) {
    runOfflineError.push_back( error );
}

double Statistics::getRunOfflineError() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runOfflineError.begin();
         it != runOfflineError.end(); it++) {
        d += *it;
    }

    return d / (runOfflineError.size());
}


double Statistics::getRunOfflinePerformance() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runOfflinePerformance.begin();
         it != runOfflinePerformance.end(); it++) {
        d += *it;
    }

    return d / (runOfflinePerformance.size());
}

void Statistics::clearRunStat() {
    runOfflinePerformance.clear();
    runOfflineError.clear();
    runBestErrorAtChange.clear();
    runAccuracy.clear();
    runAbsoluteRecoveryRate.clear();

    runTimeSpan.clear();
    diversity.clear();
    Algorithm::NumberOfFitnessEvaluation = 0;
}

double Statistics::getAccuracy() const {
    double a = 0.;
    for (vector<double>::const_iterator it = accuracy.begin(); it != accuracy.end(); it++) {
        a += *it;
    }

    return a / (accuracy.size());
}

double Statistics::getRunAccuracy() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runAccuracy.begin(); it != runAccuracy.end(); it++) {
        a += *it;
    }

    return a / (runAccuracy.size());
}

void Statistics::addRunAccuracy(double error) {
    runAccuracy.push_back( error );
}

void Statistics::addTimeSpan(clock_t start) {
    clock_t end = clock();
    double diffTime = 0.;

    diffTime = ((double)(end - start))/CLOCKS_PER_SEC;
    runTimeSpan.push_back( diffTime );
}

double Statistics::getRunTimespan() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runTimeSpan.begin(); it != runTimeSpan.end(); it++) {
        d += *it;
    }

    return d / (runTimeSpan.size());
}

void Statistics::addRunBestErrorAtChange(double error) {
    runBestErrorAtChange.push_back( error );
}

double Statistics::getBestErrorAtChange() const {
    double a = 0.;
    for (vector<double>::const_iterator it = bestErrorAtChange.begin(); it != bestErrorAtChange.end(); it++) {
        a += *it;
    }

    return a / (bestErrorAtChange.size());
}

double Statistics::getRunBestErrorAtChange() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runBestErrorAtChange.begin(); it != runBestErrorAtChange.end(); it++) {
        a += *it;
    }

    return a / (runBestErrorAtChange.size());
}

void Statistics::addBestErrorAtChangeStat(double bestFitness) {
    bestErrorAtChange.push_back( global_max - bestFitness );
}

double Statistics::getAbsoluteRecoveryRate() const {
	double result ;
    double bestAtGenZero = absoluteRecoveryRate.at(0);
    double bestAtGenLast = absoluteRecoveryRate.at(absoluteRecoveryRate.size()-1);

    double a = 0.;
    for (vector<double>::const_iterator it = absoluteRecoveryRate.begin(); it != absoluteRecoveryRate.end(); it++) {
        a += *it - bestAtGenZero;
    }

    double b = absoluteRecoveryRate.size() * ( global_max - bestAtGenZero );
    //double b = absoluteRecoveryRate.size() * ( bestAtGenLast - bestAtGenZero );


    return a / b;
}



double Statistics::getRunAbsoluteRecoveryRate() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runAbsoluteRecoveryRate.begin(); it != runAbsoluteRecoveryRate.end(); it++) {
        a += *it;
    }

    return a / (runAbsoluteRecoveryRate.size());
}

void Statistics::addRunAbsoluteRecoveryRate(double val) {
    runAbsoluteRecoveryRate.push_back( val );
}

void Statistics::addAbsoluteRecoveryRateStat(double bestFitness) {
    absoluteRecoveryRate.push_back( bestFitness );
}

void Statistics::addCalculatedAbsoluteRecoveryRate(double arr) {
    calculatedAbsoluteRecoveryRate.push_back(arr);
}

double Statistics::getCalculatedAbsoluteRecoveryRate() const {
    double a = 0.;

    for (vector<double>::const_iterator it = calculatedAbsoluteRecoveryRate.begin(); it != calculatedAbsoluteRecoveryRate.end(); it++) {
        a += *it;
    }

    return a / (calculatedAbsoluteRecoveryRate.size());
}

void Statistics::clearAbsoluteRecoveryRateStat() {
    absoluteRecoveryRate.clear();
}

double Statistics::addCalculateDiversityInmQSO(BaseSwarm* swarm, int dimention)
{
	double diversity_val = 0;
	double centroid[dimention];
	double pop_size = 0;
	double Inertia = 0;

	for(int i=0;i<dimention;i++){
		centroid[i]=0;
	}

	for(int n=0;n<swarm->getSwarmSize();n++){
		for(int i=0;i<dimention;i++){
			centroid[i] = centroid[i] + swarm->getBaseParticle(n)->getCoordinates()[i];
		}
		pop_size ++;
	}

	for(int i=0;i<dimention;i++){
		for(int n=0;n<swarm->getSwarmSize();n++){
			Inertia = Inertia + pow(swarm->getBaseParticle(n)->getCoordinates()[i] - (centroid[i]/pop_size) ,2);
		}
	}

	diversity_val = Inertia;


	return diversity_val;
}

void Statistics::addCalculateDiversity(Population pop, int dimention){
	double diversity_val = 0;
	double centroid[dimention];
	double pop_size = 0;
	double Inertia = 0;

	for(int i=0;i<dimention;i++){
		centroid[i]=0;
	}

	for(int n=0;n<pop.getPopulationSize();n++){
		for(int i=0;i<dimention;i++){
			centroid[i] = centroid[i] + pop.getIndividual(n).getValues()[i];
		}
		pop_size ++;
	}

	for(int i=0;i<dimention;i++){
		for(int n=0;n<pop.getPopulationSize();n++){
			Inertia = Inertia + pow(pop.getIndividual(n).getValues()[i] - (centroid[i]/pop_size) ,2);
		}
	}

	diversity_val = Inertia;


	diversity.push_back(diversity_val);
}

void Statistics::addCalculateDiversitySOS(ParentPopulation pop, int dimention){
	double diversity_val = 0;
	double centroid[dimention];
	double pop_size = 0;
	double Inertia = 0;

	for(int i=0;i<dimention;i++){
		centroid[i]=0;
	}

	for(int n=0;n<pop.getPopulationSize();n++){
		for(int i=0;i<dimention;i++){
			centroid[i] = centroid[i] + pop.getIndividual(n).getValues()[i];
		}
		pop_size ++;
	}

	for(int i=0;i<dimention;i++){
		for(int n=0;n<pop.getPopulationSize();n++){
			Inertia = Inertia + pow(pop.getIndividual(n).getValues()[i] - (centroid[i]/pop_size) ,2);
		}
	}

	diversity_val = Inertia;


	diversity.push_back(diversity_val);
}

double Statistics::addCalculateDiversityCPSOR(PopPSO pop, int dimention){
	double diversity_val = 0;
	double centroid[dimention];
	double pop_size = 0;
	double Inertia = 0;

	for(int i=0;i<dimention;i++){
		centroid[i]=0;
	}

	for(int n=0;n<pop.getPopulationSize();n++){
		for(int i=0;i<dimention;i++){
			centroid[i] = centroid[i] + pop.getIndividual(n).location[i];
		}
		pop_size ++;
	}

	for(int i=0;i<dimention;i++){
		for(int n=0;n<pop.getPopulationSize();n++){
			Inertia = Inertia + pow(pop.getIndividual(n).location[i] - (centroid[i]/pop_size) ,2);
		}
	}

	diversity_val = Inertia;


	return diversity_val;
}

double Statistics::calcDiversityCPSOR(vector<PopPSO> &plst, int dimention){

	double diversity_val = 0;
	PopPSO pop;

	int	s_indis=0;
	for (int i = 0; i < plst.size(); ++i) {
		for (int j = 0; j < plst[i].getPopulationSize(); ++j)
		{
			pop.addIndividual(plst[i].getIndividual(j));
		}
	}

	diversity_val = addCalculateDiversityCPSOR(pop, Algorithm::DimensionSize);
	return diversity_val;
}

double Statistics::calcDiversityAMSO(vector<Swarm> &plst, int dimention){

	double diversity_val = 0;
	Swarm pop;

	int	s_indis=0;
	for (int i = 0; i < plst.size(); ++i) {
		for (int j = 0; j < plst[i].getParticleSize(); ++j)
		{
			pop.addParticle(plst[i].getParticle(j));
		}
	}

	diversity_val = addCalculateDiversityAMSO(pop, Algorithm::DimensionSize);
	return diversity_val;
}

double Statistics::addCalculateDiversityAMSO(Swarm pop, int dimention){
	double diversity_val = 0;
	double centroid[dimention];
	double pop_size = 0;
	double Inertia = 0;

	for(int i=0;i<dimention;i++){
		centroid[i]=0;
	}

	for(int n=0;n<pop.getParticleSize();n++){
		for(int i=0;i<dimention;i++){
			centroid[i] = centroid[i] + pop.getParticle(n).location[i];
		}
		pop_size ++;
	}

	for(int i=0;i<dimention;i++){
		for(int n=0;n<pop.getParticleSize();n++){
			Inertia = Inertia + pow(pop.getParticle(n).location[i] - (centroid[i]/pop_size) ,2);
		}
	}

	diversity_val = Inertia;


	return diversity_val;
}

void Statistics::addDiversityCPSOR(double diversity_val){
	diversity.push_back(diversity_val);
}

void Statistics::addDiversityAMSO(double diversity_val){
	diversity.push_back(diversity_val);
}

void Statistics::addDiversitymQSO(double diversity_val){
	diversity.push_back(diversity_val);
}

double Statistics::getCalculateDiversity() const{
    double d = 0.;

    for (vector<double>::const_iterator it = diversity.begin(); it != diversity.end(); it++) {
        d += *it;
    }

    return d / (diversity.size());
}


