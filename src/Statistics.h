
#ifndef EVOLUTIONARYLIB_STATISTICS_H
#define EVOLUTIONARYLIB_STATISTICS_H

#include <stddef.h>
#include <ctime>
#include "PopPSO.h"
#include "Population.h"
#include "ParentPopulation.h"
#include "Swarm.h"
#include "BaseSwarm.h"

#include <vector>

class Statistics {
public:
    Statistics();
    ~Statistics();

    double getAbsoluteRecoveryRate() const;
    double getRunAbsoluteRecoveryRate() const;
    double getCalculatedAbsoluteRecoveryRate() const;

    double getAccuracy() const;
    double getRunAccuracy() const;

    double getBestErrorAtChange() const;
    double getRunBestErrorAtChange() const;

    double getOfflineError() const;
    double getRunOfflineError() const;

    double getOfflinePerformance() const;
    double getRunOfflinePerformance() const;

    double getRunTimespan() const;

    void addStat(double bestFitness);
    void addAbsoluteRecoveryRateStat(double bestFitness);
    void addBestErrorAtChangeStat( double bestFitness );
    void addCalculatedAbsoluteRecoveryRate(double arr);

    void clear();
    void clearAbsoluteRecoveryRateStat();
    void clearRunStat();

    void addRunOfflinePerformance(double performance);
    void addRunOfflineError(double error);
    void addRunAccuracy(double error);
    void addRunBestErrorAtChange(double error);
    void addRunAbsoluteRecoveryRate( double val );
    void addBestOfGeneration( double val );
    void addSwarmDiversity(vector<Swarm> &plst);
    void addTimeSpan( clock_t start );

    void addCalculateDiversity(Population pop, int dimention);
    double addCalculateDiversityInmQSO(BaseSwarm* swarm, int dimention);
    void addCalculateDiversitySOS(ParentPopulation pop, int dimention);
    double addCalculateDiversityCPSOR(PopPSO pop, int dimention);
    double calcDiversityCPSOR(vector<PopPSO> &plst, int dimention);
    double addCalculateDiversityAMSO(Swarm pop, int dimention);
    double calcDiversityAMSO(vector<Swarm> &plst, int dimention);
    void addDiversityCPSOR(double diversity_val);
    void addDiversityAMSO(double diversity_val);
    void addDiversitymQSO(double diversity_val);
    double getCalculateDiversity() const;

private:
    std::vector<double> bestErrorAtChange;
    std::vector<double> accuracy;
    std::vector<double> offlineError;
    std::vector<double> offlinePerformance;
    std::vector<double> absoluteRecoveryRate;
    std::vector<double> calculatedAbsoluteRecoveryRate;

    std::vector<double> runBestErrorAtChange;
    std::vector<double> runOfflineError;
    std::vector<double> runOfflinePerformance;
    std::vector<double> runAccuracy;
    std::vector<double> runAbsoluteRecoveryRate;

    std::vector<double> runTimeSpan;

    std::vector<double> diversity;
    std::vector<double> runbestOfGeneration;

    size_t m_generationCounter;
};


#endif //EVOLUTIONARYLIB_STATISTICS_H
