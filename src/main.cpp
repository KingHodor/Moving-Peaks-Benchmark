#include <iostream>

#include "Algorithm.h"
#include "Util.h"
#include "Test.h"

using namespace std;

static const char* filename = "mpbrc";

#ifndef TEST

int main() {
Util::Load_RC_File(filename);

	/**/
	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunSimpleGA();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "SGA \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunHyperMutation();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "HM \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunRandomImmigrants();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "RI \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunMemorySearch();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "MS \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunSelfOrganizingScouts();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "SOS \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;
	/**/

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunEIGA();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "EIGA \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunMemoryBasedImmigrantsGeneticAlgorithm();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "MIGA \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunCPSOR();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "CPSOR \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			//<< " AC: " << Algorithm::statistics.getRunAccuracy()
			//<< " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunAMSO();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "AMSO \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	Algorithm::statistics.clearRunStat();
	for (size_t i = 0; i < Algorithm::RunCount; i++) {
		Algorithm::RunMultiQuantumSwarmAlgorithm();
		Util::LoadBar(i, Algorithm::RunCount);
	}

	cout << "mQSO \tOE: " << Algorithm::statistics.getRunOfflineError()
			<< " Avr.Best: " << Algorithm::statistics.getOfflinePerformance()
			<< " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
			<< " Diversity: " << Algorithm::statistics.getCalculateDiversity()
			<< endl; //<< " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

	return 0;
}

#endif
