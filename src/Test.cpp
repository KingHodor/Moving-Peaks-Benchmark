/*
 * Test.cpp
 *
 *  Created on: 26 Mar 2016
 *      Author: alptekin
 */

#include <iostream>

#include "Test.h"
#include "Algorithm.h"
#include "Util.h"

using namespace std;

static const char* filename = "/home/dev/CLionProjects/evolutionarylib/mpbrc";

#ifdef TEST

int main() {
    Util::Load_RC_File( filename );

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSimpleGA();
        Util::LoadBar( i, Algorithm::RunCount );
        //Util::DisplayPercentage( i, Algorithm::RunCount, Util::PT_Dot );
    }

    cout << "SGA Average offline error: " << Algorithm::statistics.getRunOfflineError()
    <<  " Average running time: " << Algorithm::statistics.getRunOfflineErrorTimespan() <<  endl;


    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunHyperMutation();
        Util::LoadBar( i, Algorithm::RunCount );
        //Util::DisplayPercentage( i, Algorithm::RunCount, Util::PT_Dot );
    }

    cout << "HM Average offline error: " << Algorithm::statistics.getRunOfflineError()
    <<  " Average running time: " << Algorithm::statistics.getRunOfflineErrorTimespan() <<  endl;


    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunRandomImmigrants();
        Util::LoadBar( i, Algorithm::RunCount );
        //Util::DisplayPercentage( i, Algorithm::RunCount, Util::PT_Dot );
    }

    cout << "RI Average offline error: " << Algorithm::statistics.getRunOfflineError()
    <<  " Average running time: " << Algorithm::statistics.getRunOfflineErrorTimespan() <<  endl;


    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunMemorySearch();
        Util::LoadBar( i, Algorithm::RunCount );
        //Util::DisplayPercentage( i, Algorithm::RunCount, Util::PT_Dot );
    }

    cout << "MS Average offline error: " << Algorithm::statistics.getRunOfflineError()
    <<  " Average running time: " << Algorithm::statistics.getRunOfflineErrorTimespan() <<  endl;

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScouts();
        Util::LoadBar( i, Algorithm::RunCount );
//        Util::DisplayPercentage( i, Algorithm::RunCount, Util::PT_Dot );
    }

    cout << "SOS Average offline error: " << Algorithm::statistics.getRunOfflineError()
    <<  " Average running time: " << Algorithm::statistics.getRunOfflineErrorTimespan() <<  endl;


    return 0;
}

#endif
