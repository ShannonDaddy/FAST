//
// Created by david on 16/04/16.
//

#include "DatabaseWriter.h"
#include "prjFiles.hh"

void *DatabaseWriter::threadEntry(void *args) {
    ((DatabaseWriter *) args)->threadFunction();
    return NULL;
}

void DatabaseWriter::threadFunction() {
    /*
    while(true){
        SEM_WAIT(output);

        // Copy out the prj
        prjFiles prj;
        prj.accumulatePrj();

        // Copy out the Multisequence

        // Copy out the SubsetSuffixArray

    }
     */
}

void DatabaseWriter::startThread() {
    pthread_create(&thread, NULL, threadEntry, this);
}

void DatabaseWriter::joinThread() {
    pthread_join(thread, NULL);
}

DatabaseWriter::DatabaseWriter() { }
