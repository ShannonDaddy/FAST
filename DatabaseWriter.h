//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASEWRITER_H
#define THREADEDLAST_DATABASEWRITER_H

#include <pthread.h>

class DatabaseWriter {
		private:
		pthread_t thread;

		static void* threadEntry(void *args);

		void threadFunction();

		void startThread();

		void joinThread();

		DatabaseWriter();
};


#endif //THREADEDLAST_DATABASEWRITER_H
