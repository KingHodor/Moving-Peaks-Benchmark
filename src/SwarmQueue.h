/*
 * SwarmQueue.h
 *
 *  Created on: 30 Nis 2016
 *      Author: alptekin
 */

#ifndef SWARMQUEUE_H_
#define SWARMQUEUE_H_


typedef struct QueueNode{

	int t;
	int popSize;
    struct QueueNode *next;

} QueueNode;

class SwarmQueue{

    public:

		QueueNode *back;

		QueueNode *front;

		SwarmQueue();

    	virtual ~SwarmQueue();

        void push(int t, int popSize);

        void pop();

        void clear();

        bool isEmpty();

};
#endif /* SWARMQUEUE_H_ */
