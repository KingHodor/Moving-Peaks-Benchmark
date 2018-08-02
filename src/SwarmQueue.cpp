/*
 * SwarmQueue.cpp
 *
 *  Created on: 30 Nis 2016
 *      Author: alptekin
 */

#include "SwarmQueue.h"
#include <stddef.h>
#include <iostream>

SwarmQueue::SwarmQueue(){

    back = NULL;

    front = NULL;

}

SwarmQueue::~SwarmQueue() {
	// TODO Auto-generated destructor stub
}

void SwarmQueue::push(int t, int popSize){

    int data;

    QueueNode *temp = new QueueNode;


    temp->t = t;
    temp->popSize = popSize;

    temp->next = NULL;

    if(front == NULL){

        front = temp;

    }else{

        back->next = temp;
    }

    back = temp;

}

void SwarmQueue::pop(){

	QueueNode *temp ;

    if(front == NULL){

       // cout<<"\nQueue is Emtpty\n";

    }else{

        temp = front;

        front = front->next;

        delete temp;

        if(isEmpty())
        {
        	front = back = NULL;        // Reset everything to empty queue
        }

    }
}

void SwarmQueue::clear()
{
	QueueNode *temp;
    while(front != NULL)
    {
        temp = front;
        front = front->next;
        delete temp;
    }

    front = back = NULL; // Reset indices to start over
}

bool SwarmQueue::isEmpty()
{
	 return (front == NULL);
}
