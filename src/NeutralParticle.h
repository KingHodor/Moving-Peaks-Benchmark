/*
 * NeutralParticle.h
 *
 *  Created on: 8 May 2016
 *      Author: alptekin
 */

#ifndef NEUTRALPARTICLE_H_
#define NEUTRALPARTICLE_H_


class NeutralParticle: public BaseParticle {
public:
	NeutralParticle();
	virtual ~NeutralParticle();
	virtual void move();
};

#endif /* NEUTRALPARTICLE_H_ */
