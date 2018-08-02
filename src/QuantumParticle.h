
#ifndef EVOLUTIONARYLIB_QUANTUMPARTICLE_H
#define EVOLUTIONARYLIB_QUANTUMPARTICLE_H


class QuantumParticle : public BaseParticle {
public:
	QuantumParticle();
	virtual ~QuantumParticle();

    virtual void move();

};


#endif /* EVOLUTIONARYLIB_QUANTUMPARTICLE_H */
