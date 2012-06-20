/*
 * GeneralTypes.h
 *
 *  Created on: Apr 13, 2011
 *      Author: jwlee
 */

#ifndef GENERALTYPES_H_
#define GENERALTYPES_H_
#include <TROOT.h>
#include <CLHEP/Vector/ThreeVector.h>

struct moduleData_t {
	int ID;
	double depE;
	CLHEP::Hep3Vector pos;

};


static const short E14_NCsI  = 2716;
static const short E14_NCC03 = 32;
static const short E14_NCV    = 10;
static const short E14_NLaser = 5;
static const short E14_NCosmic= 20;
static const short E14_NOEV   = 44;
static const short E14_NEtc   = 40;
#endif /* GENERALTYPES_H_ */
