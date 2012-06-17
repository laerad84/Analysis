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
static const short E14_CV    = 10;
static const short E14_Laser = 5;
static const short E14_Cosmic= 20;
static const short E14_OEV   = 44;

#endif /* GENERALTYPES_H_ */
