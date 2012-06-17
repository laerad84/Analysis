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

#endif /* GENERALTYPES_H_ */
