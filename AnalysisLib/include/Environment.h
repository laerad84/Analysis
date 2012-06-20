#ifndef ENVIRONMENT__H__
#define ENVIRONMENT__H__
#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>

std::string sumFileDir;
std::string convFileDir;
std::string rawFileDir;
std::string expcalFileDir;
std::string simcalFileDir;

void         PrintEnvironment();
unsigned int GetEnvironment();

#endif
