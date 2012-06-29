#ifndef ENVIRONMENT__H__
#define ENVIRONMENT__H__
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

/*
std::string sumFileDir;
std::string convFileDir;
std::string rawFileDir;
std::string expcalFileDir;
std::string simcalFileDir;
*/

char sumFileDir[128];
char convFileDir[128];
char rawFileDir[128];
char simcalFileDir[128];
char expcalFileDir[128];
char waveAnaFileDir[256];

void         PrintEnvironment();
unsigned int GetEnvironment();

#endif
