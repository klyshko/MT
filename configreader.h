/*
 * configreader.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

#pragma once

#include <cstring>
#include <sstream>
#include <vector>
#ifdef CUDA
	#include <vector_types.h> // For float4
#endif
#include "mystl.h"

void parseParametersFile(const char* filename, int argc = 0, char *argv[] = NULL);

int getParameter(char* paramValue, const char* paramName, const char* defaulValue, int allowDefault);
int getIntegerParameter(const char* paramName, int defaultValue, int allowDefault);
long long int getLongIntegerParameter(const char* paramName, long defaultValue, int allowDefault);
float getFloatParameter(const char* paramName, float defaulValue, int allowDefault);
int getYesNoParameter(const char* paramName, int defaultValue, int allowDefault);
int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ, int allowDefault);
int getMaskedParameter(char* result, const char* paramName, const char* defaultValue, int allowDefault);

int addParameter(const char* paramName, const char* paramValue);
int setParameter(const char* paramName, const char* paramValue, int force);

int getParameter(char* paramValue, const char* paramName, const char* defaulValue);
int getIntegerParameter(const char* paramName, int defaultValue);
long long int getLongIntegerParameter(const char* paramName, long defaultValue);
float getFloatParameter(const char* paramName, float defaulValue);
int getYesNoParameter(const char* paramName, int defaultValue);
int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ);
int getMaskedParameter(char* result, const char* paramName, const char* defaultValue);

int getParameter(char* paramValue, const char* paramName);
int getIntegerParameter(const char* paramName);
long long int getLongIntegerParameter(const char* paramName);
float getFloatParameter(const char* paramName);
int getYesNoParameter(const char* paramName);
int getVectorParameter(const char* paramName, float* x, float* y, float* z);
#ifdef CUDA
float4 getFloat3Parameter(const char* paramName);
#endif
int getMaskedParameter(char* result, const char* paramName);

int getMaskedParameterWithReplacement(char* result, const char* paramName,
                const char* replacementString, const char* stringToReplace);
int getMaskedParameterWithReplacement(char* result, const char* paramName, const char* defaultValue,
                const char* replacementString, const char* stringToReplace);

void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);
void applyMask(char* result, const char* parameterMask);

std::vector<int> getIntegerArrayParameter(const char* paramName);

template <typename T>
inline int getMaskedParameterWithReplacementT(char* result, const char* paramName,
                        T replacementString, const char* stringToReplace) {
    return getMaskedParameterWithReplacement(result, paramName, any2str(replacementString).c_str(), stringToReplace);
}

template <typename T1, typename T2>
inline int getMaskedParameterWithReplacementT(char* result, const char* paramName, T1 defaultValue,
                        T2 replacementString, const char* stringToReplace) {
    return getMaskedParameterWithReplacement(result, paramName, any2str(defaultValue).c_str(), any2str(replacementString).c_str(), stringToReplace);
}

template <typename T>
inline int addParameterT(const char* paramName, T paramValue) {
    return addParameter(paramName, any2str(paramValue).c_str());
}

template <typename T>
inline int setParameterT(const char* paramName, T paramValue, int force) {
    return setParameter(paramName, any2str(paramValue).c_str(), force); 
}

//    0
//  _/\\_
//  `/-=|
//  O' `O
template <typename T>
inline T getParameterAs(const char* paramName) {
    char paramValue[4096]; // Just hope it will be enough
    getParameter(paramValue, paramName);
    return str2any<T>(paramValue);
}

template <>
inline std::string getParameterAs<std::string>(const char* paramName) {
    char paramValue[4096]; // Just hope it will be enough
    getParameter(paramValue, paramName);
    std::string ret(paramValue);
    ret.erase( 0                              , ret.find_first_not_of( " \t" ) );
    ret.erase( ret.find_last_not_of( " \t" ) + 1, ret.max_size()               );
    return ret;
}


template <typename T>
inline T getMaskedParameterAs(const char* paramName) {
    char paramValue[4096]; // Just hope it will be enough
    getMaskedParameter(paramValue, paramName);
    return str2any<T>(paramValue);
}

template <>
inline std::string getMaskedParameterAs<std::string>(const char* paramName) {
    char paramValue[4096]; // Just hope it will be enough
    getMaskedParameter(paramValue, paramName);
    std::string ret(paramValue);
    ret.erase( 0                              , ret.find_first_not_of( " \t" ) );
    ret.erase( ret.find_last_not_of( " \t" ) + 1, ret.max_size()               );
    return ret;
}

template <typename T>
inline T getParameterAs(const char* paramName, T defVal) {
    char paramValue[4096]; // Just hope it will be enough
    getParameter(paramValue, paramName, "{default}", 1);
    if (strcmp(paramValue, "{default}") == 0)
        return defVal;
    return str2any<T>(paramValue);
}

inline bool hasParameter(const char* paramName) {
    char paramValue[4096]; // Just hope it will be enough
    getParameter(paramValue, paramName, "{default}", 1);
    if (strcmp(paramValue, "{default}") == 0)
        return false;
    return true;
}
