/*
 * config_reader.c
 *
 *  Created on: Jan 23, 2009
 *      Author: zhmurov
 */
#define MDIS_IO_CONFIG_DEBUG

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "wrapper.h"
#include "configreader.h"
#include <algorithm>

#define BUF_SIZE 4096
#define NAME_LENGTH 100
#define VALUE_LENGTH 2048

int paramCount;
char** paramNames;
char** paramValues;

void parseParametersFile(const char* filename, int argc, char *argv[]){
	FILE* file = safe_fopen(filename, "r");
	if(file != NULL){
		printf("Parsing '%s' parameters file...\n", filename);
	} else {
		DIE("ERROR: Parameters file '%s' can not be found.", filename);
	}
	paramCount = 0;
	char buffer[BUF_SIZE];
	while(fgets(buffer, BUF_SIZE, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s\n", buffer);
			paramCount++;
		}
	}
	paramNames = (char**)malloc(paramCount * sizeof(char*));
	paramValues = (char**)malloc(paramCount * sizeof(char*));
	int i;
	for(i = 0; i < paramCount; i++){
		paramNames[i] = (char*)calloc(NAME_LENGTH,  sizeof(char));
		paramValues[i] = (char*)calloc(VALUE_LENGTH,  sizeof(char));
	}
	rewind(file);
	i = 0;
	while(fgets(buffer, BUF_SIZE, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
			strcpy(paramNames[i], pch);
			pch = strtok(NULL, "#\n");
			strcpy(paramValues[i], pch);
#ifdef MDIS_IO_CONFIG_DEBUG
			printf("%s\t%s\n", paramNames[i], paramValues[i]);
#endif
			i++;
		}
	}
	fclose(file);
    // And now patch parameters according to argc/argv (if we have ones)
    for (int i = 2; i < argc; ++i) { // argv[0] = binary name, argv[1] = config name
        for (int c = 0; argv[i][c] != '\0'; ++c) {
            if (argv[i][c] == '=') {
                char name[NAME_LENGTH], value[VALUE_LENGTH];
                strncpy(name, argv[i], c); name[c] = '\0';
                strncpy(value, &argv[i][c + 1], VALUE_LENGTH);
#ifdef MDIS_IO_CONFIG_DEBUG
                printf("Enforced via argv: %s\t=%s\n", name, value);
#endif
                setParameter(name, value, true);
                break;
            }
        }
    }
}

int addParameter(const char* paramName, const char* paramValue){
#ifdef MDIS_IO_CONFIG_DEBUG
    printf("Creating fictional parameter '%s' with value '%s'\n", paramName, paramValue);
#endif
    paramCount++;
	paramNames = (char**)realloc(paramNames,paramCount*sizeof(char*));
	paramValues = (char**)realloc(paramValues,paramCount*sizeof(char*));
    if (paramNames == NULL || paramValues == NULL) return -1;
	paramNames[paramCount-1] = (char*)calloc(NAME_LENGTH, sizeof(char));
    paramValues[paramCount-1] = (char*)calloc(VALUE_LENGTH, sizeof(char));
    if (paramNames[paramCount-1] == NULL || paramValues[paramCount-1] == NULL) return -1;
    strncpy(paramNames[paramCount-1], paramName, NAME_LENGTH-1);
    strncpy(paramValues[paramCount-1], paramValue, VALUE_LENGTH-1);
    return 0;
}

int setParameter(const char* paramName, const char* paramValue, int force){
#ifdef MDIS_IO_CONFIG_DEBUG
    printf("Setting parameter '%s' to value '%s'\n", paramName, paramValue);
#endif
	for(int i = 0; i < paramCount; i++){
		if(strcmp(paramName, paramNames[i]) == 0){
            strncpy(paramValues[i], paramValue, VALUE_LENGTH);
			return 0;
		}
	}
#ifdef MDIS_IO_CONFIG_DEBUG
    printf("Parameter '%s' not found!\n", paramName);
#endif
    if (force) {
       return addParameter(paramName, paramValue);
    }
    return 1;
}


int getParameter(char* paramValue, const char* paramName, const char* defaultValue, int allowDefault){
	int i;
	for(i = 0; i < paramCount; i++){
		if(strcmp(paramName, paramNames[i]) == 0){
			strcpy(paramValue, paramValues[i]);
#ifdef MDIS_IO_CONFIG_DEBUG
			printf("'%s' = '%s'\n", paramName, paramValue);
#endif
			return 0;
		}
	}
	if(allowDefault){
		strcpy(paramValue, defaultValue);
		printf("Using default value for parameter %s: '%s' = '%s'\n", paramName, paramName, defaultValue);
	} else {
		DIE("Parameter '%s' should be specified in a configuration file.", paramName);
	}
	return 0;
}

int getIntegerParameter(const char* paramName, int defaultValue, int allowDefault){
	char paramValue[VALUE_LENGTH];
	char defaultString[VALUE_LENGTH];
	sprintf(defaultString, "%d", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	int result = atoi(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		DIE("ERROR: Wrong value of %s in a configuration file ('%s'). Should be integer.", paramName, paramValue);
	}
	if(error != 0){
		return 0;
	}
	return result;
}

std::vector<int> getIntegerArrayParameter(const char* paramName){
    char a[VALUE_LENGTH];
    getParameter(a, paramName);
    std::vector<int> ret;
    std::istringstream iss(a);

    // Split string
    int t;
    std::string ts;
    while (iss >> ts) {
        std::transform(ts.begin(), ts.end(), ts.begin(), ::tolower ); // Convert 'ts' to lower case
        if (ts == "to") { // We're working with range
            int from = *(ret.end()-1);
            int to;
            iss >> to;
            for (int i = from + 1; i <= to; ++i)
                ret.push_back(i);
        } else {
            t = atoi(ts.c_str());
            ret.push_back(t);
        }
    }
	return ret;
}

long long int getLongIntegerParameter(const char* paramName, long defaultValue, int allowDefault){
	char paramValue[VALUE_LENGTH];
	char defaultString[VALUE_LENGTH];
	sprintf(defaultString, "%ld", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	long long int result = atol(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		DIE("ERROR: Wrong value of %s in a configuration file ('%s'). Should be long integer.", paramName, paramValue);
	}
	if(error != 0){
		return 0;
	}
	return result;
}

float getFloatParameter(const char* paramName, float defaultValue, int allowDefault){
	char paramValue[VALUE_LENGTH];
	char defaultString[VALUE_LENGTH];
	sprintf(defaultString, "%f", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	float result = atof(paramValue);
	if(result == 0.0 && strcmp(paramValue, "0") != 0 && strcmp(paramValue, "0.0") != 0 && strcmp(paramValue, "0.0f") != 0 &&
			strcmp(paramValue, "0.000000") != 0){
		DIE("ERROR: Wrong value of %s in a configuration file ('%s'). Should be float.", paramName, paramValue);
	}
	if(error != 0){
		return 0.0;
	}
	return result;
}

int getYesNoParameter(const char* paramName, int defaultValue, int allowDefault){
	char paramValue[VALUE_LENGTH];
	char defaultString[VALUE_LENGTH];
	if(defaultValue){
		sprintf(defaultString, "YES");
	} else {
		sprintf(defaultString, "NO");
	}
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	if(error != 0){
		return 0;
	}
	if(strcmp(paramValue, "YES") == 0 || strcmp(paramValue, "Yes") == 0 || strcmp(paramValue, "yes") == 0
			|| strcmp(paramValue, "Y") == 0 || strcmp(paramValue, "y") == 0
			|| strcmp(paramValue, "ON") == 0 || strcmp(paramValue, "On") == 0 || strcmp(paramValue, "on") == 0
			|| strcmp(paramValue, "TRUE") == 0 || strcmp(paramValue, "True") == 0 || strcmp(paramValue, "true") == 0){
		return 1;
	} else {
		return 0;
	}
}

int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ, int allowDefault){
		char paramValue[VALUE_LENGTH];
		char defaultString[VALUE_LENGTH];
		sprintf(defaultString, "%f %f %f", defaultX, defaultY, defaultZ);
        int error = getParameter(paramValue, paramName, defaultString, allowDefault);
        if(error != 0){
			return error;
        }
        char* pch = strtok(paramValue, " \t");
        x[0] = atof(pch);
        pch = strtok(NULL, " \t");
        y[0] = atof(pch);
        pch = strtok(NULL, " \t");
        z[0] = atof(pch);
        return 0;
}

#ifdef CUDA
float4 getFloat3Parameter(const char* paramName) {
    float4 t;
    getVectorParameter(paramName, &t.x, &t.y, &t.z);
    t.w = 0.0f;
    return t;
}
#endif

int getMaskedParameter(char* result, const char* paramName, const char* defaultValue, int allowDefault){
	char parameterMask[VALUE_LENGTH];
	int error = getParameter(parameterMask, paramName, defaultValue, allowDefault);
	char* parameterMaskTrim = strtok(parameterMask, " \t");
	if(error == 0){
		applyMask(result, parameterMaskTrim);
	}
	return error;
}

int getParameter(char* paramValue, const char* paramName, const char* defaulValue){
	return getParameter(paramValue, paramName, defaulValue, 1);
}
int getIntegerParameter(const char* paramName, int defaultValue){
	return getIntegerParameter(paramName, defaultValue, 1);
}
long long int getLongIntegerParameter(const char* paramName, long defaultValue){
	return getLongIntegerParameter(paramName, defaultValue, 1);
}
float getFloatParameter(const char* paramName, float defaulValue){
	return getFloatParameter(paramName, defaulValue, 1);
}
int getYesNoParameter(const char* paramName, int defaultValue){
	return getYesNoParameter(paramName, defaultValue, 1);
}
int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ){
	return getVectorParameter(paramName, x, y, z, defaultX, defaultY, defaultZ, 1);
}


int getMaskedParameter(char* result, const char* paramName, const char* defaultValue){
	return getMaskedParameter(result, paramName, defaultValue, 1);
}


int getParameter(char* paramValue, const char* paramName){
	return getParameter(paramValue, paramName, "", 0);
}
int getIntegerParameter(const char* paramName){
	return getIntegerParameter(paramName, 0, 0);
}
long long int getLongIntegerParameter(const char* paramName){
	return getLongIntegerParameter(paramName, 0, 0);
}
float getFloatParameter(const char* paramName){
	return getFloatParameter(paramName, 0.0f, 0);
}
int getYesNoParameter(const char* paramName){
	return getYesNoParameter(paramName, 0, 0);
}
int getVectorParameter(const char* paramName, float* x, float* y, float* z){
	return getVectorParameter(paramName, x, y, z, 0, 0, 0, 0);
}
int getMaskedParameter(char* result, const char* paramName){
	return getMaskedParameter(result, paramName, "", 0);
}

int getMaskedParameterWithReplacement(char* result, const char* paramName,
		const char* replacementString, const char* stringToReplace){
	char tempString[VALUE_LENGTH];
	int error = getMaskedParameter(tempString, paramName);
	replaceString(result, tempString, replacementString, stringToReplace);
	return error;
}

int getMaskedParameterWithReplacement(char* result, const char* paramName, const char* defaultValue,
		const char* replacementString, const char* stringToReplace){
	char tempString[VALUE_LENGTH];
	int error = getMaskedParameter(tempString, paramName, defaultValue, 1);
    if (error == 0)
    	replaceString(result, tempString, replacementString, stringToReplace);
    else
        strncpy(result, tempString, VALUE_LENGTH);
	return error;
}

void applyMask(char* result, const char* parameterMask){
	char tempstring[100];
	strcpy(tempstring, parameterMask);
	int i;
	for(i = 0; i < paramCount; i++){
		char paramName[NAME_LENGTH];
		sprintf(paramName, "<%s>", paramNames[i]);
		char replacementString[VALUE_LENGTH];
		strcpy(replacementString, paramValues[i]);
		replaceString(result, tempstring, strtok(replacementString, " \t"), paramName);
		strcpy(tempstring, result);
	}
}

void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace){
	//printf("Looking for %s in %s.\n", stringToReplace, initialString);
	int len1, len2;
	//printf("Result string: %s\n", resultString);
	if(strstr(initialString, stringToReplace) != NULL){
		len1 = strlen(initialString) - strlen(strstr(initialString, stringToReplace));
		//printf("len1: %d\n", len1);
		strncpy(resultString, initialString, len1);
		//printf("Result string: %s\n", resultString);
		strncpy(&resultString[len1], replacementString, strlen(replacementString));
		len2 = len1 + strlen(stringToReplace);
		len1 += strlen(replacementString);
		//printf("Result string: %s\n", resultString);

		strcpy(&resultString[len1], &initialString[len2]);
		strcat(resultString, "");
		//printf("Found. Result: %s\n", resultString);
	} else {
		strcpy(resultString, initialString);
	}
}
