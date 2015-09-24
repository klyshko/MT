/*
 * timer.c
 *
 *  Created on: Mar 23, 2009
 *      Author: zhmurov
 */
#include <time.h>
#include <stdio.h>

clock_t initialTime;
clock_t lastTime;

void printFormatedTime(float timer);

void initTimer(){
	initialTime = clock();
	lastTime = clock();
}

void splitTimer(){
	lastTime = clock();
}

void printTime(long long int step){
	float timer = ((float)(clock() - initialTime))/((float)CLOCKS_PER_SEC);
	printf("Computation time: ");
	printFormatedTime(timer);
	timer = ((float)(clock() - lastTime))/((float)CLOCKS_PER_SEC);
	if(step != 0){
		printf(" (~%f steps/sec)\n", ((float)step)/timer);
	} else {
		printf("\n");
	}
}

void printEstimatedTimeleft(float fractionCompleted){
	if(fractionCompleted != 0.0f){
		float timer = ((float)(clock() - initialTime))/((float)CLOCKS_PER_SEC) * (1.0f/fractionCompleted - 1.0f);
		printf("Estimated time left: ");
		printFormatedTime(timer);
		printf(" (%3.1f%% completed)\n", fractionCompleted*100.0f);
	}
}

void printFormatedTime(float timer){
	int days = (int)(timer/(3600.0f*24.0f));
	int hours = (int)(timer/3600.0f - days*24.0f);
	int minutes = (int)(timer/60.0f - hours*60.0f - days*24.0f*60.0f);
	int seconds = (int)(timer - hours*3600.0f - days*24.0f*3600.0f - minutes*60.0f);
	printf("%dd %dh %dm %ds", days, hours, minutes, seconds);
}

