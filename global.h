#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>
//#include "object_mxl.h"


using namespace std;

#include "random.h"


#define PI  3.1415926535897932384626433832795

//------------- Parameters in test instance ------------------

int     nvar,  nobj;                    //  the number of variables and objectives

vector<double>  lowBound,   uppBound;   //  lower and upper bounds of variables

char    strTestInstance[256];
char    strFunctionType[256];
int     num_max_evulation;

char	strCrossType[256];

//------------- Parameters in random number ------------------
int     seed    = 177;
long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------

vector <double> idealpoint;
vector <double> nad_point;
double   scale[100];  


int		etax    = 20, 	etam    = 20;   // distribution indexes of crossover and mutation

double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities

int num_of_random =10;
int inst;
int position_parameters = 2;
int distance_parameters = 2;
int size_exter_pop;
int num_update_weight;
double evol_rate;
double start_time_exter_pop_update;
double rate_update_weight;
int preserve_radian = 4;

//------------------Parameters in Light Beam Search of NSGA-II--------------
//------------------Parameters in Light Beam Search of NSGA-II--------------
int numMultipeLightBeams = 1;
vector< vector<double> > AspirationPoint;
vector< vector<double> > ReservationPoint;
vector< vector<double> > VetoThreshold;
vector< vector<double> > Weight;
vector< vector<double> > PreferDirection;
double EpsilonDistance;
//vector<double> CorwdDistance;


//------------------Parameters for ankang reservior flood control--------------
//double  lowBound = 0; //编码下限
//double  uppBound = 37474; //编码上限，枢纽最大下泄流量
double  ZB = 300; //水位下限
double  ZU = 330; //水位上限
double  V_0 = 1737250000;   //初始调度水库容量
double  Q_0 = 773;     //  初始下泄流量
double  ZFL = 325;    // 偏好中心
double D = 2.0;
int  NUM = 3.0;
double ZBP = ZFL - D;     //***********偏好下界*********
double  ZUP = ZFL + D;    //***********偏好上界**********
double  relaxfactor = 1.0;  //删除阈值放松系数
double  relaxfactor1 =  1.2;  //外部种群更新阈值放松系数
double  relaxfactor2 = 1.0;  //增加阈值放松系数
vector<double>  midpoint;
vector<double>  threshold;


int T=6; //调度时间间隔6小时（入库洪水量记录间隔为1小时）

vector <double> Idata; //入库洪水量记录
vector < vector<double> > L2C;//水位库容曲线

double ref_1=330;
double ref_2=37474;
//------------------Parameters for ankang reservior flood control--------------

vector<vector<double>> TOTAL;
vector<vector<double>> HV;
//vector<double> HV ;

void xboundy(void)
{
	int i =0;
	if(!strcmp("UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		return;
	}
	if(!strcmp("UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF8", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF9", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF10", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if (!strcmp("R2_DTLZ2_M5", strTestInstance))
	{
		double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
		double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=0;i<nvar;i++)
		{
			lowBound[i] = low30[i];
			uppBound[i] = upp30[i];
		}
		return;
	}
	if (!strcmp("R2_DTLZ3_M5", strTestInstance))
	{
		double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
		double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=0;i<nvar;i++)
		{
			lowBound[i] = low30[i];
			uppBound[i] = upp30[i];
		}
		return;
	}
	//if (!strcmp("WFG1_M5", strTestInstance))
	//{
	//	lowBound = vector<double>(nvar, 0);
	//	uppBound = vector<double>(nvar, 1);
	//	return;
	//}	
	if (!strcmp("DTLZ1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if(!strcmp("WFG1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG8", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG9", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("PStar", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 100.0);
		return;
	}
	if(!strcmp("LPMS", strTestInstance))
	{
		lowBound = vector<double>(nvar, -4.9);
		uppBound = vector<double>(nvar, 3.2);
		lowBound[1] = -3.5;
		uppBound[1] = 6;
		return;
	}
	if(!strcmp("DTLZ5_3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 1.0);
		return;
	}
	if (!strcmp("ankang20001012", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 37474.0);
		return;
	}
	if (!strcmp("ankang20030828", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 37474.0);
		return;
	}
	if (!strcmp("ankang20051001", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 37474.0);
		return;
	}
	if (!strcmp("ankang20100715", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0.0);
		uppBound = vector<double>(nvar, 37474.0);
		return;
	}
}

#endif