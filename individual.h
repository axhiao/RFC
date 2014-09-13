#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"
#include "cec09.h"
#include "common.h"

#define PI  3.1415926535897932384626433832795

class CIndividual{
public:
	CIndividual();
	virtual ~CIndividual();


	vector <double> y_obj;
	int    rank;
	double    crowd_rank;
	vector <double> x_var;
	vector <double> xZ;
	double EQ; //不满足约束的总量

	void   rnd_init();
	void   obj_eval();
	void   show_objective();
	void   show_variable();
	bool   constraint();

    bool   operator<(const CIndividual &ind2);
	bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};

CIndividual::CIndividual()
{
	x_var = vector<double>(nvar, 0);
    y_obj = vector<double>(nobj, 0);
	rank = 0;
	crowd_rank = 1e30;
}

CIndividual::~CIndividual()
{
	x_var.clear();
	y_obj.clear();
}

void CIndividual::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x_var[n] = lowBound[n] + rnd_uni(&rnd_uni_init)*(uppBound[n] - lowBound[n]);    

}

void CIndividual::obj_eval()
{
	/*
	if(!strcmp("UF1", strTestInstance))  CEC09::UF1(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF2", strTestInstance))  CEC09::UF2(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF3", strTestInstance))  CEC09::UF3(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF4", strTestInstance))  CEC09::UF4(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF5", strTestInstance))  CEC09::UF5(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF6", strTestInstance))  CEC09::UF6(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF7", strTestInstance))  CEC09::UF7(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF8", strTestInstance))  CEC09::UF8(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF9", strTestInstance))  CEC09::UF9(x_var, y_obj, x_var.size()); return;
	if(!strcmp("UF10", strTestInstance)) CEC09::UF10(x_var, y_obj, x_var.size()); return;
*/
	if(!strcmp("R2_DTLZ2_M5", strTestInstance))	{ CEC09::R2_DTLZ2_M5(x_var, y_obj,  x_var.size(), y_obj.size()); return;}
	if(!strcmp("R2_DTLZ3_M5", strTestInstance)) {CEC09::R2_DTLZ3_M5(x_var, y_obj,  x_var.size(), y_obj.size()); return;}
	//if(!strcmp("WFG1_M5", strTestInstance))     {CEC09::WFG1_M5(x_var,y_obj,  x_var.size(), y_obj.size()); return;}

	//-------------------DTLZ---------------------//
	if(!strcmp("DTLZ1", strTestInstance)) {DTLZ1(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("DTLZ2", strTestInstance)) {DTLZ2(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("DTLZ3", strTestInstance)) {DTLZ3(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("DTLZ4", strTestInstance)) {DTLZ4(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("DTLZ7", strTestInstance)) {DTLZ7(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	
	//-------------------WFG---------------------//
	if(!strcmp("WFG1", strTestInstance)) {WFG1(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG2", strTestInstance)) {WFG2(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG3", strTestInstance)) {WFG3(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG4", strTestInstance)) {WFG4(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG5", strTestInstance)) {WFG5(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG6", strTestInstance)) {WFG6(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG7", strTestInstance)) {WFG7(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG8", strTestInstance)) {WFG8(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("WFG9", strTestInstance)) {WFG9(x_var, y_obj, x_var.size(), y_obj.size()); return;}

	//---MOPs for Visually Examining Diversity Maintenance Behavior in a Decision Space---//
	if(!strcmp("PStar", strTestInstance)) {PStar(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("LPMS", strTestInstance)) {LPMS(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	
	//---MOPs with a low-dimensional Pareto-optimal front---//
	if(!strcmp("DTLZ5_3", strTestInstance)) {DTLZ5_3(x_var, y_obj, x_var.size(), y_obj.size()); return;}
	if(!strcmp("ankang20100715", strTestInstance)||!strcmp("ankang20001012", strTestInstance)||!strcmp("ankang20030828", strTestInstance)||!strcmp("ankang20051001", strTestInstance))
		{ankang(x_var, xZ, y_obj, EQ); return;}//objective(x_var,xZ, y_obj,EQ);	
}


void CIndividual::show_objective()
{
    for(int n=0; n<nobj; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void CIndividual::show_variable()
{
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[n]);
	printf("\n");
}

void CIndividual::operator=(const CIndividual &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	rank  = ind2.rank;
	crowd_rank = ind2.crowd_rank;
	xZ = ind2.xZ;
}

bool CIndividual::operator<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}


bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}

bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}

bool   CIndividual::constraint()
{
	/*
	if(!strcmp("CF1", strTestInstance))
	{
	}
	if(!strcmp("CF2", strTestInstance))
	{
	}
	if(!strcmp("CF3", strTestInstance))
	{
		double a= 1,N = 2;
		if (y_obj[1]+y_obj[0]*y_obj[0]-a*sin(N * PI * (y_obj[0]*y_obj[0] - y_obj[1] + 1))-1>= 0) 
			return true;
		else
			return false;
	}
	*/
	return true;
}

// defining subproblem 

class CSubproblem  
{
public:
	CSubproblem();
	virtual ~CSubproblem();

	void show();

	CIndividual     indiv;     // best solution
	CIndividual     saved;     // last solution？？？？？？？？？？？？？？？？？？？？？？？？？？？
	vector <double> namda;     // weight vector
	vector <int>    table;     // neighbourhood table
	double			utility;
	double          fitness;

    void  operator=(const CSubproblem &sub2);
};

CSubproblem::CSubproblem()
{
    namda = vector<double>(nobj, 0);
	utility = 1.0;
}

CSubproblem::~CSubproblem()
{
	namda.clear();
	table.clear();
}

void CSubproblem::show()
{
   for(int n=0; n<namda.size(); n++)
   {
       printf("%f ",namda[n]);
   }
   printf("\n");
}

void CSubproblem::operator=(const CSubproblem &sub2)
{
    indiv  = sub2.indiv;
	saved  = sub2.saved;
	table  = sub2.table;
	namda  = sub2.namda;
	utility = sub2.utility;
	fitness = sub2.fitness;
}


#endif
