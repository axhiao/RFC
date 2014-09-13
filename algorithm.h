#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"
#include "stdio.h"  
#include "math.h" 

class CMOEAD
{

public:
	CMOEAD();
	virtual ~CMOEAD();

	//void init_neighbourhood();                 // calculate the neighbourhood of each subproblem
	void update_neighbour_table();              //calculate the neighbourhood of each subproblem
	void init_population();                    // initialize the population
	void init_prefer_parameter(); 

	void load_parameter();

	void update_reference(CIndividual &ind);					// update ideal point which is used in Tchebycheff or NBI method
	void update_nadpoint();			//？？？？？？？？？？？？？ update nad point which is used in  Tchebycheff , NBI ,PBI method and so on
	void update_problem(CIndividual &ind, int &id, int &type);	// update current solutions in the neighbourhood

	void evol_population();                                      // DE-based recombination？？？？？？？？？像论文中的吗？？？？？？
	void mate_selection(vector<int> &list, int cid, int size, int type);  // select mating parents

	// execute MOEAD
	void exec_emo(int run);

	void save_front(vector <CSubproblem> pop,char savefilename[1024]);       // save the pareto front into files
	void save_exter_pop(char savefilename[1024]);
	void save_pos(vector <CSubproblem> pop,char savefilename[1024]);

	void tour_selection(int depth, vector<int> &selected);
	void comp_utility();
	void external_population_update();
	void organize_merge_split();
	//void use_exter_pop_to_update_evol_pop();
	//vector<double> calculate_crowd_degree(vector<CIndividual> nondonimate_population);
	int calculate_crowd_degree(vector<CSubproblem> nondonimate_population);
	int calculate_crowd_degree(vector<CIndividual> nondonimate_population);
	//int form_indiv_with_min_crowd_degree();
	//int form_indiv_with_min_crowd_degree2();
	void delete_undesired_subproblems();
	void add_spare_subproblems_to_prefer_area();
	int find_sparse_indiv_from_EP();
	double Calculate_HV_for_RFC();
	int Calculate_TOTAL_for_RFC();

    vector <CSubproblem> population;
//	vector <double>      utility;//individual.h  class CSubproblem public:double utility;
	vector <CIndividual> external_population;
	vector <CIndividual> new_population;

	void operator=(const CMOEAD &moea);
	void save_individual(char saveFilename_x[1024],char saveFilename_y[1024],vector<double> x,vector<double> y);
	bool boundy_checkout(CIndividual child);
	void save_HV(vector<vector<double>> hv_value, char saveFilename[1024]);
	void save_TOTAL(vector<vector<double>> total_value, char saveFilename[1024]);
	void save_Z(vector <CSubproblem> pop,char saveFilename[1024]);
public:

	// algorithm parameters
    int		max_gen;       //  the maximal number of generations
	int     pops;          //  the population size
    int	    niche;         //  the neighborhood size
	int     limit;         //  the maximal number of solutions replaced
	double  prob;          //  the neighboring selection probability
	double  rate;          //  the differential rate between solutions

	int     nfes;          //  the number of function evluations

};

CMOEAD::CMOEAD()
{

}

CMOEAD::~CMOEAD()
{
	population.clear();
	//utility.clear();
}

void CMOEAD::init_population()
{
	int i,j;
	idealpoint = vector<double>(nobj, 1.0e+30);
	nad_point = vector<double>(nobj, -1.0e+30);
	//utility    = vector<double>(pops, 1.0);

	char filename[1024];
	// Read weight vectors from a data file
	sprintf(filename,"ParameterSetting/Weight/W%dD_%d.dat", nobj, pops);//?????????
	std::ifstream readf(filename);

    for(i=0; i<pops; i++)
	{
		CSubproblem sub;
		// Randomize and evaluate solution
		sub.indiv.rnd_init();
		sub.indiv.obj_eval();		
		sub.saved = sub.indiv;
		// Initialize the reference point
		update_reference(sub.indiv);
		// Load weight vectors
		for(int j=0; j<nobj; j++)	readf>>sub.namda[j];
		//sub.fitness = fitnessfunction(sub.indiv.y_obj, sub.namda);
		// Save in the population
		population.push_back(sub);//vector <CSubproblem> population; //CSubproblem sub;
		nfes++;//函数评价加1，计数
		//new_population.push_back(sub.indiv);
	}	

	if (!strcmp("_TCH2",strFunctionType)) 
		update_nadpoint();

	readf.close();
}
void CMOEAD::init_prefer_parameter()
{
	if(!strcmp("ankang20001012", strTestInstance))
	{
		midpoint.push_back(328.0);
		midpoint.push_back(6600.0);
		threshold.push_back(4.0);
		threshold.push_back(2000);
	}
	else if (!strcmp("ankang20030828", strTestInstance))
	{
		midpoint.push_back(325.0);
		midpoint.push_back(3800.0);
		threshold.push_back(4.0);
		threshold.push_back(2000);
	}
	else if(!strcmp("ankang20051001", strTestInstance))
	{
		midpoint.push_back(325.5);
		midpoint.push_back(12000.0);
		threshold.push_back(4.0);
		threshold.push_back(2000.0);
	}
	else if(!strcmp("ankang20100715", strTestInstance))
	{
		midpoint.push_back(326.0);
		midpoint.push_back(5500.0);
		threshold.push_back(4.0);
		threshold.push_back(2000.0);
	}
	else {}
}

void CMOEAD::update_neighbour_table()
{
	int i,j,k;
	//1. update the neighbour individual
	vector<double> dist   = vector<double>(population.size(), 0);
	vector<int>    indxxx   = vector<int>(population.size(), 0);

	for (i = 0; i < population.size(); i++)		population[i].table.clear();

	for (i = 0; i < population.size(); i++)
	{
		// calculate the distances based on weight vectors
		for(j=0; j<population.size(); j++)
		{
			dist[j]    = dist_vector(population[i].namda,population[j].namda);
			indxxx[j]  = j;
		}
		// find 'niche' nearest neighboring subproblems
		minfastsort(dist,indxxx,population.size(),niche);
		// save the indexes of the nearest 'niche' neighboring weight vectors
		for(int k=0; k<niche; k++)
			population[i].table.push_back(indxxx[k]);
	}
	
	dist.clear();
	indxxx.clear();
}

void CMOEAD::operator=(const CMOEAD &alg)
{
	//population = alg.population;
}

void CMOEAD::tour_selection(int depth, vector<int> &selected)
{
	// selection based on utility
	vector<int> candidate;
	//for(int k=0;    k<nobj; k++)    selected.push_back(k);   // select first m weights
	for(int n=0; n<pops; n++)    candidate.push_back(n);  // set of unselected weights

	while(selected.size()<int(pops/5.0))
	{
	    int best_idd = int(rnd_uni(&rnd_uni_init)*candidate.size()), i2;
		int best_sub = candidate[best_idd], s2;
		for(int i=1; i<depth; i++)
		{
		    i2  = int(rnd_uni(&rnd_uni_init)*candidate.size());
			s2  = candidate[i2];
			if(population[s2].utility>population[best_sub].utility)
			{
				best_idd = i2;
			    best_sub = s2;
			}
		}
		selected.push_back(best_sub);
		candidate.erase(candidate.begin()+best_idd);
	}

	candidate.clear();
}

void CMOEAD::comp_utility()
{
	double f1, f2, uti, delta;
    for(int n=0; n<pops; n++)
	{
		f1 = fitnessfunction(population[n].indiv.y_obj, population[n].namda);
		f2 = fitnessfunction(population[n].saved.y_obj, population[n].namda);
		delta = f2 - f1;
        if(delta>0.001)  population[n].utility = 1.0;
		else{
/*            uti        = 0.95*(1.0+delta/0.001)*utility[n];
			utility[n] = uti<1.0?uti:1.0;
			*/
			population[n].utility = (0.95+0.05*delta/0.001)*population[n].utility;
		}
        population[n].saved = population[n].indiv;
	}
}

void CMOEAD::update_problem(CIndividual &indiv, int &id, int &type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int i,j,k;
	int size, time = 0;	
	vector<bool> flag=vector<bool>(population.size(), true);
	//---------the subproblems in boundy have poiry------------//
	for (i = 0; i < nobj; i++)
	{
		double f1, f2;
		f1 = fitnessfunction(population[i].indiv.y_obj, population[i].namda);
		f2 = fitnessfunction(indiv.y_obj, population[i].namda);
		if(f2<f1)
		{
			population[i].indiv = indiv;
			flag[i] = false;
			time++;
		}
		if(time>=limit)
		{
			flag.clear();
			return;
		}
	}

    //-----------------------------------------------------------//
	if(type==1)	size = population[id].table.size();  // from neighborhood
	else        size = population.size();            // from whole population

	// a random order to update
	std::vector<int> perm(std::vector<int>(size, 0));
	for(k=0; k<size; k++) perm[k] = k;
	random_shuffle(perm.begin(), perm.end());

    for(i=0; i<size; i++)
	{
		// Pick a subproblem to update
		if(type==1) k = population[id].table[perm[i]];
		else        k = perm[i];

		flag[k] = false;
		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda);
		if(f2<f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if(time>=limit)
		{
			flag.clear();
			perm.clear();
			return;
		}
	}

	//------consider such case: if &indiv can not update its neighbour,while it can update other subpreblems-----------//
	if (type == 1)
	{	// a random order to update
		std::vector<int> perm2(std::vector<int>(population.size(), 0));
		for(k=0; k<size; k++) perm2[k] = k;
		random_shuffle(perm2.begin(), perm2.end());

		for(i=0; i<population.size(); i++)
		{
			int k = perm2[i];
			if (flag[k] == true)
			{
				double f1, f2;
				f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda);
				f2 = fitnessfunction(indiv.y_obj, population[k].namda);
				if(f2<f1)
				{
					population[k].indiv = indiv;
					time++;
				}
				// the maximal number of solutions updated is not allowed to exceed 'limit'
				if(time>=limit)
				{
					flag.clear();
					perm2.clear();
					return;
				}
			}
		}

		perm2.clear();
	}

	flag.clear();
	perm.clear();
}

void CMOEAD::update_reference(CIndividual &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
		}
	}
}

void CMOEAD::update_nadpoint()
{
	int i,j;
	double temp_low  = -1.0e30;

	for (j=0; j < nobj; j++)
	{
		temp_low  = -1e30;
		for (i=0; i < population.size(); i++)
		{
			if (temp_low < population[i].indiv.y_obj[j])
				temp_low = population[i].indiv.y_obj[j];
		}
		nad_point[j] = temp_low;
	}
}


void CMOEAD::mate_selection(vector<int> &list, int cid, int size, int type){
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss   = population[cid].table.size(), id, parent;
    while(list.size()<size)
	{
		if(type==1){
		    id      = int(ss*rnd_uni(&rnd_uni_init));
			parent  = population[cid].table[id];
		}
		else
			parent  = int(population.size()*rnd_uni(&rnd_uni_init));

		// avoid the repeated selection
		bool flag = true;
		for(int i=0; i<list.size(); i++)
		{
			if(list[i]==parent) // parent is in the list
			{
				flag = false;
				break;
			}
		}

		if(flag) list.push_back(parent);
	}
}


int CMOEAD::calculate_crowd_degree(vector<CIndividual> nondonimate_population)
{
	int i,j,k;
	vector<double> crowd_distance;	
	//1.calculate the crowd degree
	int min_index = -1;
	double min_crowd_distance = 1e20, tmp_crowd_degree;
	for (i = 0; i < nondonimate_population.size(); i++)
	{
		crowd_distance=vector<double>(nobj+1, 1.0e20);
		for (j = 0; j < nondonimate_population.size(); j++)
		{	//calculate distance to nobj-nearest neighbour
			double distance_i_j = dist_vector(nondonimate_population[i].y_obj, nondonimate_population[j].y_obj);
			if (distance_i_j < crowd_distance[nobj])
			{
				for (k=nobj; k >=1 && distance_i_j < crowd_distance[k-1]; k--)
				{	
					crowd_distance[k] = crowd_distance[k-1];
				}
				crowd_distance[k] = distance_i_j;//????????????????????????????????????????????
			}
		}
		//adapt multiple as crowd_degree
		tmp_crowd_degree = 1.0;
		for (k = 1; k <= nobj; k++)		tmp_crowd_degree = tmp_crowd_degree * crowd_distance[k];

		if (tmp_crowd_degree < min_crowd_distance)
		{
			min_index = i;
			min_crowd_distance = tmp_crowd_degree;
		}

		crowd_distance.clear();
	}

	return min_index;
}

int CMOEAD::calculate_crowd_degree(vector<CSubproblem> nondonimate_population)
{
	int i,j,k;
	vector<double> crowd_distance;	
	//1.calculate the crowd degree
	int min_index = -1;
	double min_crowd_distance = 1e20, tmp_crowd_degree;
	for (i = 0; i < nondonimate_population.size(); i++)
	{
		crowd_distance=vector<double>(nobj+1, 1.0e20);
		for (j = 0; j < nondonimate_population.size(); j++)
		{	//calculate distance to nobj-nearest neighbor
			double distance_i_j = dist_vector(nondonimate_population[i].indiv.y_obj, nondonimate_population[j].indiv.y_obj);
			if (distance_i_j < crowd_distance[nobj])
			{
				for (k=nobj; k >=1 && distance_i_j < crowd_distance[k-1]; k--)
				{	
					crowd_distance[k] = crowd_distance[k-1];
				}
				crowd_distance[k] = distance_i_j;
			}
		}
		//adapt multiple as crowd_degree
		tmp_crowd_degree = 1.0;
		for (k = 1; k <= nobj; k++)		tmp_crowd_degree = tmp_crowd_degree * crowd_distance[k];

		if (tmp_crowd_degree < min_crowd_distance)
		{
			min_index = i;
			min_crowd_distance = tmp_crowd_degree;
		}

		crowd_distance.clear();
	}

	return min_index;
}

void CMOEAD::external_population_update()
{
	int i,j,k;
	//1. get non dominate solution from the union external_population and new_population
	for (i = 0;  i < new_population.size(); i++)
	{
		int doni_relation = 0;
		for (j = 0; j < external_population.size(); j++)
		{
			doni_relation = donimate_judge(external_population[j].y_obj, new_population[i].y_obj);
			if (doni_relation == -1)//external_population[j]  dominate new_population[i]
				break;
			else if (doni_relation == 1)//new_population[i]  dominate external_population[j]
			{
				external_population.erase(external_population.begin() + j);
				j--;
			}
		}
		if (doni_relation != -1) external_population.push_back(new_population[i]);
	}

	//*************************8.7*********************************
		//**************step1选择范围内的解****************//
			int sum=0;double  xx; double   tmp, tmp_threshold;
			vector< vector<double> > InAreaPoint;
			vector<double>  point1,point;
			vector<double>  pointtoo;
			double min_asf =1e30;   int min_asf_index; double max_threshold = 0.0;  
			//vector<double>  midpoint;//731*********************************
			vector<bool> InArea = vector<bool>(external_population.size(), false);
			vector<double>  inxZ;
			for (int m = 0; m < external_population.size(); m++)
			{    
				if((ZUP-external_population[m].xZ[nvar-1] >= 0)&&(external_population[m].xZ[nvar-1]-ZBP >= 0)  )
				{    InArea[m] = true;     sum++;  }
			}	
			if (sum >= NUM ) 
			{ 
				midpoint.clear();threshold.clear();
				for (int i = 0; i < external_population.size(); i++)
				{   
					if (InArea[i] == true)	
					{
						tmp = fabs(external_population[i].xZ[nvar-1] - ZFL);
						if (tmp < min_asf)      {   min_asf = tmp ; min_asf_index = i;	}
						for(int p=0 ; p<nobj ; p++)  
						{
							point1.push_back(external_population[i].y_obj[p]);
						}
						InAreaPoint.push_back(point1); 
						point1.clear();				
					}
				}
			for(int p = 0 ; p < nobj ; p++)  point.push_back(external_population[min_asf_index].y_obj[p]);
			midpoint = point;    
			point.clear(); 
			    //cout<<"midpoint[0] =  "<<midpoint[0]<<"  midpoint[1] =  "<<midpoint[1]<<endl;	
			for(int p = 0 ; p<nobj ; p++)
			{   
				for (int q = 0; q < sum; q++)
				{
					tmp_threshold = fabs( InAreaPoint[q][p] - midpoint[p]);
					if ( tmp_threshold >= max_threshold )   {  max_threshold  = tmp_threshold ; }
				} 
				xx = relaxfactor1*max_threshold;
				threshold.push_back(xx);
			}
		}
			
			//*************************8.7********************************

	//2. use prefer information to delete the individual in EP
	if (external_population.size() > size_exter_pop)//moea.cpp  size_exter_pop=1.5*MOEAD.pops=30
	{
		vector<bool> InPreferArea = vector<bool>(external_population.size(), false);
		for (i = 0; i < external_population.size(); i++) external_population[i].crowd_rank = 1e30;
		int min_asf_index;
		double tmp_distance;
		double min_asf =1e30;
		//*************************** //
		for (i=0; i<external_population.size(); i++) external_population[i].crowd_rank = 1e30;
	    //*************************** //
		for (j = 0; j < numMultipeLightBeams; j++)
		{	//find the middle individual in the external population
			min_asf_index = -1;
			min_asf = 1.0e30;   double max_threshold = 0.0;   			
			//*********************中心点**************************//
             //**************step1选择范围内的解****************//
			int sum=0; 
			vector< vector<double> > InAreaPoint;			
		//vector<double>  midpoint;//731*********************************
			vector<bool> InArea = vector<bool>(external_population.size(), false);
			for (int m = 0; m < external_population.size(); m++)
			{    
				if((ZUP-external_population[m].xZ[nvar-1] >= 0)&&(external_population[m].xZ[nvar-1]-ZBP >= 0)  )
				{    InArea[m] = true;     sum++;  }
			}	
			if (sum < NUM )  break;
				/*{
				if(!strcmp("ankang20001012", strTestInstance))
				{
				midpoint.push_back(328.0);
				midpoint.push_back(6600.0);
				threshold.push_back(4.0);
				threshold.push_back(2000);
				}
				else if (!strcmp("ankang20030828", strTestInstance))
				{
				midpoint.push_back(325.0);
				midpoint.push_back(3800.0);
				threshold.push_back(4.0);
				threshold.push_back(2000);
				}
				else if(!strcmp("ankang20051001", strTestInstance))
				{
				midpoint.push_back(325.5);
				midpoint.push_back(12000.0);
				threshold.push_back(4.0);
				threshold.push_back(2000.0);
				}
				else if(!strcmp("ankang20100715", strTestInstance))
				{
				midpoint.push_back(326.0);
				midpoint.push_back(5500.0);
				threshold.push_back(4.0);
				threshold.push_back(2000.0);
				}
				else {}
				}*/
			//mark the outrank relationship between the middle individual
			for (i = 0; i < external_population.size(); i++)
			{
				tmp_distance = 0.0;
				for (k = 0; k < nobj; k++)
				{
					//double tmep = (external_population[i].y_obj[k] - external_population[min_asf_index].y_obj[k])/ VetoThreshold[j][k];
					//double tmep = (external_population[i].y_obj[k] - midpoint[k])/ VetoThreshold[j][k];
					double tmep = (external_population[i].y_obj[k] - midpoint[k])/threshold[k];
					tmp_distance +=  tmep * tmep;
				}
				if (tmp_distance <= 1) InPreferArea[i] = true;
				if (tmp_distance <= external_population[i].crowd_rank) 
					external_population[i].crowd_rank = tmp_distance;
			}
			//threshold.clear();//************************812
		}
		vector<int> index;
		for (i = 0; i < external_population.size(); i++)
		{	if (InPreferArea[i] == false)		index.push_back(i);}

		{//Select the solutions whose preference neighbor distances are top size_exter_pop small number 
			//from the preference neighborhoods set  as new external population.
			vector<double> CorwdRank;
			for (i = 0; i < external_population.size(); i++)	CorwdRank.push_back(external_population[i].crowd_rank);
			index.clear();
			for (i = 0; i < external_population.size(); i++)	index.push_back(i);
			minfastsort(CorwdRank, index, index.size(), size_exter_pop);//moea.cpp 160 size_exter_pop =  1.5 * MOEAD.pops;
			sort(index.begin() + size_exter_pop, index.end());
			for (i = index.size() - 1; i >= size_exter_pop; i--) external_population.erase(external_population.begin() + index[i]);
			CorwdRank.clear();
		}
		index.clear();
		InPreferArea.clear();
	}
	//3.diversity preserve to external_population
	if (external_population.size() > size_exter_pop)
	{
		double min_distance = 0.0, dist_i_j = 0.0;
		for (i = external_population.size() - 2; i >= 0; i--)
		{	//get the min distance of individual in external population to the selected individuals
			min_distance = 1e30;
			for (j = external_population.size() - 1; j > i; j--)
			{//calculate the distance of indiv[i] and indiv[j] in external population
				dist_i_j = dist_vector_square(external_population[i].y_obj, external_population[j].y_obj);
				if (dist_i_j < min_distance)	min_distance = dist_i_j;
			}

			if (min_distance < EpsilonDistance * EpsilonDistance)
			{
				external_population.erase(external_population.begin() + i);
				i++;
			}
			
			if (external_population.size() <= size_exter_pop) break;
		}
	}
	//cout<<"size of EP: "<<external_population.size()<<endl;

	//4.diversity preserve to external_population
	int index_min;
	if (external_population.size() > size_exter_pop)
	{	
		while (external_population.size() > size_exter_pop)
		{	//calculate crowd degree
			index_min = calculate_crowd_degree(external_population);			
			//2.delete the individual
			external_population.erase(external_population.begin() + index_min);
		}
	}

	new_population.clear();
}

int CMOEAD::find_sparse_indiv_from_EP()
{
	int i,j,k;
	//3. calculate the crowd degree of each individual in exter_pop if it is added to evol_pop
	vector<double> crowd_distance;
	int index = 0;
	double min_crowd_distance = -1.0,tmp,distance_i_j;	
	//3.1 calculate the crowd degree
	for (j = 0; j < external_population.size(); j++)
	{
		crowd_distance=vector<double>(nobj+1, 1.0e20);
		for (i = 0; i < population.size(); i++)
		{
			//calculate distance to nobj-nearest neighbor
			distance_i_j = dist_vector(external_population[j].y_obj, population[i].indiv.y_obj);
			if (distance_i_j < crowd_distance[nobj])
			{
				for (k=nobj; k >=1 && distance_i_j < crowd_distance[k-1]; k--)
				{	
					crowd_distance[k] = crowd_distance[k-1];
				}
				crowd_distance[k] = distance_i_j;
			}
		}
		//adapt multiple as crowd_degree
		tmp = 1.0;
		for (k = 1; k <= nobj; k++)		tmp *= crowd_distance[k];
		external_population[j].crowd_rank = tmp;
		crowd_distance.clear();
		if (external_population[j].rank < external_population[index].rank)
		{
			index = j;//一个在一个不在，选在范围内的那个
		}
		else if ( (external_population[j].rank == 0) && (external_population[index].rank == 0) )
		{
			if ( external_population[j].crowd_rank > external_population[index].crowd_rank ) 
			index = j;
			//都在范围内，选最稀疏的
		}
		else if ( (external_population[j].rank == 1) && (external_population[index].rank == 1) )
		{
			if ( external_population[j].crowd_rank < external_population[index].crowd_rank )
			index = j;
			//都不在时选密集处的那个
		}
		//else if ( (external_population[j].rank == 2) && (external_population[index].rank == 2) )
		//{
		//	if ( external_population[j].crowd_rank > external_population[index].crowd_rank ) index = j;
		//}
	}
	return index;
}

void CMOEAD::delete_undesired_subproblems()
{
	int i,j,k;
	// 1. prior avoid visiting the  undesired areas; 
	vector<bool> InPreferArea = vector<bool>(population.size(), false);  //record the solution whether in prefer area or not 
	for (i = 0; i < population.size(); i++) population[i].indiv.crowd_rank = 1e30; //record the distance to the middle individual
	double  tmp_distance;
	//*****************
	for (i = 0; i < population.size(); i++)  population[i].indiv.crowd_rank = 1e30;
	//*****************

	for (j = 0; j < numMultipeLightBeams; j++)
	{
		
//*****************************************选择范围内的解************************************************//
		int sum=0;  
		vector< vector<double> > InAreaPoint;
		//vector<double>  midpoint;//731*********************************

		vector<bool> InArea = vector<bool>(external_population.size(), false);
		for (int m = 0; m < external_population.size(); m++)
		{    
			if(  (ZUP - external_population[m].xZ[nvar-1] >= 0)&&(external_population[m].xZ[nvar-1]-ZBP >= 0)  )
			{    InArea[m] = true;     sum++;  }
		}
		cout<<"__________waibu______________total="<<sum<<endl;	
		if (sum < NUM )  break;

		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		vector<bool> InA = vector<bool>(population.size(), false);
		int num = 0;
		for (int m = 0; m < population.size(); m++)
		{    
			if(  (ZUP - population[m].indiv.xZ[nvar-1] >= 0)&&(population[m].indiv.xZ[nvar-1]-ZBP >= 0)  )
			{    InA[m] = true;     num++;  }
		}	
		cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXtotal="<<num<<endl;
		for (int i = 0; i < population.size(); i++)
		{   
		 if (InA[i] == true)	
			{
						cout<<i<<"  mm  ";
			}
		}
		//**********&&&&&&&&&&&&&&&&&&&&&&**********812&&&&&&&&&&&&&&&&&&&

		cout<<endl;
		for(int s = 0 ; s < nobj ; s++)
			cout<<"midpoint"<<" = "<<midpoint[s]<<"          ";
		cout<<"threshold[0] = "<<threshold[0]<<"   threshold[1] = "<<threshold[1]<<endl;
		
		//midpoint.clear();//731*********************************
//******************************************中心点******************************************************//
		
		for (i = 0; i < population.size(); i++)
		{
			tmp_distance = 0.0;
			for (k = 0; k < nobj; k++)
			{
				//double tmep = (population[i].indiv.y_obj[k] - population[min_asf_index].indiv.y_obj[k])/ VetoThreshold[j][k];
				//double tmep = (population[i].indiv.y_obj[k] - midpoint[k])/ VetoThreshold[j][k];
				double tmep = (population[i].indiv.y_obj[k] - midpoint[k])/ (threshold[k]*relaxfactor);
				tmp_distance +=  tmep * tmep;
			}
			tmp_distance = sqrt(tmp_distance);
			if (tmp_distance <= 1) InPreferArea[i] = true;
			if (tmp_distance < population[i].indiv.crowd_rank)
				population[i].indiv.crowd_rank = tmp_distance;
		}
		//threshold.clear();//*****************************812
	}
	//save the subproblems with weight vectors (1,0,...,0),(0,1,0,...,0),...,(0,...,0,1)
	//for (i = 0; i < nobj; i++) InPreferArea[i] = true;
	vector<int> index;
	vector<double> dist;//???????????????????????????????????????????????????
	for (i = 0; i < population.size(); i++)
	{
		if (InPreferArea[i] == false)
		{
			index.push_back(i);
			dist.push_back(population[i].indiv.crowd_rank);
		}
	}
	
	//for (i = 0; i < index.size(); i++) cout<<population[index[i]].indiv.crowd_rank<<" ";
	//cout<<endl;

	if (index.size() >= num_update_weight)
	{	//prior avoid visiting the  undesired areas;
		minfastsort(dist, index, index.size(), index.size() - num_update_weight);
		sort(index.begin()+(index.size() - num_update_weight), index.end());
		int num = (index.size() - num_update_weight);
		//for (i = index.size() - 1; i >= (index.size() - num_update_weight); i--) cout<<population[index[i]].indiv.crowd_rank<<" ";
		//cout<<endl;
		for (i = index.size() - 1; i >= num; i--)		population.erase(population.begin() + index[i]);
	}
	else 
	{	//prior avoid visiting the  undesired areas;
		//for (i = index.size() - 1; i >= 0; i--) cout<<population[index[i]].indiv.crowd_rank<<" ";
		sort(index.begin(), index.end());
		for (i = index.size() - 1; i >= 0; i--) population.erase(population.begin() + index[i]);
		//then delete the subproblem whose optimal solutions are in the prefer area but have high crowd;
		int num = num_update_weight - index.size();
		for (i = 0; i < num; i++)
		{	//find the subproblem with highest crowd degree
			int min_index = calculate_crowd_degree(population);	
			//cout<<population[min_index].indiv.crowd_rank<<" ";
			population.erase(population.begin() + min_index);			
		}
		//cout<<endl;
	}
	InPreferArea.clear();
	index.clear();
	dist.clear();
}

void CMOEAD::add_spare_subproblems_to_prefer_area()
{
	int i,j,k,m;
	//1 calculate outrank relation
	//*********************
	for (i = 0 ; i < external_population.size(); i++) external_population[i].crowd_rank = 1e30;
	//************************
	for (i = 0 ; i < external_population.size(); i++)
	{
		int min_asf_index;
		double min_asf, tmp_distance;
		for (j = 0; j < numMultipeLightBeams; j++)
		{	//find the middle individual in the external population
			min_asf_index = -1;
			min_asf = 1.0e30;   double max_threshold = 0.0;   
//******************************************选择范围内的解***********************************************//
			int sum=0; 
			vector< vector<double> > InAreaPoint;		
			vector<bool> InArea = vector<bool>(external_population.size(), false);
			vector<double>  inxZ;
			for (int m = 0; m < external_population.size(); m++)
			{    
				if(  (ZUP - external_population[m].xZ[nvar-1] >= 0)&&(external_population[m].xZ[nvar-1]-ZBP >= 0))
				{    InArea[m] = true;     sum++;  }
			}	
		
			if(sum < NUM ) break;
//***************************************中心点***********************************************************//
			//mark the outrank relationship between the middle individual
			for (int m = 0; m < external_population.size(); m++)
			{
				tmp_distance = 0.0;
				for (k = 0; k < nobj; k++)
				{
					//double tmep = (external_population[m].y_obj[k] - external_population[min_asf_index].y_obj[k])/ VetoThreshold[j][k];
					//double tmep = (external_population[m].y_obj[k] - midpoint[k])/ VetoThreshold[j][k];
					double tmep = (external_population[m].y_obj[k] - midpoint[k])/(threshold[k]*relaxfactor2);
					tmp_distance  +=  tmep * tmep;
				}
				if (tmp_distance <= external_population[m].crowd_rank) 
					external_population[m].crowd_rank = tmp_distance;
			}
			//threshold.clear();//*********************************812
		}
		if (external_population[i].crowd_rank <= 1) external_population[i].rank = 0;
		else external_population[i].rank = 1;
	}
	//2. reprocess---  not consider the boundary remote solution in exter_pop  ----//
	//for (i = 0 ; i < external_population.size(); i++)
	//{
	//	for (j = 0; j < nobj; j++) 
	//	{
	//		if ( fabs(external_population[i].y_obj[j] - idealpoint[j]) < 1e-6)
	//			
	//	}external_population[i].rank = 2; 
	//}	
	vector<double> direction=vector<double>(nobj, 0.0);
	//use weight solution (WS)- transformation to calculate the optimal subproblem to the most sparse solution in EP
	for (i = 0; i < num_update_weight; i++)
	{	//3.1 calculate the min crowd_degree of individual in external_population to evol_pop
		int index_min_crowd_degree = find_sparse_indiv_from_EP();		
		//3.2 add new subproblem
		CSubproblem sub;
		int rnd_index = int(rnd_uni(&rnd_uni_init)*population.size());
		sub.saved = population[rnd_index].indiv;//let init delta max 
		sub.indiv = external_population[index_min_crowd_degree];//	
		sub.utility = 1.0;
		//the direction of landa_t is same as the f-z* (external_population[index_min_crowd_degree].y_obj - idealpoint)
		for (m = 0; m < nobj; m++) direction[m] = external_population[index_min_crowd_degree].y_obj[m] - idealpoint[m];
		for (m = 0; m < nobj; m++)
		{	if (direction[m] < 1e-20)	direction[m] = 1e-7;}
		double sum = 0.0;
		for (m = 0; m < nobj; m++)	sum += 1.0 / direction[m];
		for (m = 0; m < nobj; m++)	sub.namda[m] = (1.0 / direction[m]) / sum;	

		population.push_back(sub);		
		external_population[index_min_crowd_degree].rank = 3;
	}
	direction.clear();
	update_neighbour_table();
}

void CMOEAD::organize_merge_split()//组织 合并 分离
{
	external_population_update();
	cout<<midpoint[0]<<"      "<<midpoint[1]<<endl;
	delete_undesired_subproblems();
	add_spare_subproblems_to_prefer_area();	
}

void CMOEAD::evol_population()
{
	// random order of subproblems at each generation
	//vector<int> order(vector<int>(pops,0));
	//for(int i=0; i<pops; i++)  order[i] = i;
	//random_shuffle(order.begin(), order.end());shuffle：洗牌
	vector<int> order;	this->tour_selection(10, order);

    for(int sub=0; sub<order.size(); sub++)
	{
		int c_sub = order[sub];    // random order
		int type;
        double rnd = rnd_uni(&rnd_uni_init);
		// mating selection based on probability
		if(rnd<prob)     type = 1;   // from neighborhood
		else             type = 2;   // from population
		// select the indexes of mating parents
		vector<int> plist;
		mate_selection(plist, c_sub, 2, type);  // neighborhood selection
		// produce a child solution
		CIndividual child;
		//double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
		//diff_evo_xoverB(population[c_sub].indiv,population[plist[0]].indiv,population[plist[1]].indiv, child, rate2);
		if (!strcmp("SBX",strCrossType))
		{
			real_sbx_xoverB(population[plist[0]].indiv,population[plist[1]].indiv,child);
		}
		else
		{
			double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
			diff_evo_xoverB(population[c_sub].indiv,population[plist[0]].indiv,population[plist[1]].indiv, child, rate2);
		}		
		plist.clear();	
		// apply polynomial mutation
		realmutation(child, 1.0/nvar);
		// evaluate the child solution
		child.obj_eval();
		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		update_problem(child, c_sub, type);
		nfes++;//************************************评价次数加一&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
		//if (nfes >= start_time_exter_pop_update)	
			new_population.push_back(child);
	}

	if (!strcmp("_TCH2",strFunctionType)) 
		update_nadpoint();
	order.clear();
}


void CMOEAD::exec_emo(int run)
{
	int i,j; int frequence_update_weight;
    char filename[1024];
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;
	// initialization

	nfes = 0;
	prob = 1;//?

	rate_update_weight = 0.2;//
	evol_rate = 0.4;//

	vector <double> hv;
	vector <double> total;
	start_time_exter_pop_update = num_max_evulation*evol_rate*0.8;//20,0000*0.4*0.8=6,4000代开始调整权向量
	num_update_weight = int(pops*rate_update_weight+1e-20);//20*0.1
	
	init_population();
    //init_neighbourhood();
	init_prefer_parameter(); 
	update_neighbour_table();

	int gen = 5;//begin  at  5
	//max_gen = int(5.0*num_max_evulation/pops);
	//int frequence_update_weight = int( (1.0 - evol_rate) * max_gen * rate_update_weight / 5);//(1-0.4)*5,0000*0.1/2=1500调整权向量的代数间隔
	
	if (max_gen < 100) max_gen =100;

	num_max_evulation = int(max_gen*pops/5);//最大评价次数=最大进化代数*每代数目即25000*20/5=10,0000
	
	 frequence_update_weight = 5000;

	int cur_gen = 0; //??????????????????????????????????????????????????????????????????????????????????

	bool mark1 = false; 
	while(nfes<num_max_evulation)
	{
		evol_population();	
		gen++;
		if(gen%30==0)    comp_utility();		
		external_population_update();
		threshold.clear();//************************812
		midpoint.clear();//731*********************************
		//if ( external_population.size() >= population.size() )  mark1 = true;
		//if (nfes >=  num_max_evulation*evol_rate )//评价次数大与20,0000*0.4=8,0000
		//{
			//if (nfes%frequence_update_weight==0 && mark1)//每隔1500代而且|ExP|>N
			if (nfes%frequence_update_weight==0 )//每隔5000次函数评价
			{
				organize_merge_split();
				midpoint.clear();//731********************************* 
				threshold.clear();//************************812
			}
			
		//}
		/*int total = 0;
		for (int m = 0; m < population.size(); m++)
		{    
		if(  (ZUP - population[m].indiv.xZ[nvar-1] > 0)&&(population[m].indiv.xZ[nvar-1]-ZBP > 0)  )
		{       total++;  }
		}		
		if(total > 19) break;*/

		//if (nfes >=  num_max_evulation*evol_rate && gen%frequence_update_weight==0 && mark1)
		//if (nfes >=  num_max_evulation*evol_rate && nfes%frequence_update_weight==0)

		//if (nfes%frequence_update_weight==0)
		//{
		//	sprintf(filename,"POFmid/POF_MOEAD_%s_OBJ%d_RUN%d_%d.txt",strTestInstance, nobj, run, nfes);
		//	save_front(population,filename);
		//	sprintf(filename,"POSmid/POS_MOEAD_%s_OBJ%d_RUN%d_%d.txt",strTestInstance, nobj, run, nfes);
		//	save_pos(population,filename);
		//	sprintf(filename,"UpstreamWaterLevelmid/UWL_MOEAD_%s_OBJ%d_RUN%d_%d.txt",strTestInstance, nobj, run,nfes);//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd
		//	save_Z(population,filename);//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd
		//	
		//}
		//if ((nfes < frequence_update_weight && nfes%(frequence_update_weight/5)==0)||nfes%frequence_update_weight==0)  
		if ( nfes%(frequence_update_weight/5)==0)
		
			hv.push_back(this->Calculate_HV_for_RFC());
		if (nfes%frequence_update_weight==0)   total.push_back(this->Calculate_TOTAL_for_RFC());
	}
		
	
	// save the final population - X space
	sprintf(filename,"POS/POS_MOEAD_%s_OBJ%d_RUN%d.txt",strTestInstance, nobj, run);
	save_pos(population, filename);
	sprintf(filename,"POF/POF_MOEAD_%s_OBJ%d_RUN%d.txt",strTestInstance, nobj, run);
	save_front(population, filename);
	sprintf(filename,"UpstreamWaterLevel/UWL_MOEAD_%s_OBJ%d_RUN%d.txt",strTestInstance, nobj, run);
	save_Z(population,filename);//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd

	HV.push_back(hv);
	hv.clear();	
	//HV.push_back(this->Calculate_HV_for_RFC());//addaddaddaddaddaddaddaddaddaddaddaddadd
	TOTAL.push_back(total);
	total.clear();

	external_population_update();
	midpoint.clear();//731*********************************
	threshold.clear();//*************************812
	organize_merge_split();
	midpoint.clear();//731*********************************
	threshold.clear();//**********************812
	new_population.clear();

	////addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd
	//vector <CSubproblem>       Newpopulation;
	//for(int n=0; n<pops; n++)
	//{
	//	CSubproblem sub;
	//	double a = 1.0*n/(pops - 1);
	//	sub.namda.push_back(a);
	//	sub.namda.push_back(1-a);

	//	Newpopulation.push_back(sub);
	//}

	//for(i=0; i<Newpopulation.size(); i++)
	//{
	//	for (j=0; j < nvar; j++) Newpopulation[i].indiv.x_var[j] = population[i].indiv.y_obj[1];//y_obj[1]=max(x_var)最大下泄流量cec09中ankang
	//	Newpopulation[i].indiv.obj_eval();
	//}
	////save data
	//sprintf(filename,"PPOF/POF_MOEAD_%s_OBJ%d_RUN%d.txt", strTestInstance, nobj, run);
	//save_front(Newpopulation,filename);
	//sprintf(filename,"PPOS/POS_MOEAD_%s_OBJ%d_RUN%d.txt", strTestInstance, nobj, run);
	//save_pos(Newpopulation,filename);
	//sprintf(filename,"PUpstreamWaterLevel/UWL_MOEAD_%s_OBJ%d_RUN%d.txt",strTestInstance, nobj, run);
	//save_Z(Newpopulation,filename);//addaddaddaddaddaddaddaddaddaddaddadda
	//Newpopulation.clear();//addaddadd
	population.clear();
	idealpoint.clear();
	new_population.clear();
	external_population.clear();
}

void CMOEAD::load_parameter()
{
	char filename[1024];
	sprintf(filename,"ParameterSetting/\%s_%d.txt", strTestInstance, nobj);

	char temp[1024];
	std::ifstream readf(filename);
	readf.getline(temp, 1024);
	//puts(temp);

	readf>>pops;
	readf>>max_gen;
	readf>>niche;
	readf>>limit;
	readf>>prob;
	readf>>rate;

	readf.close();

	int i,j;
	for (i = 0; i < AspirationPoint.size(); i++)		AspirationPoint[i].clear();
	AspirationPoint.clear();
	for (i = 0; i < ReservationPoint.size(); i++)		ReservationPoint[i].clear();
	ReservationPoint.clear();
	for (i = 0; i < VetoThreshold.size(); i++)		VetoThreshold[i].clear();
	VetoThreshold.clear();
	for (i = 0; i < Weight.size(); i++)		Weight[i].clear();
	Weight.clear();

	//read the preference information in lbs-NSGA-II from file
	sprintf(filename,"PreferParameterSetting/OnePreferArea/\%s_%d.txt", strTestInstance, nobj);
	std::ifstream readf2(filename);
	readf2.getline(temp, 1024);//filter the explain
	
	readf2>>numMultipeLightBeams;
	//allocate the Storage Resources 
	for (i = 0; i < numMultipeLightBeams; i++)
	{
		AspirationPoint.push_back(vector <double> (0.0));
		for (j=0; j<nobj; j++) AspirationPoint[i].push_back(0.0);
		ReservationPoint.push_back(vector<double> (0.0));
		for (j=0; j<nobj; j++) ReservationPoint[i].push_back(0.0);
		VetoThreshold.push_back(vector <double> (0.0));
		for (j=0; j<nobj; j++) VetoThreshold[i].push_back(0.0);
		Weight.push_back(vector<double> (0.0));
		for (j=0; j<nobj; j++) Weight[i].push_back(0.0);
	}

	for (i = 0; i < numMultipeLightBeams; i++)
	{
		for (j = 0; j < nobj; j++)
			readf2>>AspirationPoint[i][j];
		for (j = 0; j < nobj; j++)
			readf2>>ReservationPoint[i][j];
		for (j = 0; j < nobj; j++)
			readf2>>VetoThreshold[i][j];
	}
	readf2>>EpsilonDistance;

	readf2.close();
	//calculate the preference direction
	vector<double> tmp = vector<double>(nobj, 0.0);
	double total = 0.0;

	for (i = 0; i < numMultipeLightBeams; i++)
	{
		for (j = 0; j < nobj; j++)
		{
			tmp[j] = 1.0 / (ReservationPoint[i][j] - AspirationPoint[i][j]+1e-6);
			total +=  tmp[j];
		}

		for (j = 0; j < nobj; j++)
			Weight[i][j] = tmp[j] / total;
	}
}


void CMOEAD::save_front(vector <CSubproblem> pop,char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);

	for(int n=0; n<pop.size(); n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<pop[n].indiv.y_obj[k]<<"  ";		      
		   fout<<"\n";	
	}

	fout.close();
}

void CMOEAD::save_exter_pop(char savefilename[1024])
{
	std::fstream fout;
	fout.open(savefilename,std::ios::out);

	for(int n=0; n<external_population.size(); n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<external_population[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEAD::save_pos(vector <CSubproblem> pop,char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pop.size(); n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<pop[n].indiv.x_var[k]<<"  ";
		fout<<"\n";
	}
	
	fout.close();
}

void CMOEAD::save_individual(char saveFilename_x[1024],char saveFilename_y[1024],vector<double> x,vector<double> y)
{
	std::fstream fout;
	fout.open(saveFilename_x,std::ios::out);
	int n;
	for(n=0; n<nvar; n++)
	{
		fout<<x[n]<<"  ";
	}
	fout.close();

	//std::fstream fout;
	fout.open(saveFilename_y,std::ios::out);
	for(n=0; n<nobj; n++)
	{
		fout<<y[n]<<"  ";
	}
	fout.close();
}
bool CMOEAD::boundy_checkout(CIndividual child)
{
	int i;
	for (i=0; i<nvar; i++)
	{
		if (child.x_var[i]<lowBound[i]||child.x_var[i]>uppBound[i]) return true;
	}
	for (i=0; i<nobj; i++)
	{
		if (child.y_obj[i]<0) return true;
	}
	return false;
}

void CMOEAD::save_HV(vector<vector<double>> hv_value, char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int m = 0;m < hv_value.size();m++)
	{
	    for(int n = 0; n< hv_value[0].size(); n++)
		fout<<hv_value[m][n]<<"  ";	
	    fout<<endl;  
	}

	fout.close();
}
void CMOEAD::save_TOTAL(vector<vector<double>> total_value, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int m = 0;m < total_value.size();m++)
	{
		for(int n = 0; n< total_value[0].size(); n++)
			fout<<total_value[m][n]<<"  ";	
		fout<<endl;  
	}

	fout.close();
}

double CMOEAD::Calculate_HV_for_RFC()
{
	int i,j;
	double Thershold1,Thershold2;
	vector<double> ReferencePoint=vector<double>(nobj,0.0);
	if(!strcmp("ankang20001012", strTestInstance))
	{
		/*Thershold1 = 320;
		Thershold2 = 330;*/
		Thershold1 = ZBP;
		Thershold2 = ZUP;
		ReferencePoint[0] = 328.0;//ReferencePoint[0] = 330.0;//old: ReferencePoint[0] = 326.0;
		ReferencePoint[1] = 11000.0;//ReferencePoint[1] = 11000.0;//old: ReferencePoint[1] = 6000.0;
	}
	else if (!strcmp("ankang20030828", strTestInstance))
	{
		/*Thershold1 = 320;
		Thershold2 = 330;*/
		Thershold1 = ZBP;
		Thershold2 = ZUP;
		ReferencePoint[0] = 328.0;//ReferencePoint[0] = 330.0;
		ReferencePoint[1] = 6000.0;//ReferencePoint[1] = 7826.21;
	}
	else if(!strcmp("ankang20051001", strTestInstance))
	{
		/*Thershold1 = 320;
		Thershold2 = 330;*/
		Thershold1 = ZBP;
		Thershold2 = ZUP;
		ReferencePoint[0] = 328.0;//ReferencePoint[0] = 330.0;
		ReferencePoint[1] = 14000.0;//ReferencePoint[1] = 17000.0;
	}
	else if(!strcmp("ankang20100715", strTestInstance))
	{
		/*Thershold1 = 320;
		Thershold2 = 330;*/
		Thershold1 = ZBP;
		Thershold2 = ZUP;
		ReferencePoint[0] = 328.0;//ReferencePoint[0] = 330.0;
		ReferencePoint[1] = 7000.0;//ReferencePoint[1] = 7784.36;
	}
	else
	{}
	//get nondominate solution from the evoluationary population
	vector <CIndividual> TempData;
	for (i = 0;  i < population.size(); i++)
	{
		int doni_relation = 0;
		for (j = 0; j < TempData.size(); j++)
		{
			doni_relation = donimate_judge(TempData[j].y_obj, population[i].indiv.y_obj);
			if (doni_relation == -1)//TempData[j].y_obj  dominate population[i].indiv.y_obj
				break;
			else if (doni_relation == 1)//population[i].indiv.y_obj  dominate TempData[j].y_obj
			{
				TempData.erase(TempData.begin() + j);
				j--;//???????????????   ?????????????????????????????????  ??????????????????????????????
			}
		}
		if (doni_relation != -1) TempData.push_back(population[i].indiv);
	}
	//remove the nonfeasible solution
	for (i = 0; i < TempData.size(); i++)
	{
		if (TempData[i].xZ[nvar-1] < Thershold1 || TempData[i].xZ[nvar-1] > Thershold2 || TempData[i].y_obj[0] >= ReferencePoint[0]||TempData[i].y_obj[1] >= ReferencePoint[1])//prefer solution
		{
			TempData.erase(TempData.begin() + i);
			i--;//移除当前位置，当前位置已是之前的下一个位置，还没有考察，--之后运行循环的++，i不变即可考察
		}
	}
	// sort the nondominating solution according to the f1.
	CIndividual temp;
	for (i = 0; i < TempData.size() - 1; i++)
	{
		for (j = i+1; j < TempData.size(); j++)
		{
			if (TempData[j].y_obj[0]<TempData[i].y_obj[0])
			{//swap(TempData[i],TempData[j]);
				temp = TempData[i];
				TempData[i] = TempData[j];
				TempData[j] = temp;
			}
		}
	}
	//calculate the HV
	double hv1 = 0.0;
	if (TempData.size() == 1) hv1 = (ReferencePoint[1]-TempData[0].y_obj[1])*(ReferencePoint[0]-TempData[0].y_obj[0]);
	else if (TempData.size() > 1)
	{
		hv1 = (ReferencePoint[1]-TempData[0].y_obj[1])*(ReferencePoint[0]-TempData[0].y_obj[0]);
		for (i = 1; i < TempData.size(); i++)
		{
			hv1 = hv1 + (TempData[i-1].y_obj[1]-TempData[i].y_obj[1]) * (ReferencePoint[0]-TempData[i].y_obj[0]);
		}
	}
	ReferencePoint.clear();
	TempData.clear();
	return hv1;
}

int CMOEAD::Calculate_TOTAL_for_RFC()
{  
	int total = 0;
	for (int m = 0; m < population.size(); m++)
	{    
		if(  (ZUP - population[m].indiv.xZ[nvar-1] >= 0)&&(population[m].indiv.xZ[nvar-1]-ZBP >= 0)  )
		{       total++;  }
	}	
	return total;
}
void CMOEAD::save_Z(vector <CSubproblem> pop,char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pop.size(); n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<pop[n].indiv.xZ[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

#endif
