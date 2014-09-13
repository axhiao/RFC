/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/


#include "algorithm.h"


int main()
{
	int i,j;
	int pop;
	std::ifstream readf("TestInstance.txt");


	//num_max_evulation = 300000;
    // "_TCH1": Chebyshev, "_TCH2": normalized Tchebycheff, "_PBI": Penalty-based BI, 
	// "_WS": Weight Sum
    //strcpy(strFunctionType,"_TCH1"); 
	strcpy(strFunctionType,"_TCH1"); 
	
//	int numOfInstance;
//	readf>>numOfInstance;

//	printf("\n -- %d Instances are being tested ---\n\n", numOfInstance);
	//set needed test case index from begin to end 
	int fun_start_num = 52;
	int fun_end_num = 55;
	int num_run = 40;
	//skip the front test instance
	int seq;
	for (i =1;i<fun_start_num;i++)
	{
		readf>>seq;
		readf>>strTestInstance;
		readf>>nvar;
		readf>>nobj;

	}
	//
	char outFilename[1024];
	sprintf(outFilename, "LOG/LOG_MOEAD_IGD.dat", strTestInstance);
	std::fstream fout2;
	fout2.open(outFilename,std::ios::out);
	
	for(inst=fun_start_num; inst<=fun_end_num; inst++)
	{
		// the parameter setting of test instance
		readf>>seq;
		readf>>strTestInstance;
		readf>>nvar;
		readf>>nobj;

		preprocess_for_reservior(); //ÔØÈëÊý¾Ý

		position_parameters = position_parameters * (nobj -1);
		distance_parameters = distance_parameters * 2;

		printf("\n -- Instance: %s, %d variables %d objectives \n\n", strTestInstance, nvar, nobj);

		xboundy();//get  the domain of x_var.
		/*
		if (inst >= 11 && inst <= 23||inst >= 34 && inst <= 40)
		{
			if (nobj == 2)
			{
				num_max_evulation = 30000;
				strcpy(strCrossType,"SBX");
				//pops = 100;
			}
			else if (nobj == 3)
			{
				num_max_evulation = 60000;
				strcpy(strCrossType,"SBX");
				//pops = 100;
			}
		}
		else if (inst >= 1 && inst <= 10 || inst >= 24 && inst <= 33|| inst >= 41 && inst <= 133)
		{
			if (nobj == 2)
			{
				num_max_evulation = 300000;
				strcpy(strCrossType,"DE");
				//pops = 600;
			}
			else if (nobj == 3)
			{
				num_max_evulation = 300000;
				strcpy(strCrossType,"DE");
				//pops = 1000;
			}
		}
		else if (inst >= 134 && inst <= 140)
		{
			if (nobj == 2)
			{
				num_max_evulation = 20000;
				strcpy(strCrossType,"SBX");
				//pops = 100;
			}
			else if (nobj == 3)
			{
				num_max_evulation = 30000;
				strcpy(strCrossType,"SBX");
			}
			else if (nobj == 4)
			{
				num_max_evulation = 75000;
				strcpy(strCrossType,"SBX");
				//pops = 100;
			}
			else if (nobj == 6)
			{
				num_max_evulation = 100000;
				strcpy(strCrossType,"SBX");
			}
		}
		else if (inst >= 141 && inst <= 145)
		{
			num_max_evulation = 70000;
			strcpy(strCrossType,"SBX");
		}
*/

		clock_t start, temp, finish;
		double  duration, last = 0;
		start = clock();

		std::fstream fout;
		char logFilename[1024];
		sprintf(logFilename, "LOG/LOG_MOEAD_%s.dat", strTestInstance);
		fout.open(logFilename,std::ios::out);

		fout<<"Time: "<<endl;
		fout<<"Inst: "<<strTestInstance<<endl;

		CMOEAD MOEAD;
		MOEAD.load_parameter();
		pop = MOEAD.pops;

		size_exter_pop = 4.0 * MOEAD.pops;
		
		for(int run=1; run<=num_run; run++)
		{
			printf("\n -- %d-th run  -- \n", run);
			
			MOEAD.exec_emo(run);

			temp = clock();
			duration = (double)(temp - start) / CLOCKS_PER_SEC;
			fout<<run<<":"<<duration - last<<" ";
			last = duration;
			if(run%10==0) fout<<"\n";
		}

		fout<<"\n\n";

		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;

		fout<<"Mean  CPU Time  "<<duration/num_run<<" seconds"<<endl;
		fout<<"Total CPU Time  "<<duration<<" seconds"<<endl;
		fout.close();

		char savefilename[1024];//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadd
		sprintf(savefilename,"HV/HV_MOEAD_%s_OBJ%d.txt",strTestInstance, nobj);//addaddaddaddaddadd
		MOEAD.save_HV(HV, savefilename);//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddadda
		HV.clear();//addaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddaddad
		char saveFilename[1024];//addaddaddadd830addaddaddaddaddaddaddaddaddaddadd
		sprintf(saveFilename,"TOTAL/TOTAL_MOEAD_%s_OBJ%d.txt",strTestInstance, nobj);//addad830ddaddaddd
		MOEAD.save_TOTAL(TOTAL, saveFilename);//addaddaddaddaddad830daddaddaddaddaddaddaddaddaddaddadd
		TOTAL.clear();//addaddaddaddaddaddaddaddadda830ddaddaddaddaddaddaddaddaddaddadaddddadd
/*

		//calculate the merit of IGD
		double * * pf=new double *[pop]; 
		for(i=0;i<pop;i++) 	pf[i]=new double[nobj]; 

		int pop_test=100;
		if (nobj == 2)	pop_test = 100;
		else if (nobj == 3) pop_test = 150;
		double * * pf_find=new double *[pop_test]; 
		for(i=0;i<pop_test;i++) 	pf_find[i]=new double[nobj]; 

		double *temp_pf = new double [nobj];
		//load the pf_find and ideal PF	
		double IGD = 0;
		double IGD_TOTAL = 0;
		int PFtrue_number = 0;

		std::fstream fout_igd;
		char igd_Filename[1024];
		sprintf(igd_Filename, "IGD/IGD_MOEAD_%s.txt", strTestInstance);
		fout_igd.open(igd_Filename,std::ios::out);
		
		for (run=1; run<=num_run; run++)
		{
			IGD = 0;
			//load the  pf_find
			sprintf(logFilename,"POF/POF_MOEAD_%s_RUN%d.txt",strTestInstance,run);
			std::ifstream readf(logFilename);			
			for (int i=0; i<pop; i++)
			{
				for (int j=0; j<nobj; j++)	readf>>pf[i][j];
			}
			readf.close();

			//get nondominate solution in pf
			int index_pop = pop;
			for (i=0; i<index_pop; i++)
			{
				for (int j=i+1; j<index_pop; j++)
				{
					double tt = donimate_judge(pf[i],pf[j]);
					if (tt>0)
					{
						for (int k=0; k<nobj; k++) 
						{
							double temp = pf[index_pop-1][k];
							pf[index_pop-1][k] = pf[i][k];
							pf[i][k] = temp;
						}
						i--;
						index_pop--;
						continue;
					}
					else if (tt<0)
					{
						for (int k=0; k<nobj; k++) 
						{
							double temp = pf[index_pop-1][k];
							pf[index_pop-1][k] = pf[j][k];
							pf[j][k] = temp;
						}
						j--;
						index_pop--;
					}					
				}			
			}
			char saveFilename[100];
			sprintf(saveFilename,"POF_NONDONIMATE/POF_MOEAD_%s_RUN%d.txt",strTestInstance,run);
			std::fstream fout4;
			fout4.open(saveFilename,std::ios::out);
			for(int n=0; n<index_pop; n++)
			{
				for(int k=0;k<nobj;k++)
					fout4<<pf[n][k]<<"  ";
				fout4<<"\n";
			}
			fout4.close();
			//random select start point,and run 10 times
			vector<bool> flag=vector<bool>(index_pop,true);
			if (nobj == 2||nobj == 3)
			{
				double igd_random = 0;
				for (int random_count=0; random_count<num_of_random; random_count++)
				{
					igd_random = 0;

					int count_index = 0;
					flag=vector<bool>(index_pop,true);
					//random select start point
					int rand_start_index = int(index_pop*rnd_uni(&rnd_uni_init));
					flag[rand_start_index] = false;
					for (int j=0; j<nobj; j++) pf_find[count_index][j] = pf[rand_start_index][j];
					count_index++;
					int point_index;
					for (int num_point=2; num_point<=100; num_point++)
					{
						point_index = find_point(pf,flag,pf_find,count_index,index_pop);
						for (int j=0; j<nobj; j++) pf_find[count_index][j] = pf[point_index][j];
						count_index++;
						flag[point_index] = false;
					}
					//load the ideal PF
					sprintf(logFilename,"ParameterSetting/PFture/%s.dat",strTestInstance);
					std::ifstream readf2(logFilename);
					
					//if (nobj == 2)
					//{
					//	PFtrue_number = 1000;		
					//}
					//else if (nobj == 3)
					//{
					//	PFtrue_number = 10000;
					//}
					double *temp=new double[nobj];
					double min_dist= 0,temp_ela_dist = 0,total_dist = 0;

					PFtrue_number = 0;
					while (!readf2.eof())
					{
						for (int n=0; n<nobj; n++)	readf2>>temp[n];
						PFtrue_number++;

						min_dist = 1e8;
						for (int j=0; j<pop_test; j++)
						{
							temp_ela_dist = 0;
							for (int n=0; n<nobj; n++)	temp_ela_dist += (temp[n]-pf_find[j][n])*(temp[n]-pf_find[j][n]);
							temp_ela_dist = sqrt(temp_ela_dist);
							if (temp_ela_dist < min_dist)	min_dist = temp_ela_dist;
						}
						//calculate the IGD
						total_dist +=  min_dist;
					}
					delete temp;
					readf2.close();
					//calculate the IGD
					igd_random = total_dist/(PFtrue_number-1);		
					IGD += igd_random/num_of_random;			
				}								
			}
			flag.clear();
			fout_igd<<IGD<<endl;
			IGD_TOTAL += IGD;			
		}
		printf("\n -- Instance: %s, %d variables %d objectives the average of IGD is %f\n\n", strTestInstance, nvar, nobj,IGD_TOTAL/num_run);
		fout_igd<<endl<<"Instance: "<<strTestInstance<<", "<<nvar<<" variables "<<nobj<<" objectives the average of IGD is "<<IGD_TOTAL/num_run<<endl;

		for(i=0;i<pop_test;i++) 	delete pf_find[i];
		delete pf_find;
		for(i=0;i<pop;i++) 	delete pf[i];
		delete pf;
		delete temp_pf;
		//		
		fout2<<"-- Instance: "<<strTestInstance<<", "<<"the average of IGD is "<<IGD/num_run<<endl;
		*/
	}
	fout2.close();

	lowBound.clear();
	uppBound.clear();
	idealpoint.clear();
}
