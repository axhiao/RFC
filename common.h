#ifndef __COMMON_H_
#define __COMMON_H_

#include "global.h"
#define PI  3.1415926535897932384626433832795

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
	 // 数组x的前n个元素参与比较，选出其中前m小，标号的重新排列放入数组idx中
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}

double dist(double *pf,double *pf_find)
{
	double temp  = 0;
	for (int i=0; i<nobj; i++)
	{
		temp += (pf[i]-pf_find[i])*(pf[i]-pf_find[i]);
	}
	return temp;
}

int find_point(double **pf,vector<bool> flag,double **pf_find,int num,int pop)
{
	int i,j;
	double max_distance = -1e9;
	int max_index = -1;
	for (i=0; i<pop; i++)
	{
		if (flag[i] == true)
		{
			double distance_point = 1e10;
			for (j=0; j<num; j++)
			{
				double distance = dist(pf[i],pf_find[j]);
				if (distance <distance_point)
				{
					distance_point = distance;
				}
			}
			if (distance_point > max_distance)
			{
				max_distance = distance_point;
				max_index = i;
			}
		}
	}
	return max_index;
}

int donimate_judge(vector<double> pf1, vector<double> pf2)
{
	//return 1:  pf2 donimate pf1;
	//return -1: pf1 donimate pf2;
	//return 0:  pf1,pf2 is non donimate;
	double big=0,small=0;
	for (int i=0; i<nobj; i++)
	{
		if (pf1[i]>=pf2[i])	big++;
		if (pf1[i]<=pf2[i]) small++;		
	}
	if (small == nobj)	return -1;
	else if (big == nobj) return 1;
	else return 0;
}

int donimate_relation(vector<double> pf1, vector<double> pf2)
{
	//return 2:  pf1 = pf2
	//return 1:  pf1 be strong donimated pf2;
	//return -1: pf1 strong donimate pf2;
	//return 0:  pf1,pf2 is non donimate;
	double big=0,small=0;
	for (int i=0; i<nobj; i++)
	{
		if (pf1[i]>=pf2[i])	big++;
		if (pf1[i]<=pf2[i]) small++;		
	}
	if (small == nobj && big == nobj) return 2;
	else if (small == nobj)	return -1;
	else if (big == nobj) return 1;
	else return 0;
}

double dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double dist_vector_square(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sum;
}

double norm(vector <double> vec)
{
	int dim = vec.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum += vec[n] * vec[n];
	return sqrt(sum);
}

double innerproduct(vector <double> x, vector <double> y)
{
	int dim = x.size();
	double sum = 0;
	for (int n = 0; n < dim; n++)	sum = sum + x[n] * y[n];
	return sum;
}

bool   dominate(vector<double> &u, vector<double> &v, double epsilon)
{
    int dim = u.size();
	for(int i=0; i<dim; i++)
	{
	    if(u[i]<v[i]-epsilon)
		{
		    return false;
		}
	}
	return true;
}

double fitnessfunction(vector <double> &y_obj, vector <double> &namda)
{
    double fvalue = 0;    
	// Tchebycheff approach
	if(!strcmp(strFunctionType,"_TCH1"))
	{
		double max_fun = -1.0e+30;
		for(int n=0; n<nobj; n++)
		{
			double diff = fabs(y_obj[n] - idealpoint[n] );
			double feval;
			if(namda[n]==0) 
				feval = 0.00001*diff;
			else
			    feval = diff*namda[n];
			if(feval>max_fun) max_fun = feval;

		}
		
		fvalue = max_fun;
		return fvalue;
	}

	// normalized Tchebycheff approach
	if(!strcmp(strFunctionType,"_TCH2"))
	{
		double max_fun = -1.0e+30;
		for(int n=0; n<nobj; n++)
		{
			double diff = (y_obj[n] - idealpoint[n])/(nad_point[n]- idealpoint[n]);
			double feval;
			if(namda[n]==0) 
				feval = 0.0001*diff;
			else
			    feval = diff*namda[n];
			if(feval>max_fun) max_fun = feval;

		}
		fvalue = max_fun;
		return fvalue;
	}
/*	// normalized Tchebycheff approach
	if(!strcmp(strFunctionType,"_TCH2"))
	{
		vector <double> scale;
		for(int i=0; i<nvar; i++)
		{
			double min = 1.0e+30, max = -1.0e+30;
			for(int j=0; j<nvar; j++)
			{
				double tp = nbi_node[j].y_obj[i];
				if(tp>max) max = tp;
				if(tp<min) min = tp;
			}
			scale.push_back(max-min);
			if(max-min==0) return 1.0e+30;
		}

		double max_fun = -1.0e+30;
		for(int n=0; n<nvar; n++)
		{
			double diff = (y_obj[n] - idealpoint[n])/scale[n];
			double feval;
			if(namda[n]==0) 
				feval = 0.0001*diff;
			else
			    feval = diff*namda[n];
			if(feval>max_fun) max_fun = feval;

		}
		fvalue = max_fun;
		return fvalue;
	}
*/

	//* Boundary intersection approach
	if(!strcmp(strFunctionType,"_PBI"))
	{
		// normalize the weight vector (line segment)
		double nd = norm(namda);
		int n;
		for(int i=0; i<nobj; i++)
			namda[i] = namda[i]/nd;

		vector <double> realA(nobj,0);
		vector <double> realB(nobj,0);

		// difference beween current point and reference point
		for(n=0; n<nobj; n++)
			realA[n] = (y_obj[n] - idealpoint[n]);

		// distance along the line segment
		double d1 = fabs(innerproduct(realA,namda));

		// distance to the line segment
		for(n=0; n<nobj; n++)
			realB[n] = (y_obj[n] - (idealpoint[n] + d1*namda[n]));
		double d2 = norm(realB);

		fvalue = d1 + 25*d2;
		realA.clear();
		realB.clear();
		return fvalue;
	}

	// Weight Sum approach
	if(!strcmp(strFunctionType,"_WS"))
	{		
		for(int n=0; n<nobj; n++)
		{
			fvalue += namda[n]*y_obj[n];
		}	
		
		return fvalue;
	}	
}

double achievement_scalarizing_functions(vector<double> weight, vector<double> f, vector<double> reference)
{
	int i;
	double maxsum = -1e30, sum = 0.0, tmp = 0.0;
	double lou = 1e-6;
	for (i=0; i < nobj; i++)
	{
		tmp = weight[i]*(f[i] - reference[i]);
		if (tmp > maxsum)	maxsum = tmp;
		sum += f[i] - reference[i];
	}
	return maxsum + lou * sum;
}

double InformationLoss(vector<double> x)
{
	double u1, u2, u3;
	u1 = 3 * (1 - x[0]) * (1-x[0]) * exp(-x[0] * x[0] - (x[1] + 1) * (x[1] + 1));
	u2 = -10 * (0.25 * x[0] - pow(x[0],3) - pow(x[1],5)) * exp(-x[0] * x[0] - x[1] * x[1]);
	u3 = 1.0 / 3 * exp(-(x[0] + 1) * (x[0] + 1) - x[1] * x[1]);
	return -u1 - u2 - u3 + 10;
}

void load_data(char filename[1024], vector< vector<double> > &data, int dim)
{
	std::ifstream readf(filename);
	vector<double> vect = vector<double>(dim, 0);
	while(!readf.eof())
	{
        for(int i=0; i<dim; i++)
		{
			readf>>vect[i];
			//printf(" %f ", vect[i]);
		}
		data.push_back(vect);
		//vect.clear();    // bug here. 
	}
	vect.clear();
	readf.close();    
}

int loadFData(char * a)
{
	//从文件读入库洪水量记录**************************
	Idata.resize(0);  //入库洪水量记录
	int dim=0;
	char filename[1024];
    sprintf(filename,a); 
    std::fstream fin;
	fin.open(filename,std::ios::in);
	if(fin.is_open())
	{
		const char* pch2;
		char  q[20];
		std::string str;
		int cur=0;
		double Q=0;
		while(!fin.eof())
		{
			std::getline(fin,str,'\n');
			if(str=="")
				continue;
			pch2 = str.c_str();
			sscanf(pch2,"%s",q);
			cur++;
			if(T==1)
			{
				Q=atof(q);	
				Idata.push_back(Q);
				dim++;	
			}	
			else
			{
				if(cur%T==1)
				{
					Q=atof(q);	
					Idata.push_back(Q);
					dim++;			
				}
			}		
		}
	} //end if
	else
		std::cout<<"failed to open "<<filename<<endl;
    fin.close();

	//文件输出入库流量数据
	std::fstream fout;
	char saveQdataFilename[100];
	sprintf(saveQdataFilename,"Qdata_Sample/%s",a);
	fout.open(saveQdataFilename,std::ios::out);
	for(int i=0; i<dim; i++)
	{
		fout<<Idata[i]<<endl;
	}

	return dim;
}

void load_L2C_Curve(char * a)
{
	L2C.resize(2); //存储水位容量曲线： L2C[0]:水位  L2C[1]:对应该水位的水库容量
	L2C[0].resize(0);
	L2C[1].resize(0);
	char filename[1024];
    sprintf(filename,a); 
    std::fstream fin;
	fin.open(filename,std::ios::in);
	if(fin.is_open())
	{
		const char* pch2;
		char  L[20],C[20];
		std::string str;
		double wl=0,cp=0;
		while(!fin.eof())
		{
			std::getline(fin,str,'\n');
			if(str=="")
				continue;
			pch2 = str.c_str();
			//cout<<pch2<<endl;
			sscanf(pch2,"%s %s",L,C);
            wl=atof(L);
			//cout<<wl<<endl;
			L2C[0].push_back(wl);
			cp=1000000*atof(C);     //?????????????????????????????????????????????
			L2C[1].push_back(cp);
		}
	}
}


void preprocess_for_reservior()
{
	if (strcmp("ankang20100715", strTestInstance) && strcmp("ankang20001012", strTestInstance) && strcmp("ankang20030828", strTestInstance) && strcmp("ankang20051001", strTestInstance))
		return;

	if(!strcmp("ankang20100715", strTestInstance)) T=6;
	if(!strcmp("ankang20001012", strTestInstance)) T=6;
	if(!strcmp("ankang20030828", strTestInstance)) T=3;
	if(!strcmp("ankang20051001", strTestInstance)) T=4;

	char QFile[100];
	sprintf(QFile,"Qdata/%s.txt",strTestInstance);	
	//char QFile[100]="Qdata/ankang20100715.txt";
	int Qsize = loadFData(QFile);         // 载入入库流量 quantity of flow 数据
	char L2C_File[20]="Qdata/L2C.txt";
	load_L2C_Curve(L2C_File);             // 载入水位level容量capacity曲线数据
	
	nvar  = Qsize - 1;  //dimensionality of search space 是对moea.cpp中nvar的修订
	nobj  = 2;  //number of objective functions
}


double max(vector<double> & var) //求向量最大值
{
	int dim=var.size();
	double max=var[0];
	for(int i=1; i<dim; i++)
	{
		if(max<var[i])
			max=var[i];
	}
	return max;
}


double Get_WL_from_CP(double cp) // get water level from capacity 从水库容量获得水位
{
	//-0.0000037331 x.^2 + 0.031249 x +274.49  x(百万立方米)
	double wl;
	double x;
	if(L2C.size()==0)
	{
		cout<<"无法获得水库水位！"<<endl;
		return -1;
	}
	if(cp < L2C[1][0])
	{
		//cout<<"水库容量过小，无法获得水位值！"<<endl;
		//return -1;
		x=cp/1000000;
		wl=-0.0000037331*x*x + 0.031249* x +274.49;
		return wl;
	}
	int size=L2C[1].size();
	int cur=0;
	while(cur<=size-1 && cp>L2C[1][cur])
	{
		cur++;
	}
	if(cur>=size)
	{
	    //cout<<"水库容量过大，无法获得水位值！"<<endl;
		//return -1;
		x=cp/1000000;
		wl=-0.0000037331*x*x + 0.031249* x +274.49;
	}
	else
	{
		wl=L2C[0][cur];
	}
	return wl;
}

double Get_CP_from_WL(double wl) //get capacity from water level 从水库水位获得容量
{
	if(L2C.size()==0)
	{
		cout<<"无法获得水库容量！"<<endl;
		return -1;
	}
	if(wl<L2C[0][0])
	{
		cout<<"水库水位过低，无法获得容量值！"<<endl;
		return -1;
	}
	int size=L2C[0].size();
	int cur=0;
	while(wl>L2C[0][cur] && cur<size-1)
	{
		cur++;
	}
	if(cur>size)
	{
	    cout<<"水库水位过高，无法获得容量值！"<<endl;
		return -1;
	}
	double cp=L2C[1][cur];
	return cp;
}

double GetmaxZ(vector<double> &x_var, vector<double> &xZ,double & EQ)//最大水位
//x_var 下泄流量 x_var：下泄流量 我们唯一可以调节的，能动地改变之
//xZ    水位向量
//EQ    超出水位约束的总量
{
	//根据下泄流量计算水位
	int j;
	EQ=0;
	int dim = x_var.size();  //dim = x_var.size() = nvar (individual.h 51)  nvar  = Qsize - 1(411); 
	if(dim != Idata.size()-1 )
	{
		cout<<"x_var.size()错误！"<<endl;
	    return -1;
	}
	xZ.clear();
	xZ.resize(dim);
	double Vt1=V_0, Vt2;        // double  V_0 = 1737250000; 初始调度水库容量  
	double I1=Idata[0], I2=Idata[1];   //入库洪水流量
	double Q1=Q_0, Q2=x_var[0]; //double  Q_0 = 773; 初始下泄流量
	Vt2=Vt1+0.5*(I2+I1)*(T*3600)-0.5*(Q1+Q2)*(T*3600); //0.5???????????????????????????
	xZ[0]=Get_WL_from_CP(Vt2);
//cout <<xZ[0]<<" & ";
	if( xZ[0] > ZU )
		EQ+=xZ[0]-ZU;
	if( xZ[0] < ZB )
		EQ+=ZB-xZ[0];
	for( int i=1; i<dim; i++)
	{
		I1=Idata[i];
		I2=Idata[i+1];
		Q1=x_var[i-1];
		Q2=x_var[i];
		Vt1=Vt2;
		Vt2=Vt1+0.5*(I2+I1)*(T*3600)-0.5*(Q1+Q2)*(T*3600);
		if (i==12)
		{
			j=0;
		}
		xZ[i]=Get_WL_from_CP(Vt2);
//	cout <<xZ[i]<<" & ";
		if( xZ[i] > ZU )
			EQ+=xZ[i]-ZU;
		if( xZ[i] < ZB )
			EQ+=ZB-xZ[i];
	}	
	//调度期末水位约束
	double dz=xZ[dim-1]-ZFL;
	if(dz<0)
		dz=-dz;
	
	if(dz<10)
		dz=0;
	//	else
	//		dz=dz-5;
	
	if(EQ>0)
	{
		//		cout <<EQ<<" & ";
		return ZU+EQ/nvar;//+1.0*dz; // dz的系数越大，PF越窄//???????????????????????????????????????????????????
	}
	else
		return max(xZ);//+1.0*dz;//?????????????????????????????????????????????????????????????????????????
}

//----------------------------------------------//
/* the following function designed for  WFG*/
//----------------------------------------------//
// the following function designed for shape for object-function//
double Convex(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= 1 - cos(x_var[i-1]*PI/2); 
	if (m == 1)			return result;
	else				return result *= 1-sin(x_var[nobj - m + 1 - 1]*PI/2);		
}

double Linear(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= x_var[i-1]; 
	if (m == 1)			return result;
	else				return result *= 1-x_var[nobj - m + 1 - 1];
}

double Concave(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= sin(x_var[i-1]*PI/2); 
	if (m == 1)			return result;
	else				return result *= cos(x_var[nobj - m + 1 - 1]*PI/2);		
}

double Mixed(vector<double> x_var, double alfa, int A)
{//alfa > 0,A subject to {1,2,3,...}
	return pow(1 - x_var[0] - (cos(2*A*PI + PI/2)) / (2*A*PI), alfa);	
}

double Disc(vector<double> x_var, double alfa, double belta, int A)
{//alfa > 0, belta > 0, A subject to {1,2,3,...}
	return 1 - pow(x_var[0], alfa) * cos(A * pow(x_var[0], belta) * PI)* cos(A * pow(x_var[0], belta) * PI);	
}

//----------------------------------------------//
// the following function designed for Transformation functions//
double b_poly(double y, double alfa)
{//alfa > 0, alfa not equal to 1
	return pow(y, alfa);	
}

double b_flat(double y, double A, double B, double C)
{//A,B,C subject to [0,1]
	double temp1,temp2;
	if (floor(y - B) <= 0)	temp1 = floor(y - B);
	else					temp1 = 0;
	if (floor(C - y) <= 0)	temp2 = floor(C - y);
	else					temp2 = 0;
	return A + temp1 * A* (B - y) / B - temp2 * (1 - A) * (y - C) / (1 - C);	
}

double b_param(double y, double u_y, double A, double B, double C)
{//A subject to [0,1], 0 < B < C	
	return pow(y, B + (C - B) * (A - (1 - 2 * u_y) * fabs(floor(0.5 - u_y) + A) ));	
}

double s_linear(double y, double A)
{//A subject to [0,1]	
	return fabs(y - A) / fabs(floor(A - y) + A);	
}

double s_multi(double y, int A, double B, double C)
{//A subject to {1,2,...}, B>=0, (4*A+2)*PI>=4*B, C subject to [0,1]	
	return (1 + cos((4*A + 2) * PI * (0.5 - fabs(y-C)/(2*(floor(C - y) + C)))) + 4 * B * fabs(y-C)/(2*(floor(C-y)+C)) *fabs(y-C)/(2*(floor(C-y)+C)) ) / (B + 2);	
}

double s_decept(double y, double A, double B, double C)
{//A subject to [0,1], 0<=B<<1, 0<=C<<1, A-B>0,A+B<1	
	return 1+(fabs(y-A)-B)*((floor(y-A+B)*(1-C+(A-B)/B))/(A-B)+(floor(A+B-y)*(1-C+(1-A-B)/B))/(1-A-B)+1.0/B);
}

double r_sum(vector<double> y, vector<double> w)
{//|w| = |y|, w1,w2,...w|y| > 0
	double total_w = 0;
	double weight_sum = 0;
	int i;
	for (i = 0; i < w.size(); i++)
	{
		total_w += w[i];
		weight_sum += w[i] * y[i];
	}
	return weight_sum/total_w;	
}

double r_nonsep(vector<double> y, double A)
{//A subject to {1,...,|y|}, |y| mod A = 0
	double numerator = 0.0;
	
	for( int j = 0; j < y.size(); j++ )
	{
		numerator += y[j];

		for( int k = 0; k <= A-2; k++ )
		{
			numerator += fabs( y[j] - y[( j+k+1 ) % y.size()] );
		}
	}

	const double tmp = ceil( A/2.0 );
	const double denominator = y.size()*tmp*( 1.0 + 2.0*A - 2.0*tmp )/A;

	return  numerator / denominator;
}

vector<double>  WFG_normalise_z(vector<double> x)
{	
	int i;
	vector<double> z(nvar, 0);
	for (i = 1; i <= x.size(); i++)		z[i-1] = x[i-1]*1.0 / (2*i);	
	return z;
}

vector<double> WFG1_t1(vector<double> y_old, int k)
{	
	int i;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)	y_new[i-1] = y_old[i-1];
	for (i = k + 1; i <= nvar; i++)		y_new[i-1] = s_linear(y_old[i-1],0.35);
	
	return y_new;
}

vector<double> WFG1_t2(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)	y_new[i-1] = y_old[i-1];
	for (i = k + 1; i <= nvar; i++)		y_new[i-1] = b_flat(y_old[i-1],0.8,0.75,0.85);
	
	return y_new;
}

vector<double> WFG1_t3(vector<double> y_old)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= nvar; i++)	y_new[i-1] = b_poly(y_old[i-1],0.02);	
	
	return y_new;
}

vector<double> WFG1_t3_var(vector<double> y_old)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= nvar; i++)	y_new[i-1] = b_poly(y_old[i-1],0.8);	
	
	return y_new;
}

vector<double> WFG1_t4(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;
	vector<double> w;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(2*j);
		}
		y_new[i-1] = r_sum(temp_y, w);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j < nvar; j++)
	{
		temp_y.push_back(y_old[j-1]);
		w.push_back(2*j);
	}
	y_new[nobj - 1] = r_sum(temp_y, w);

	return y_new;
}

vector<double> WFG1_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Mixed(x,1,5);

	return y_new;
}
vector<double> WFG1_VAR_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Mixed(x,7,8);

	return y_new;
}


vector<double> WFG2_t2(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new;
	int l = nvar - k;
	
	for (i = 1; i <= k ; i++) y_new.push_back( y_old[i-1] );
	for (i = k + 1;  i <= k+l/2; i++)
	{
		vector<double> temp;
		temp.push_back(y_old[k+2*(i-k) - 1-1]);
		temp.push_back(y_old[k+2*(i-k)-1]);
		y_new.push_back(r_nonsep(temp, 2));

		temp.clear();
	}

	return y_new;
}

vector<double> WFG2_t3(vector<double> y_old, int k)
{	
	int i,j;
	int l = nvar - k;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;
	vector<double> w;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = r_sum(temp_y, w);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j <= k + l/2; j++)
	{
		temp_y.push_back(y_old[j-1]);
		w.push_back(1);
	}
	y_new[nobj - 1] = r_sum(temp_y, w);

	return y_new;
}

vector<double> WFG2_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Disc(x,1,1,5);

	return y_new;
}

vector<double> WFG3_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	vector<double> A(nobj, 0);
	A[0] = 1;
	
	//turn to x
	for (i = 1; i < nobj; i++)
	{
		double max;
		if (y_old[nobj - 1] < A[i-1])	max = 1;
		else						max = y_old[nobj - 1];
		x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	}
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i <= nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Linear(x,i);

	return y_new;
}

vector<double> WFG4_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i < nvar; i++) y_new[i - 1] = s_multi(y_old[i - 1], 30, 10, 0.35);

	return y_new;
}

vector<double> WFG4_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	vector<double> A(nobj, 1);
	
	//turn to x
	for (i = 1; i < nobj; i++)
	{
		double max;
		if (y_old[nobj - 1] < A[i-1])	max = A[i-1];
		else							max = y_old[nobj - 1];
		x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	}
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i <= nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Concave(x,i);

	return y_new;
}

vector<double> WFG5_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i < nvar; i++) y_new[i - 1] = s_decept(y_old[i - 1], 0.35, 0.001, 0.05);

	return y_new;
}

vector<double> WFG6_t2(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
		}
		y_new[i-1] = r_nonsep(temp_y, k*1.0/(nobj - 1));
		
		temp_y.clear();
	}
	for (j = k+1; j <= nvar; j++)		temp_y.push_back(y_old[j-1]);

	y_new[nobj - 1] = r_nonsep(temp_y, nvar - k);

	return y_new;
}

vector<double> WFG7_t1(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	for (i = 1; i <= k; i++)	
	{
		for (j = i + 1; j <= nvar; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j <= nvar; j++)		y_new[j - 1] = y_old[j-1];	

	return y_new;
}

vector<double> WFG8_t1(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	
	for (j = 1; j <= k; j++)		y_new[j - 1] = y_old[j-1];
	for (i = k+1; i <= nvar; i++)	
	{
		for (j = 1; j <= i - 1; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}	

	return y_new;
}

vector<double> WFG9_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	for (i = 1; i <= nvar - 1; i++)	
	{
		for (j = i + 1; j <= nvar; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}
	y_new[nvar - 1] = y_old[nvar - 1];	

	return y_new;
}

vector<double> WFG9_t2(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)		y_new[i - 1] = s_decept(y_old[ i - 1 ], 0.35, 0.001, 0.05);
	for (i = k + 1; i <= nvar; i++)	y_new[i - 1] = s_multi(y_old[i - 1],30, 95, 0.35);	

	return y_new;
}

#endif