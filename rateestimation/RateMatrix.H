#ifndef _RATEMATRIX_
#define _RATEMATRIX_
//This class is very specific to the particular problem we are working with
#include <vector>
#include <map>
#include <string>
using namespace std;
class Matrix;
class RateMatrix
{
	public:
	RateMatrix();
	~RateMatrix();
	int initMatrices(double a, double b);
	int setPis(vector<double>&);
	int setPi(int,double);
	Matrix* getQ();
	Matrix* getEigenVectors();
	Matrix* getInvEigenVectors();
	Matrix* getEigenValues();
	double getPi(int);
	typedef struct
	{
		map<int,double> T;
		map<int,map<int,double>*> N;
	} SuffStat;
	int estimateQ(map<string,SuffStat*>&);
	inline double getRate()
	{
		return a_param;
	}
	private:
	Matrix* Q;
	Matrix* V;
	Matrix* Vinv;
	Matrix* L;
	vector<double> pis;
	double a_param, b_param;
};
#endif
