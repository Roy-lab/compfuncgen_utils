#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
using namespace std;
class RateMatrix;
class RateEstimator;
class SpeciesDistManager;
class Framework
{
	public:
	Framework();
	~Framework();
	int readBranchLengths(const char*);
	int initSpeciesTree(const char*);
	int readInputMapping(const char*);
	int readInputData(const char*);
	int otherInits(double a, double b);
	int start();
	private:
	RateMatrix m;
	RateEstimator e;
	SpeciesDistManager sdMgr;
	map<int,int> inputMapping;
	map<string,double> branchlen;
	map<int,map<string,int>*> data;
};

#endif
