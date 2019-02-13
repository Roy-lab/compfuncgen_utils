#ifndef _RATEESTIMATOR_
#define _RATEESTIMATOR_

#include <map>
#include <string>
#include <gsl/gsl_randist.h>

using namespace std;

class RateMatrix;
class Gamma;
class SpeciesDistManager;
class RateEstimator
{
	public:
	RateEstimator();
	~RateEstimator();
	int setInitParams(RateMatrix*,map<string,double>& blen,map<int,int>&);
	int setInitProbs(vector<double>&);
	int setData(map<int,map<string,int>*>* dataset);
	int setSpeciesDistManager(SpeciesDistManager*);
	int learnParams();
	int sampleData(int datapoints,const char*);
	private:
	//Sampling functions
	int initTransProb();
	int initChildTransitionProb(SpeciesDistManager::Species* species);
	int sampleDataFromNode(gsl_rng* r,int parentID, SpeciesDistManager::Species* node,map<string,int>& assignment);
	int sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs);
	int sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs);

	int initGammas();
	int expectationStep();
	int maximizationStep();
	double estimateTransProb(int i, int j, double t, RateMatrix*);
	int estimateGammas();
	int estimateGammaForDataPt(Gamma::Node* g, double ll);
	int estimateAlphas(Gamma::Node* g,map<string,int>* dataVal);
	int estimateBetas(Gamma::Node* g,map<string,int>* dataVal);
	int clearSuffStats();
	double computeLikelihood();
	double getRootContrib(Gamma*,string&,int);
	double getBranchContrib(Gamma*g,string& ,string& species,int k,int l,double t);
	double computeSufficientStat(int a, int b, int c, int d, double t);
	int clearSufficientStats();
	RateMatrix* Qmat;
	map<int,Gamma*> gammaSet;	
	map<string, double> branchlength;
	map<int,map<string,int>*>* data;
	//Sufficient statistics are calculated for each non-root node and stored
	//here indexed on a species.
	//The final sufficient stats that is used in the M step will require us to sum
	//over all these values.
	map<string,RateMatrix::SuffStat*> sufficientStats;
	int statecnt;
	map<int,int> mappingStateToIndex;
	SpeciesDistManager* sdMgr;
	map<int,double> initprobs;
	int showAlpha;
	gsl_rng* r;
};
#endif
