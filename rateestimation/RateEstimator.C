#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "Matrix.H"
#include "SpeciesDistManager.H"
#include "Gamma.H"
#include "RateMatrix.H"
#include "RateEstimator.H"


RateEstimator::RateEstimator()
{
}

RateEstimator::~RateEstimator()
{
}


int
RateEstimator::setInitProbs(vector<double>& pis)
{
	for(int i=0;i<pis.size();i++)
	{
		initprobs[i]=pis[i];
	}
	return 0;
}

int 
RateEstimator::setInitParams(RateMatrix* m, map<string,double>& blen, map<int,int>& mappings)
{
	Qmat=m;
	statecnt=Qmat->getQ()->getRowCnt();
	for(map<string,double>::iterator bIter=blen.begin();bIter!=blen.end();bIter++)
	{
		branchlength[bIter->first]=bIter->second;
	}
	for(map<int,int>::iterator mIter=mappings.begin();mIter!=mappings.end();mIter++)
	{
		mappingStateToIndex[mIter->first]=mIter->second;
	}
	return 0;
}

int 
RateEstimator::setData(map<int,map<string,int>*>* dataset)
{
	data=dataset;
	return 0;
}
int 
RateEstimator::setSpeciesDistManager(SpeciesDistManager*distMgr)
{
	sdMgr=distMgr;
	return 0;
}

int 
RateEstimator::learnParams()
{
	initGammas();
	int iter=0;
	double currlikelihood=0;
	bool converged=false;
	int maxIter=100;
	while(!converged && iter<maxIter)
	{
		expectationStep();
		maximizationStep();
		double ll=computeLikelihood();
		double diff=0;
		if(iter==0)
		{
			currlikelihood=ll;
		}
		else
		{
			diff=ll-currlikelihood;
			if(fabs(diff)<1e-2)
			{
				converged=true;
			}
			currlikelihood=ll;
		}
		cout <<"Iter: " << iter <<" LL: " << currlikelihood << " "<< diff<< endl;
		iter++;
	}
	return 0;
}


int 
RateEstimator::initTransProb()
{
	SpeciesDistManager::Species* root=sdMgr->getRoot();
	Matrix* conditional=root->getParams();
	int colcnt=conditional->getColCnt();
	for(map<int,double>::iterator aIter=initprobs.begin();aIter!=initprobs.end();aIter++)
	{
		double aval=aIter->second;
		conditional->setValue(aval,0,aIter->first);
	}
	initChildTransitionProb(root->leftchild);
	initChildTransitionProb(root->rightchild);	
	return 0;
}

int
RateEstimator::initChildTransitionProb(SpeciesDistManager::Species* species)
{
	double blen=branchlength[species->name];
	Matrix* param=species->getParams();
	for(map<int,double>::iterator aIter=initprobs.begin();aIter!=initprobs.end();aIter++)
	{
		for(map<int,double>::iterator bIter=initprobs.begin();bIter!=initprobs.end();bIter++)
		{
			double condprob=estimateTransProb(aIter->first,bIter->first,blen,Qmat);
			param->setValue(condprob,aIter->first,bIter->first);
		}
	}
	cout <<"Transition matrix for " << species->name << endl;
	param->showMatrix();
	if(species->leftchild!=NULL)
	{
		initChildTransitionProb(species->leftchild);
	}
	if(species->rightchild!=NULL)
	{
		initChildTransitionProb(species->rightchild);
	}
	return 0;
}

int
RateEstimator::sampleData(int datapoints, const char* aFName)
{
	initTransProb();
	r=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(r,rseed);
	ofstream oFile(aFName);
	
	for(int i=0;i<datapoints;i++)
	{
		double pval=gsl_ran_flat(r,0,1);
		int clustId=-1;
		SpeciesDistManager::Species* root=sdMgr->getRoot();
		map<string,int> assignment;
		Matrix* params=root->getParams();
		vector<int>* sortedClustIDs=root->getSortedClusterIDs(0);
		if(sortedClustIDs==NULL)
		{
			sortedClustIDs=new vector<int>;
			for(int i=0;i<params->getColCnt();i++)
			{
				sortedClustIDs->push_back(i);
			}
			for(int i=0;i<params->getColCnt();i++)
			{
				for(int j=i+1;j<params->getColCnt();j++)
				{
					double v1=params->getValue(0,(*sortedClustIDs)[i]);
					double v2=params->getValue(0,(*sortedClustIDs)[j]);
					if(v1<v2)
					{
						int oldval=(*sortedClustIDs)[i];
						(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
						(*sortedClustIDs)[j]=oldval;
					}
				}
			}
			root->setSortedClusterIDs(0,sortedClustIDs);
		}
		double cdf=params->getValue(0,(*sortedClustIDs)[0]);
		clustId=0;
		while(pval>cdf)
		{
			clustId++;
			cdf=cdf+params->getValue(0,(*sortedClustIDs)[clustId]);
		}
		int actualClustID=(*sortedClustIDs)[clustId];
		assignment[root->name]=actualClustID;
		sampleDataFromNode(r,actualClustID,root->leftchild,assignment);
		sampleDataFromNode(r,actualClustID,root->rightchild,assignment);
		if(i==0)
		{
			oFile <<"Gene";
			for(map<string,int>::iterator aIter=assignment.begin();aIter!=assignment.end();aIter++)
			{
				if(strstr(aIter->first.c_str(),"Anc")!=NULL)
				{
					continue;	
				}
				oFile <<"\t"<< aIter->first;
			}
			oFile << endl;
		}
		oFile<<"Gene"<<i;
		for(map<string,int>::iterator aIter=assignment.begin();aIter!=assignment.end();aIter++)
		{
			if(strstr(aIter->first.c_str(),"Anc")!=NULL)
			{
				continue;	
			}
			oFile <<"\t"<< aIter->second;
		}
		oFile << endl;
		assignment.clear();
	}
	oFile.close();
	return 0;
}

int
RateEstimator::sampleDataFromNode(gsl_rng* r,int parentID, SpeciesDistManager::Species* child,map<string,int>& assignment)
{
	vector<int>* sortedClustIDs=child->getSortedClusterIDs(parentID);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		sortIndices(child->getParams(),parentID,sortedClustIDs);
		child->setSortedClusterIDs(parentID,sortedClustIDs);
	}
	int childID=sampleChildCluster(r,parentID,child->getParams(),sortedClustIDs);
	assignment[child->name]=childID;
	if(child->leftchild!=NULL)
	{
		sampleDataFromNode(r,childID,child->leftchild,assignment);
	}
	if(child->rightchild!=NULL)
	{
		sampleDataFromNode(r,childID,child->rightchild,assignment);
	}
	return childID;
}

int 
RateEstimator::sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs)
{
	int clustId=-1;
	double pval=gsl_ran_flat(r,0,1);
	double cdf=params->getValue(parentClusterId,(*sortedClustIDs)[0]);
	clustId=0;
	while(pval>cdf)
	{
		clustId++;
		cdf=cdf+params->getValue(parentClusterId,(*sortedClustIDs)[clustId]);
	}
	int actualClustID=(*sortedClustIDs)[clustId];
	return actualClustID;
}

int
RateEstimator::sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs)
{
	int colCnt=params->getColCnt();
	for(int i=0;i<colCnt;i++)
	{
		sortedClustIDs->push_back(i);
	}
	for(int i=0;i<params->getColCnt();i++)
	{
		for(int j=i+1;j<params->getColCnt();j++)
		{
			double v1=params->getValue(parentID,(*sortedClustIDs)[i]);
			double v2=params->getValue(parentID,(*sortedClustIDs)[j]);
			if(v1<v2)
			{
				int oldval=(*sortedClustIDs)[i];
				(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
				(*sortedClustIDs)[j]=oldval;
			}
		}
	}
	return 0;
}




int
RateEstimator::initGammas()
{
	SpeciesDistManager::Species* node=sdMgr->getRoot();
	//Have one genetree to initialize all gammas since we will be using essentially the species tree
	for(map<int,map<string,int>*>::iterator dIter=data->begin();dIter!=data->end();dIter++)
	{
		map<string,int>* dpt=dIter->second;
		Gamma* g=new Gamma;
		g->setMaxClusterCnt(statecnt);
		g->initUsingTree(node);
		gammaSet[dIter->first]=g;
	}
	return 0;
}

int 
RateEstimator::expectationStep()
{
	estimateGammas();
	clearSuffStats();
	//Now using each data point we need to estimate the sufficient statistics per data point
	for(map<int,map<string,int>*>::iterator dIter=data->begin();dIter!=data->end();dIter++)
	{
		map<string,int>* dPt=dIter->second;
		Gamma* g=gammaSet[dIter->first];
		int statecnt=Qmat->getQ()->getRowCnt();
		//cout <<"Gamma for " << dIter->first << endl;
		for(map<string,double>::iterator bIter=branchlength.begin();bIter!=branchlength.end();bIter++)
		{
			SpeciesDistManager::Species* species=sdMgr->getSpecies((string&)bIter->first);
			RateMatrix::SuffStat* ss=NULL;
			if(sufficientStats.find(species->name)==sufficientStats.end())
			{
				ss=new RateMatrix::SuffStat;
				sufficientStats[species->name]=ss;
			}
			else
			{
				ss=sufficientStats[species->name];
			}
			//Seems like this is never called, so we will just do this outside of the loop
			if(species->parent==NULL)
			{
				for(int k=0;k<statecnt;k++)
				{
					double v=getRootContrib(g,species->name,k);
					if(isnan(v))
					{
						cout <<"Root is nan for " << dIter->first << endl;
					}
					if(ss->T.find(k)==ss->T.end())
					{
						ss->T[k]=v;
					}
					else
					{
						ss->T[k]=ss->T[k]+v;
					}
				}	
			}
			else
			{
				double t=bIter->second;
				//cout <<"Branch length "<< bIter->first <<"\t" << t << endl;
				//Calculate the N_{k,l}'s which are the expected number of transitions between k and l, for each branch
			//	cout <<"Species " << species->name<< endl;
				//Number of states
				for(int k=0;k<statecnt;k++)
				{
					for(int l=0;l<statecnt;l++)
					{
						double q=Qmat->getQ()->getValue(k,l);
						//This will compute the integrals we need
						double v=getBranchContrib(g,species->name,species->parent->name,k,l,t);
						if(isnan(v))
						{
							cout <<"Value is nan for "<< species->name <<"\t" << species->parent->name <<"\t" <<dIter->first << endl;
						}
						map<int,double>* cnts=NULL;
						if(ss->N.find(k)==ss->N.end())
						{
							cnts=new map<int,double>;
							ss->N[k]=cnts;
						}
						else
						{
							cnts=ss->N[k];
						}
						if(k!=l)
						{
							v=v*q;
						}
						if(cnts->find(l)==cnts->end())
						{
							(*cnts)[l]=v;
						}
						else
						{
							(*cnts)[l]=(*cnts)[l]+v;
						}
					}
				}
			}
		}
		for(int k=0;k<statecnt;k++)
		{
			SpeciesDistManager::Species* root=sdMgr->getRoot();
			RateMatrix::SuffStat* ss=NULL;
			if(sufficientStats.find(root->name)==sufficientStats.end())
			{
				ss=new RateMatrix::SuffStat;
				sufficientStats[root->name]=ss;
			}
			else
			{
				ss=sufficientStats[root->name];
				double v=getRootContrib(g,root->name,k);
				if(isnan(v))
				{
					cout <<"Root is nan for " << dIter->first << endl;
				}
				if(ss->T.find(k)==ss->T.end())
				{
					ss->T[k]=v;
				}
				else
				{
					ss->T[k]=ss->T[k]+v;
				}
			}
		}
	}
	//Now we have computed all our sufficient statistics
	return 0;
}

int 
RateEstimator::maximizationStep()
{
	initprobs.clear();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		Gamma*g=gIter->second;
		double lll=0;
		for(int k=0;k<statecnt;k++)
		{
			if(initprobs.find(k)==initprobs.end())
			{
				initprobs[k]=g->getRoot()->gamma->getValue(0,k);
			}
			else
			{
				initprobs[k]=initprobs[k]+g->getRoot()->gamma->getValue(0,k);
			}
		}
	}
	for(int i=0;i<statecnt;i++)
	{
		initprobs[i]=initprobs[i]/gammaSet.size();
		cout <<"Pi(" <<i  <<") " << initprobs[i] << endl;
	//	Qmat->setPi(i,initprobs[i]);
	}
	
	Qmat->estimateQ(sufficientStats);
	return 0;
}

int
RateEstimator::clearSufficientStats()
{
	for(map<string,RateMatrix::SuffStat*>::iterator sIter=sufficientStats.begin();sIter!=sufficientStats.end();sIter++)
	{
		RateMatrix::SuffStat* ss=sIter->second;
		ss->T.clear();
		for(map<int,map<int,double>*>::iterator qIter=ss->N.begin();qIter!=ss->N.end();qIter++)
		{
			qIter->second->clear();
			delete qIter->second;
		}
	}
	return 0;
}
//This is the probability of transitioning from a to b
double 
RateEstimator::estimateTransProb(int a, int b, double t, RateMatrix* r)
{
	double p=0;
	Matrix* V=r->getEigenVectors();
	Matrix* Vinv=r->getInvEigenVectors();
	Matrix* L=r->getEigenValues();
	for (int i=0;i<V->getRowCnt();i++)
	{
		p=p+(V->getValue(a,i)*Vinv->getValue(i,b)*exp(L->getValue(0,i)*t));
	}
	/*cout << "estimateTransProb values" << endl;
	cout << "V" << endl;
	V->showMatrix(cout);
	cout << "Vinv" << endl;
        Vinv->showMatrix(cout);
	cout << "L" << endl;
	L->showMatrix(cout);
	cout << p << endl;*/
	return p;
}

int
RateEstimator::estimateGammas()
{
	//estimateAlphas
	showAlpha=0;
	for(map<int,map<string,int>*>::iterator dIter=data->begin();dIter!=data->end();dIter++)
	{
		Gamma* g=gammaSet[dIter->first];
		estimateAlphas(g->getRoot(),dIter->second);
		showAlpha=1;
		Gamma::Node* r=g->getRoot();
		//cout <<"gene"<<dIter->first <<"\t" << r->alpha->getValue(0,0) <<" " << r->alpha->getValue(0,1) << " " << r->alpha->getValue(0,2) << " " << r->alpha->getValue(0,3) << endl; 
		estimateBetas(g->getRoot(),dIter->second);
		double ll=0;
		for(int k=0;k<statecnt;k++)
		{
			ll=ll+(g->getRoot()->alpha->getValue(0,k)*Qmat->getPi(k));
		}
		if(g->getRoot()->gamma==NULL)	
		{
			g->getRoot()->gamma=new Matrix(1, statecnt);
		}
		g->getRoot()->gamma->setAllValues(0);
		estimateGammaForDataPt(g->getRoot(),ll);
	}
	return 0;
}

int
RateEstimator::estimateGammaForDataPt(Gamma::Node* g, double ll)
{
	if(g->leftchild!=NULL)
	{
		estimateGammaForDataPt(g->leftchild,ll);
	}
	if(g->rightchild!=NULL)
	{
		estimateGammaForDataPt(g->rightchild,ll);
	}
	if(g->parent==NULL)
	{
		return 0;
	}
	if(g->gamma==NULL)
	{
		g->gamma=new Matrix(statecnt,statecnt);
	}
	double t=branchlength[g->name];
	for(int k=0;k<statecnt;k++)
	{
		double rootp=0;
		for(int l=0;l<statecnt;l++)
		{
			double condprob=estimateTransProb(k,l,t,Qmat);
			double qab=(g->beta->getValue(0,k)*g->alpha->getValue(0,l)*condprob)/ll;
			g->gamma->setValue(qab,k,l);
			rootp=rootp+qab;
		}
		if(g->parent->parent==NULL)
		{
			if(g->parent->gamma->getValue(0,k)>0)
			{
				if(fabs(g->parent->gamma->getValue(0,k)-rootp)>1e-4)
				{
					cout <<"Value mismatch of root  " << g->parent->gamma->getValue(0,k) << " vs " << rootp <<endl;
				}
			}
			g->parent->gamma->setValue(rootp,0,k);	
		}
	}
	//g->gamma->showMatrix();
	return 0;
}

int
RateEstimator::estimateAlphas(Gamma::Node* g,map<string,int>* dataVal)
{
	if(g->alpha==NULL)
	{
		g->alpha=new Matrix(1,statecnt);
	}
	if(g->leftchild==NULL && g->rightchild==NULL)
	{
		int val=(*dataVal)[g->name];
		int mappedVal=mappingStateToIndex[val];
		for(int i=0;i<statecnt;i++)
		{
			if(i==mappedVal)
			{
				g->alpha->setValue(1,0,i);
			}
			else
			{
				g->alpha->setValue(0,0,i);
			}
		}
		return 0;
	}
	if(g->leftchild!=NULL)
	{
		estimateAlphas(g->leftchild,dataVal);
	}
	if(g->rightchild!=NULL)
	{
		estimateAlphas(g->rightchild,dataVal);
	}
	for(int i=0;i<statecnt;i++)
	{
		double leftprod=0;
		double rightprod=0;
		double s=0;
		for(int j=0;j<statecnt;j++)
		{
			double tl=branchlength[g->leftchild->name];
			double transprob=estimateTransProb(i,j,tl,Qmat);
			if(showAlpha==0)
			{
			//	cout <<"Species " << g->leftchild->name <<"\t"<< i << " to " << j << "=" << transprob << endl;
			}
			if(transprob>1.0 || transprob<0)
			{
				cout <<"Wierd prob calculation " << transprob << " for " << i << "->" << j << "\t" << tl << endl;
			}
			s=s+transprob;
			double lval=g->leftchild->alpha->getValue(0,j);
			leftprod=leftprod+(transprob*lval);
			double tr=branchlength[g->rightchild->name];
			transprob=estimateTransProb(i,j,tr,Qmat);
			double rval=g->rightchild->alpha->getValue(0,j);
			rightprod=rightprod+(transprob*rval);
		}
		if(fabs(s-1)>1e-4)
		{
			cout <<"There is a problem for transition probabilities " << s << " for state "<<i << endl;
		}
	//	cout <<"Now stting gamma for " <<g->name << endl;
		g->alpha->setValue(leftprod*rightprod,0,i);
	}
	return 0;
}


int
RateEstimator::estimateBetas(Gamma::Node* g,map<string,int>* dataVal)
{
	if(g->beta==NULL)
	{
		g->beta=new Matrix(1,statecnt);
	}
	if(g->parent==NULL)
	{
		estimateBetas(g->leftchild,dataVal);
		estimateBetas(g->rightchild,dataVal);
		return 0;
	}
	Gamma::Node* sibling=g->parent->leftchild;
	if(sibling==g)
	{
		sibling=g->parent->rightchild;
	}
	if(g->parent->parent==NULL)
	{
		double t=branchlength[sibling->name];
		for(int i=0;i<statecnt;i++)
		{
			double bval=0;
			double pi=Qmat->getPi(i);
			for(int j=0;j<statecnt;j++)
			{
				double pval=estimateTransProb(i,j,t,Qmat);
				bval=bval+(pi*sibling->alpha->getValue(0,j)*pval);
			}
			g->beta->setValue(bval,0,i);
		}
	}
	else
	{
		double t=branchlength[sibling->name];
		double t1=branchlength[g->parent->name];
		for(int i=0;i<statecnt;i++)
		{
			double bval=0;
			for(int j=0;j<statecnt;j++)
			{
				double pval=estimateTransProb(i,j,t,Qmat);
				double v=pval*sibling->alpha->getValue(0,j);
				double s=0;
				for(int k=0;k<statecnt;k++)
				{
					double pval=estimateTransProb(k,i,t1,Qmat);
					s=s+(g->parent->beta->getValue(0,k)*pval);
				}
				bval=bval+(v*s);
			}
			g->beta->setValue(bval,0,i);
		}
	}
	if(g->leftchild!=NULL)
	{
		estimateBetas(g->leftchild,dataVal);
	}
	if(g->rightchild!=NULL)
	{
		estimateBetas(g->rightchild,dataVal);
	}
	return 0;
}

double 
RateEstimator::computeLikelihood()
{
	double ll=0;
	Matrix allSum(Qmat->getQ()->getRowCnt(),Qmat->getQ()->getColCnt());
	SpeciesDistManager::Species* root=sdMgr->getRoot();
	for(int k=0;k<Qmat->getQ()->getRowCnt();k++)
	{
		for(int l=0;l<Qmat->getQ()->getColCnt();l++)
		{
			//Now sum of over the branches in suffstat
			double s=0;
			double s1=0;
			for(map<string,RateMatrix::SuffStat*>::iterator sIter=sufficientStats.begin();sIter!=sufficientStats.end();sIter++)
			{
				if(strcmp(sIter->first.c_str(),root->name.c_str())==0)
				{
					continue;
				}
				double aval=0;
				if(sIter->second->N.find(k)==sIter->second->N.end())
				{
					cout <<"No N stats for " << k << endl;	
					exit(0);
				}
				map<int,double>* vals=sIter->second->N[k];
				s1=s1+(*vals)[l];
				if(k==l)
				{
					//aval=sIter->second->T[l]*Qmat->getQ()->getValue(k,l);
					aval=(*vals)[l]*Qmat->getQ()->getValue(k,l);
				}
				else
				{
					aval=(*vals)[l]*log(Qmat->getQ()->getValue(k,l));
				}
				//This sum is the sum of all avals for a given k,l combination.
				s=s+aval;
			}
			//ll=ll+(s*Qmat->getQ()->getValue(k,l));
			ll=ll+s;
			allSum.setValue(s1,k,l);
		}
	}
	double ll2=0;
	for(int k=0;k<Qmat->getQ()->getRowCnt();k++)
	{
		for(int l=0;l<Qmat->getQ()->getColCnt();l++)
		{
			double sval=allSum.getValue(k,l);
			if(k==l)
			{
				ll2=ll2+(sval*Qmat->getQ()->getValue(k,l));
			}
			else
			{
				ll2=ll2+(sval*log(Qmat->getQ()->getValue(k,l)));
			}
		}
		RateMatrix::SuffStat* rss=sufficientStats[root->name];
		ll2=ll2+rss->T[k];
	}
	

	double ll_alpha=0;
	double ll_gamma=0;
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		Gamma*g=gIter->second;
		double lll_a=0;
		double lll_g=0;
		for(int k=0;k<statecnt;k++)
		{
			lll_a=lll_a+(g->getRoot()->alpha->getValue(0,k)*Qmat->getPi(k));
			lll_g=lll_g+(g->getRoot()->gamma->getValue(0,k)*log(Qmat->getPi(k)));
			///lll=lll+(g->getRoot()->gamma->getValue(0,k)*log(Qmat->getPi(k)));
			if(isnan(lll_a) || isnan(lll_g))
			{
				cout <<"Likelihood became nan for " << gIter->first << "\tk=" << k << endl;
			}
		}
		ll_alpha=ll_alpha+log(lll_a);
		ll_gamma=ll_gamma+lll_g;
		//cout <<"data " << gIter->first << " " << ll_alpha << endl;
	}
	cout <<"Likelihood From Alpha: " << ll_alpha << " Likelihood from suff stat1: " << (ll) << "suff stat2: " << ll2 << endl;
	return ll_alpha;
	//return ll+ll_gamma;
}

int
RateEstimator::clearSuffStats()
{
	for(map<string,RateMatrix::SuffStat*>::iterator sIter=sufficientStats.begin();sIter!=sufficientStats.end();sIter++)
	{
		sIter->second->T.clear();
		map<int,map<int,double>*>& N=sIter->second->N;
		for(map<int,map<int,double>*>::iterator nIter=N.begin();nIter!=N.end();nIter++)
		{
			nIter->second->clear();
			delete nIter->second;
		}
		N.clear();
		delete sIter->second;
	}
	return 0;
}


double 
RateEstimator::getRootContrib(Gamma* g,string& spec,int stateid)
{	
	//cout <<"In root " <<endl;
	Matrix* m=g->getRoot()->alpha;
	double pval=m->getValue(0,stateid);
	pval=pval*Qmat->getPi(stateid);
	pval=log(pval);
	return pval;
}

double 
RateEstimator::getBranchContrib(Gamma*g,string& spec,string& parentspec,int k,int l, double t)
{
	//Note k and l are the intermediates.
	Matrix* Q=Qmat->getQ();
	Matrix* m=g->getGamma(spec);
	double pval=0;
	for(int a=0;a<Q->getRowCnt();a++)
	{
		for(int b=0;b<Q->getColCnt();b++)
		{
			double val=m->getValue(a,b);
			//double R=computeSufficientStat(k,l,a,b,t);
			double R=computeSufficientStat(a,k,l,b,t);
			double transprob=estimateTransProb(a,b,t,Qmat);
			pval=pval+((R*val)/transprob);
		}
	}	
	return pval;
}

//This is the computeSS function from mathematica
/*double
RateEstimator::computeSufficientStat(int a, int b, int c, int d, double t)
{
	Matrix* V=Qmat->getEigenVectors();
	Matrix* L=Qmat->getEigenValues();
	double ss=0;
	for(int i=0;i<V->getRowCnt();i++)
	{
		double f=V->getValue(a,i)*V->getValue(b,i);
		double rowsum=0;
		for(int j=0;j<V->getColCnt();j++)
		{
			double e=0;
			double l1=L->getValue(0,i);
			double l2=L->getValue(0,j);
			if(fabs(l1-l2)<1e-10)
			{
				e=t*exp(L->getValue(0,i)*t);
			}
			else
			{
				e=(exp(l1*t)-exp(l2*t))/(l1-l2);
				if(isnan(e))
				{
					cout <<"Found nan for suff stat " << i <<"\t" << j << endl;
				}
			}
			double g=V->getValue(c,j)*V->getValue(d,j)*e;
			rowsum=rowsum+g;
		}
		ss=ss+(f*rowsum);
	}
	double transprob=estimateTransProb(a,d,t,Qmat);
	//ss=ss/transprob;
	return ss;
}*/

 
//This is the computeSS function from mathematica
//a is a, b is i (and also k from above), c is j and also l, and d is b
double
RateEstimator::computeSufficientStat(int a, int i, int j, int b, double t)
{
	Matrix* V=Qmat->getEigenVectors();
	Matrix* Vinv=Qmat->getInvEigenVectors();
	Matrix* L=Qmat->getEigenValues();
	double ss=0;
	for(int k=0;k<V->getRowCnt();k++)
	{
		double f=V->getValue(a,k)*Vinv->getValue(k,i);
		double rowsum=0;
		for(int l=0;l<V->getColCnt();l++)
		{
			double e=0;
			double l1=L->getValue(0,k);
			double l2=L->getValue(0,l);
			if(fabs(l1-l2)<1e-10)
			{
				e=t*exp(L->getValue(0,k)*t);
			}
			else
			{
				e=(exp(l1*t)-exp(l2*t))/(l1-l2);
				if(isnan(e))
				{
					cout <<"Found nan for suff stat " << i <<"\t" << j << endl;
				}
			}
			double g=V->getValue(j,l)*Vinv->getValue(l,b)*e;
			rowsum=rowsum+g;
		}
		ss=ss+(f*rowsum);
	}
	return ss;
} 
