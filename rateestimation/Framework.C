#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <map>
#include "SpeciesDistManager.H"
#include "RateMatrix.H"
#include "Gamma.H"
#include "RateEstimator.H"
#include "Framework.H"


Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readBranchLengths(const char* blenFName)
{
	ifstream inFile(blenFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		double l=0;
		string species;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				species.append(tok);
			}
			else if(tokCnt==1)
			{
				l=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		branchlen[species]=l;
	}
	inFile.close();
	return 0;
}

int 
Framework::initSpeciesTree(const char* treeFName) 
{
	sdMgr.setMaxClusters(inputMapping.size());
	sdMgr.readSpeciesTree(treeFName);
	return 0;
}

int 
Framework::readInputMapping(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)	
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int inputsymbol=0;
		int matrixind=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				inputsymbol=atoi(tok);
			}
			else if(tokCnt==1)
			{
				matrixind=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		inputMapping[inputsymbol]=matrixind;
	}
	inFile.close();
	return 0;

}

int
Framework::readInputData(const char* inFName)
{
	/*int specID=0;
	vector<string> specOrder;
	//sdMgr.getSpeciesListPrefix(specOrder);
	for(int s=0;s<specOrder.size();s++)
	{
		if(strstr(specOrder[s],"Anc")!=NULL)
		{
			continue;
		}
		specIDMap[specID]=specOrder[s];
		specID++;
	}
	while(mFile.good())
	{
		mFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
	}
	mFile.close();*/
	map<int,string> specIDMap;
	int dId=0;
	char buffer[1024];
	ifstream inFile(inFName);
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(dId==0)
		{
			char* tok=strtok(buffer,"\t");	
			int specID=0;
			while(tok!=NULL)
			{
				if(specID>0)
				{
					string specName(tok);
					specIDMap[specID]=specName;
				}
				tok=strtok(NULL,"\t");
				specID++;
			}
		}
		else
		{
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			map<string,int>* dpt=new map<string,int>;
			while(tok!=NULL)
			{
				if(tokCnt>0)
				{
					(*dpt)[specIDMap[tokCnt]]=atoi(tok);
				}	
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
			data[dId-1]=dpt;
		}
		dId++;
	}
	inFile.close();
	return 0;
}

int 
Framework::otherInits(double a, double b)
{
	m.initMatrices(a,b);
	vector<double> pis;
	pis.push_back(2.0/3.0);
	pis.push_back(1.0/3.0);
	m.setPis(pis);
	e.setInitParams(&m,branchlen,inputMapping);
	e.setInitProbs(pis);
	e.setSpeciesDistManager(&sdMgr);
	e.setData(&data);
	
	return 0;
}

int 
Framework::start()
{
	e.learnParams();
	//e.sampleData(100,"samples.txt");
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=7)
	{	
		cout <<"Usage: estimateMatrix blen inputmapping speciesdist data param1 param2" << endl;
		return 0;
	}
	Framework fw;
	fw.readBranchLengths(argv[1]);
	fw.readInputMapping(argv[2]);
	fw.initSpeciesTree(argv[3]);
	fw.readInputData(argv[4]);
	fw.otherInits(atof(argv[5]), atof(argv[6]));
	fw.start();
	return 0;
}
