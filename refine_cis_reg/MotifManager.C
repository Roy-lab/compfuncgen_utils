#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "MotifManager.H"

MotifManager::MotifManager()
{
}

MotifManager::~MotifManager()
{
}

int 
MotifManager::readMotifFile(const char* aFName)
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
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string motifName;
		string geneName;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				motifName.append(tok);
			}
			else if(tokCnt==1)
			{
				geneName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* geneSet=NULL;
		if(motifGeneInstance.find(motifName)==motifGeneInstance.end())
		{
			geneSet=new map<string,int>;
			motifGeneInstance[motifName]=geneSet;
		}
		else
		{
			geneSet=motifGeneInstance[motifName];
		}
		(*geneSet)[geneName]=0;
		map<string,int>* motifSet=NULL;
		if(geneMotifInstance.find(geneName)==geneMotifInstance.end())
		{
			motifSet=new map<string,int>;
			geneMotifInstance[geneName]=motifSet;
		}
		else
		{
			motifSet=geneMotifInstance[geneName];
		}
		(*motifSet)[motifName]=0;
	}
	inFile.close();
	return 0;
}

map<string,map<string,int>*>&
MotifManager::getAllMotifs()
{
	return motifGeneInstance;
}

map<string,int>* 
MotifManager::getGenesForMotif(const char* aKey)
{
	string key(aKey);
	if(motifGeneInstance.find(key)==motifGeneInstance.end())
	{
		return NULL;
	}	

	return  motifGeneInstance[key];
}

map<string,int>* 
MotifManager::getMotifsForGene(const char* aKey)
{
	string key(aKey);
	if(geneMotifInstance.find(key)==geneMotifInstance.end())
	{
		return NULL;
	}
	return geneMotifInstance[key];
}

int
MotifManager::addFilteredMotif(const char* gene, const char* motif)
{
	string genename(gene);
	string motifname(motif);
	map<string,int>* genes=NULL;
	if(filteredMotifInstance.find(motifname)==filteredMotifInstance.end())
	{
		genes=new map<string,int>;
		filteredMotifInstance[motifname]=genes;
	}
	else
	{
		genes=filteredMotifInstance[motifname];
	}
	(*genes)[genename]=0;
	return 0;
}

map<string,map<string,int>*>&
MotifManager::getFilteredMotifs()
{
	return filteredMotifInstance;
}
