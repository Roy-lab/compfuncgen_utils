#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>

#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "SpeciesNode.H"
#include "ParsimonyInferrer.H"
#include "MotifManager.H"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readMotifs(const char* aFName)
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
		char fName[1024];
		string specName;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				strcpy(fName,tok);
			}
			else
			{
				specName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		MotifManager* mm=new MotifManager;
		mm->readMotifFile(fName);
		motifCollection[specName]=mm;
		map<string,map<string,int>*>& motifs=mm->getAllMotifs();
		for(map<string,map<string,int>*>::iterator mIter=motifs.begin();mIter!=motifs.end();mIter++)
		{
			motifNames[mIter->first]=0;
		}
	}
	inFile.close();
	return 0;
}

int 
Framework::readTree(const char* aFName)
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
		string child;
		string parent;
		string branch;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				child.append(tok);
			}
			if(tokCnt==1)
			{
				branch.append(tok);
			}
			else if(tokCnt==2)
			{
				parent.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		SpeciesNode* parentNode=NULL;
		if(nodeSet.find(parent)==nodeSet.end())
		{
			parentNode=new SpeciesNode;
			parentNode->name.append(parent);
			parentNode->left=NULL;
			parentNode->right=NULL;
			parentNode->parent=NULL;
			nodeSet[parent]=parentNode;
		}	
		else
		{
			parentNode=nodeSet[parent];
		}
		SpeciesNode* childNode=NULL;
		if(nodeSet.find(child)==nodeSet.end())
		{
			childNode=new SpeciesNode;
			childNode->name.append(child);
			childNode->left=NULL;
			childNode->right=NULL;
			childNode->parent=NULL;
			nodeSet[child]=childNode;
		}	
		else
		{
			childNode=nodeSet[child];
		}
		childNode->parent=parentNode;
		if(strcmp(branch.c_str(),"left")==0)
		{
			parentNode->left=childNode;
		}
		else if(strcmp(branch.c_str(),"right")==0)
		{
			parentNode->right=childNode;
		}
	}
	for(map<string,SpeciesNode*>::iterator nIter=nodeSet.begin();nIter!=nodeSet.end();nIter++)
	{
		SpeciesNode* n=nIter->second;
		if(n->parent==NULL)
		{
			root=n;
		}
	}
	inFile.close();
	return 0;
}

int 
Framework::readSpeciesOrthology(const char* speciesOrder, const char* orthoMapFName)
{
	mor.readSpeciesMapping(speciesOrder);
	mor.readFile(orthoMapFName);
	return 0;
}

int 
Framework::setScore(const char* sFName)
{
	pi.setScore(sFName);
	return 0;
}

int 
Framework::refineTargets(const char* aFName)
{
	//Start from the set of orthogroups
	map<int,MappedOrthogroup*>& orthogroupSet=mor.getMappedOrthogroups();
	//Do this for all motifs
	pi.setSpeciesTree(root);
	for(map<string,int>::iterator moIter=motifNames.begin();moIter!=motifNames.end();moIter++)
	{
		for(map<int,MappedOrthogroup*>::iterator mIter=orthogroupSet.begin();mIter!=orthogroupSet.end();mIter++)
		{
			MappedOrthogroup* mog=mIter->second;
			map<int,map<string,string>*>& geneSets=mog->getGeneSets();
			map<string,int> motifassign;
			for(map<int,map<string,string>*>::iterator sIter=geneSets.begin();sIter!=geneSets.end();sIter++)
			{
				map<string,string>* genesinspecies=sIter->second;
				for(map<string,string>::iterator aIter=genesinspecies->begin();aIter!=genesinspecies->end();aIter++)
				{
					MotifManager* mm=motifCollection[aIter->first];	
					map<string,int>* motifSet=mm->getGenesForMotif(moIter->first.c_str());
					if(motifSet==NULL)
					{
						motifassign[aIter->first]=0;
						continue;
					}
					if(motifSet->find(aIter->second)==motifSet->end())
					{
						motifassign[aIter->first]=0;
						continue;
					}
					motifassign[aIter->first]=1;
				}
				map<string,int> ancestralassign;
				pi.inferAncestralAssignment(motifassign,ancestralassign);
				//cout <<mIter->first << endl;
				//showAssign(ancestralassign,root);
				map<string,int> filteredMotifs; 
				getAssignFiltered(ancestralassign,root,false,filteredMotifs);
				for(map<string,int>::iterator fIter=filteredMotifs.begin();fIter!=filteredMotifs.end();fIter++)
				{
					MotifManager* mm=motifCollection[fIter->first];	
					mm->addFilteredMotif((*genesinspecies)[fIter->first].c_str(),moIter->first.c_str());
				}
			}
		}
	}
	for(map<string,MotifManager*>::iterator sIter=motifCollection.begin();sIter!=motifCollection.end();sIter++)
	{
		MotifManager* mm=sIter->second;
		map<string,map<string,int>*>& motifSet=mm->getFilteredMotifs();
		char outFName[1024];
		sprintf(outFName,"%s_%s_regnet.txt",aFName,sIter->first.c_str());
		ofstream oFile(outFName);
		for(map<string,map<string,int>*>::iterator mIter=motifSet.begin();mIter!=motifSet.end();mIter++)
		{
			map<string,int>* genes=mIter->second;
			for(map<string,int>::iterator gIter=genes->begin();gIter!=genes->end();gIter++)
			{
				oFile <<gIter->first<<"\t" <<mIter->first << endl;
			}
		}
		oFile.close();
	}
	return 0;
}

int
Framework::showAssign(map<string,int>& assignment,SpeciesNode* tree)
{
	if(tree->left==NULL && tree->right==NULL)
	{
		if(assignment.find(tree->name)==assignment.end())
		{
			cout <<tree->name <<"\t" << 0 << endl;;
		}
		else
		{
			cout <<tree->name <<"\t"<< assignment[tree->name] << endl;
		}
	}
	else
	{
		if(tree->left!=NULL)
		{
			showAssign(assignment,tree->left);
		}
		if(tree->right!=NULL)
		{
			showAssign(assignment,tree->right);
		}
		if(assignment.find(tree->name)==assignment.end())
		{
			cout <<tree->name <<"\t" << 0 << endl;;
		}
		else
		{
			cout <<tree->name <<"\t"<< assignment[tree->name] << endl;
		}
	}
	return 0;
}

int
Framework::getAssignFiltered(map<string,int>& assignment,SpeciesNode* tree, bool show,map<string,int>& genes)
{
	//This is the leaf node. We will "show" the leaf node only if the show variable is set to true.
	//show is set to true only if the ancestor of this node is set to true.
	if(tree->left==NULL && tree->right==NULL)
	{
		if(show)
		{
			if(assignment.find(tree->name)!=assignment.end() && assignment[tree->name]==1)
			{
				genes[tree->name]=0;
			}
		}
	}
	else
	{
		int assign=assignment[tree->name];
		bool showFlag=false;
		if(assign==1)
		{
			showFlag=true;
		}
		if(tree->left!=NULL)
		{
			getAssignFiltered(assignment,tree->left,showFlag,genes);
		}
		if(tree->right!=NULL)
		{
			getAssignFiltered(assignment,tree->right,showFlag,genes);
		}
	}
	return 0;
}


int
main(int argc, const char** argv)
{
	if(argc!=7)
	{
		cout <<"Usage: refineMotifs speciesorder orthogroups spectiestree motifcollection scorefile outputprefix" << endl;
		return 0;
	}
	Framework fw;
	fw.readSpeciesOrthology(argv[1],argv[2]);
	fw.readTree(argv[3]);
	fw.readMotifs(argv[4]);
	fw.setScore(argv[5]);
	fw.refineTargets(argv[6]);
	return 0;
}
