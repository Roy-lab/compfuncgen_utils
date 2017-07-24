#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include "SpeciesNode.H"
#include "ParsimonyInferrer.H"

ParsimonyInferrer::ParsimonyInferrer()
{
}

ParsimonyInferrer::~ParsimonyInferrer()
{
}

int 
ParsimonyInferrer::setSpeciesTree(SpeciesNode* aroot)
{
	root=aroot;
	return 0;
}
	
int 
ParsimonyInferrer::setScore(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	map<int,int> letters;
	int linecnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		map<int,double>* valset=NULL;
		while(tok!=NULL)
		{
			if(linecnt==0)
			{
				int l=atoi(tok);
				letters[tokCnt]=l;
			}
			else 
			{
				if(tokCnt==0)
				{
					valset=new map<int,double>;
					int base=atoi(tok);
					scoreMatrix[base]=valset;				
				}	
				else
				{
					double v=atof(tok);
					(*valset)[tokCnt-1]=v;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		linecnt++;
	}
	inFile.close();
	return 0;
}

int 
ParsimonyInferrer::inferAncestralAssignment(map<string,int>& assignment, map<string,int>& fullassignment)
{
	GeneNode* genetree=populateChild(root);
	getParsimoniousAssignment(genetree,assignment);
	populateAssignment(genetree,fullassignment);
	freeMe(genetree);
	return 0;
}

ParsimonyInferrer::GeneNode*
ParsimonyInferrer::populateChild(SpeciesNode* speciesNode)
{
	ParsimonyInferrer::GeneNode* gnode=new ParsimonyInferrer::GeneNode;
	gnode->left=NULL;
	gnode->right=NULL;
	gnode->name.append(speciesNode->name);
	if(speciesNode->left!=NULL)
	{
		gnode->left=populateChild(speciesNode->left);
	}
	if(speciesNode->right!=NULL)
	{
		gnode->right=populateChild(speciesNode->right);
	}
	return gnode;
}

int
ParsimonyInferrer::freeMe(ParsimonyInferrer::GeneNode* gNode)
{
	if(gNode->left!=NULL)
	{
		freeMe(gNode->left);
	}
	if(gNode->right!=NULL)
	{
		freeMe(gNode->right);
	}
	gNode->values.clear();
	for(map<int,map<int,int>*>::iterator kIter=gNode->leftchildminimizer.begin();kIter!=gNode->leftchildminimizer.end();kIter++)
	{
		kIter->second->clear();
		delete kIter->second;
	}
	gNode->leftchildminimizer.clear();
	for(map<int,map<int,int>*>::iterator kIter=gNode->rightchildminimizer.begin();kIter!=gNode->rightchildminimizer.end();kIter++)
	{
		kIter->second->clear();
		delete kIter->second;
	}
	gNode->rightchildminimizer.clear();
	delete gNode;
	return 0;
}

int
ParsimonyInferrer::getParsimoniousAssignment(GeneNode* n,map<string,int>& assignment)
{
	if(n->left==NULL && n->right==NULL)
	{
		for(map<int,map<int,double>*>::iterator sIter=scoreMatrix.begin();sIter!=scoreMatrix.end();sIter++)
		{
			if(assignment[n->name]==sIter->first)
			{
				n->values[sIter->first]=0;
			}
			else
			{
				n->values[sIter->first]=1000;
			}
		}
	}
	else if(n->left!=NULL && n->right!=NULL)
	{
		getParsimoniousAssignment(n->left,assignment);
		getParsimoniousAssignment(n->right,assignment);
		for(map<int,map<int,double>*>::iterator sIter=scoreMatrix.begin();sIter!=scoreMatrix.end();sIter++)
		{
			double minval_left=1000;
			int minassign_left;
			double minval_right=1000;
			int minassign_right;
			map<int ,double>* scorevals=scoreMatrix[sIter->first];
			map<int,double> leftMinimizer;
			map<int,double> rightMinimizer;
			for(map<int,map<int,double>*>::iterator tIter=scoreMatrix.begin();tIter!=scoreMatrix.end();tIter++)
			{
				double subvalues=(*scorevals)[tIter->first];
				double childval_left=n->left->values[tIter->first];
				double myval=childval_left+subvalues;
				leftMinimizer[tIter->first]=myval;
				if(myval<minval_left)
				{
					minval_left=myval;
					minassign_left=tIter->first;
				}
				
				double childval_right=n->right->values[tIter->first];
				myval=childval_right+subvalues;
				if(myval<minval_right)
				{
					minval_right=myval;
					minassign_right=tIter->first;
				}
				rightMinimizer[tIter->first]=myval;
			}
			n->values[sIter->first]=minval_right+minval_left;
			map<int,int>* childState_left=new map<int,int>;
			n->leftchildminimizer[sIter->first]=childState_left;
			map<int,int>*childState_right=new map<int,int>;
			n->rightchildminimizer[sIter->first]=childState_right;
			double mval=leftMinimizer[minassign_left];
			for(map<int,double>::iterator aIter=leftMinimizer.begin();aIter!=leftMinimizer.end();aIter++)
			{
				if(aIter->second==mval)
				{
					(*childState_left)[aIter->first]=mval;
				}
			}
			mval=rightMinimizer[minassign_right];
			for(map<int,double>::iterator aIter=rightMinimizer.begin();aIter!=rightMinimizer.end();aIter++)
			{
				if(aIter->second==mval)
				{
					(*childState_right)[aIter->first]=mval;
				}
			}
		}
	}
	else
	{
		cout <<"Tree is messed up" << endl;
		exit(0);
	}
	return 0;
}

int
ParsimonyInferrer::populateAssignment(ParsimonyInferrer::GeneNode* genenode, map<string,int>& assignment)
{
	double minResidue=1000;
	int minVal;
	for(map<int,double>::iterator sIter=genenode->values.begin();sIter!=genenode->values.end();sIter++)
	{
		//We will disambiguate based on the loss. That is if both absence or presence are giving us
		//equal valuess we will take the one where the motif is lost. The loss in turn is the smallest value
		//so minResidue is updated if there is another value that has a lower values than the "loss" character
		if(sIter->second<=minResidue)
		{
			minResidue=sIter->second;
			minVal=sIter->first;
		}
	}
	cout <<"Score of tree " << minResidue << endl;
	int tie=0;
	for(map<int,double>::iterator sIter=genenode->values.begin();sIter!=genenode->values.end();sIter++)
	{
		if(sIter->second!=minResidue)
		{
			continue;
		}
		tie++;
	}
	if(tie>1)
	{
		cout <<"There is a tie for " << minVal<< endl;
	}
	getOptimalAssign(genenode,minVal,assignment);
	return 0;
}

int
ParsimonyInferrer::getOptimalAssign(ParsimonyInferrer::GeneNode* node,int residue, map<string,int>& value)
{
	value[node->name]=residue;
	if(node->left!=NULL)
	{
		int leftchildresidue=node->leftchildminimizer[residue]->begin()->first;
		getOptimalAssign(node->left,leftchildresidue,value);
	}
	if(node->right!=NULL)
	{
		int rightchildresidue=node->rightchildminimizer[residue]->begin()->first;
		getOptimalAssign(node->right,rightchildresidue,value);
	}
	return 0;
}
