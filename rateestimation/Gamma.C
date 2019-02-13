
#include <string.h>
#include "Matrix.H"
#include "SpeciesDistManager.H"
#include "Gamma.H"

Gamma::Gamma()
{
	dupAncName.clear();	
	nodeID=0;
}

Gamma::~Gamma()
{
}


int 
Gamma::setGeneClusterID(string& geneName,string& specName,int clusterID)
{
	if(geneMap.find(geneName)==geneMap.end())
	{
		cout <<"No gene " << geneName <<  " in species " << specName << endl;
		exit(0);
	}
	Gamma::Node* gene=geneMap[geneName];
	Matrix* gamma=gene->gamma;
	gamma->setValue(1.0,0,clusterID);
	return 0;
}

//Here we will assume that there is a single duplication even at most. In the future, we want
//to initialize this directly from the gene tree.
int 
Gamma::initUsingTree(SpeciesDistManager::Species* specnode)
{
	initSubtree(specnode,specnode->name);
	return 0;
}


Gamma::Node*
Gamma::initSubtree(SpeciesDistManager::Species* specnode,string& rootname)
{
	Gamma::Node* node=new Gamma::Node;
	node->parent=NULL;
	node->name.append(specnode->name);
	if(strcmp(specnode->name.c_str(),rootname.c_str())==0)
	{
		root=node;
		node->gamma=new Matrix(1,maxClusterCnt);
	}
	else
	{
		node->gamma=new Matrix(maxClusterCnt,maxClusterCnt);
	}
	node->gamma->setAllValues(0);
	node->normTerm=NULL;
	node->alpha=NULL;
	node->beta=NULL;
	geneMap[node->name]=node;
	if(specnode->leftchild!=NULL)
	{
		Gamma::Node* leftchild=initSubtree(specnode->leftchild,rootname);
		node->leftchild=leftchild;
		if(strcmp(leftchild->species.c_str(),rootname.c_str())!=0)
		{
			leftchild->parent=node;
		}
	}
	else
	{
		node->leftchild=NULL;
	}
	if(specnode->rightchild!=NULL)
	{
		Gamma::Node* rightchild=initSubtree(specnode->rightchild,rootname);
		node->rightchild=rightchild;
		if(strcmp(rightchild->species.c_str(),rootname.c_str())!=0)
		{
			rightchild->parent=node;
		}
	}
	else
	{
		node->rightchild=NULL;
	}
	return node;
}

Matrix* 
Gamma::getGamma(string& speciesName)
{
	if(geneMap.find(speciesName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[speciesName];
	return node->gamma;
}

Matrix* 
Gamma::getNormTerm(string& speciesName)
{
	if(geneMap.find(speciesName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[speciesName];
	return node->normTerm;
}


Gamma::Node* 
Gamma::getGeneNode(string& speciesName)
{
	if(geneMap.find(speciesName)==geneMap.end())
	{
		return NULL;
	}
	Gamma::Node* node=geneMap[speciesName];
	return node;
}

Gamma::Node* 
Gamma::getRoot()
{
	return root;
}


int 
Gamma::showTree()
{
	showTree(root);
	return 0;
}



string&
Gamma::getDupAncestor()
{
	return dupAncName;
}

int
Gamma::showTree(Gamma::Node* node)
{
	if(node->leftchild==NULL && node->rightchild==NULL)
	{
		return 0;
	}
	if(node->leftchild!=NULL)
	{
		cout << node->name <<"->left "<< node->leftchild->name << endl;
	}
	else
	{
		cout << node->name <<"->left null"<< endl;
	}
	if(node->rightchild!=NULL)
	{
		cout << node->name <<"->right "<< node->rightchild->name << endl;
	}
	else
	{
		cout << node->name <<"->right null"<< endl;
	}
	if(node->leftchild!=NULL)
	{
		showTree(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		showTree(node->rightchild);
	}

	return 0;
}
