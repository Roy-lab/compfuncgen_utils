#ifndef _GAMMA_
#define _GAMMA_
#include <map>
#include <string>
using namespace std;
class Matrix;
class SpeciesDistManager;
class Gamma
{
	public:
	Gamma();
	~Gamma();
	inline int setMaxClusterCnt(int k)
	{
		maxClusterCnt=k;
	}
	int addGene(string& geneName,string& specName,int maxClusterCnt);
	int setGeneClusterID(string&,string&,int);
	int initUsingTree(SpeciesDistManager::Species*);
	Matrix* getGamma(string&);
	Matrix* getNormTerm(string&);
	double getScore();
	struct Node;
	struct Node
	{
		string name;
		string species;
		Node* leftchild;
		Node* rightchild;
		Node* parent;
		//This is 2 for duplication and 1 for speciation
		int nodeType;
		Matrix* normTerm;
		Matrix* gamma;
		Matrix* alpha;
		Matrix* beta;
	};
	Node* getRoot();
	Node* getGeneNode(string&);
	Node* root;
	int showTree();
	string& getDupAncestor();
	private:
	Gamma::Node* initSubtree(SpeciesDistManager::Species* node,string& rootname);
	int showTree(Node*);
	map<string,Node*> geneMap;
	map<string,Matrix*> gammas;
	int maxClusterCnt;
	int nodeID;
	string dupAncName;
	
};
#endif
