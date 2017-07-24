#ifndef _PARSIMONYINFERRER
#define _PARSIMONYINFERRER
#include <map>
#include <string>
using namespace std;
class SpeciesNode;
class ParsimonyInferrer
{
	public:
	ParsimonyInferrer();
	~ParsimonyInferrer();
	int setSpeciesTree(SpeciesNode*);
	int setScore(const char*);
	int inferAncestralAssignment(map<string,int>&, map<string,int>&);
	struct GeneNode
	{
		GeneNode* left;
		GeneNode* right;
		GeneNode* parent;
		string name;
		map<int,double> values;
		map<int,map<int,int>*> leftchildminimizer;
		map<int,map<int,int>*> rightchildminimizer;
	};
	ParsimonyInferrer::GeneNode* populateChild(SpeciesNode* speciesNode);
	private:
	int getParsimoniousAssignment(GeneNode* n,map<string,int>& assignment);
	int populateAssignment(ParsimonyInferrer::GeneNode* genenode, map<string,int>& assignment);
	int getOptimalAssign(ParsimonyInferrer::GeneNode* node,int residue, map<string,int>& value);
	map<int,map<int,double>*> scoreMatrix;
	SpeciesNode* root;
	int freeMe(ParsimonyInferrer::GeneNode*);
};
#endif
