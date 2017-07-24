#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
using namespace std;

class MappedOrthogroupReader;
class ParsimonyInferrer;
class SpeciesNode;
class Framework
{
	public:
	Framework();
	~Framework();
	int readMotifs(const char*);
	int readTree(const char*);
	int readSpeciesOrthology(const char*, const char*);
	int setScore(const char*);
	int refineTargets(const char*);
	private:
	int showAssign(map<string,int>&, SpeciesNode* s);
	int getAssignFiltered(map<string,int>&, SpeciesNode* s,bool,map<string,int>&);
	MappedOrthogroupReader mor;
	SpeciesNode* root;
	ParsimonyInferrer pi;
	map<string,MotifManager*> motifCollection;
	map<string,SpeciesNode*> nodeSet;
	map<string,int> motifNames;
};
#endif
