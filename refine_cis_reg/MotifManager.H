#ifndef _MOTIF_MANAGER_
#define _MOTIF_MANAGER_
#include <map>
#include <string>
using namespace std;
class MotifManager
{
	public:
	MotifManager();
	~MotifManager();
	int readMotifFile(const char*);
	map<string,map<string,int>*>& getAllMotifs();
	map<string,int>* getGenesForMotif(const char*);
	map<string,int>* getMotifsForGene(const char*);
	int addFilteredMotif(const char*, const char*);
	map<string,map<string,int>*>& getFilteredMotifs();
	private:
	map<string,map<string,int>*> geneMotifInstance;
	map<string,map<string,int>*> motifGeneInstance;
	map<string,map<string,int>*> filteredMotifInstance;
};
#endif
