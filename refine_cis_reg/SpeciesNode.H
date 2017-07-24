#ifndef _SPECIESNODE_
#define _SPECIESNODE_
#include <string>
using namespace std;
class SpeciesNode
{
	public:
	SpeciesNode()
	{
	}
	~SpeciesNode()
	{
	}
	SpeciesNode* left;
	SpeciesNode* right;
	SpeciesNode* parent;
	string name;
};
#endif
