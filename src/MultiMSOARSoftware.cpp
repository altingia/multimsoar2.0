#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include "Hungarian.h"
#include "TreeCentric.h"
#include "NodeCentric.h"
#include "TreeAnalysis.h"

using namespace std;

// Species tree and the number of species
int S;
string speciesTree;
int maximumN;

// Vertex and Adjacency list
map<string,int> species;
map<string, vector<string> > adjacency;
map<pair<string,string>, double> edges;


map<string,int> visited;

vector<string> AllTrees;
vector<vector<string> > AllTreeGeneName;

set<string> AllGeneBirth;
set<string> AllGeneDuplication;
map<int, int> AllGeneLoss;

// Group connected by single linkage of ortholog pairs
vector<string> group;

// Single linkage of connected component
void DFS(string cur)
{
	group.push_back(cur);
	visited[cur]=1;
	for(int i=0; i<adjacency[cur].size(); i++)
	{
		string next=adjacency[cur][i];
		if(visited.count(next)==0)
			DFS(next);
	}
}

// Partition each group into N layers
void Partition()
{
	vector<vector<vector<string> > > v(S);

	for(int i=0; i<group.size(); i++)
	{
		int sp=species[group[i]];
		vector<string> tmp;
		tmp.push_back(group[i]);
		v[sp].push_back(tmp);
	}
	
	// N: the number of layers
	int N=0;
	for(int i=0; i<S; i++) N=max(N, (int)(v[i].size()));

	//cout<<"N="<<N<<endl;

	// Pad each species in a group with dummy vertices
	for(int i=0; i<S; i++)
	{
		vector<string> dummy;
		dummy.push_back("");

		for(int j=v[i].size(); j<N; j++) 
			v[i].push_back(dummy);
	}

	//////////////////////////////////
	
	vector<vector<vector<string> > > stack;

	int index=0;

	for(int i=0; i<speciesTree.size(); i++)
	{
		if(speciesTree[i]!='N')
		{
			stack.push_back(v[index++]);
		}
		else
		{
			vector<vector<string> > v2=stack.back(); stack.pop_back();
			vector<vector<string> > v1=stack.back(); stack.pop_back();

			vector<vector<int> > matrix(N, vector<int> (N) );

			// Calculate the added weight for an edge in Bipartite Graph
			for(int j=0; j<N; j++) for(int k=0; k<N; k++) matrix[j][k]=0;

			for(int j=0; j<N; j++) for(int k=0; k<N; k++)
			{
				for(int jj=0; jj<v1[j].size(); jj++) for(int kk=0; kk<v2[k].size(); kk++)
				{
					string gene1=v1[j][jj];
					string gene2=v2[k][kk];
					if(gene1=="" or gene2=="") matrix[j][k]+=0;
					if(edges.count(make_pair(gene1, gene2))) matrix[j][k]+=(int)edges[make_pair(gene1, gene2)];
				}
			}

			// Run the Hungarian maximum matching algorithm for weighted bipartite graph
			//cout<<"Running Hungarian ..."<<endl;
			Hungarian H(matrix);
			//cout<<"Hungarian Done."<<endl;

			// Merge the vertices after matching
			for(int j=0; j<N; j++)
			{
				int p=H.matchingX[j];
				for(int k=0; k<v2[p].size(); k++)
					v1[j].push_back(v2[p][k]);
			}

			stack.push_back(v1);
		}
	}
	
	// Cout the partition for each group
	//cout<<N<<" layers: "<<endl;
	vector<string> trees(N, string(S,'0'));
	vector<vector<string> > treeGeneName(N, vector<string> (S));

	for(int i=0; i<N; i++) for(int j=0; j<S; j++) treeGeneName[i][j]="";

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<stack[0][i].size(); j++)
		{
			if(stack[0][i][j]!="") 
			{
				trees[i][species[stack[0][i][j]]]='1';
				treeGeneName[i][species[stack[0][i][j]]]=stack[0][i][j];
			}
		}
		for(int j=0; j<speciesTree.size(); j++) if(speciesTree[j]=='N')
			trees[i]=trees[i].substr(0, j)+"N"+trees[i].substr(j);
		//cout<<trees[i]<<endl;
	}

	// Store the results;
	for(int i=0; i<N; i++) AllTrees.push_back(trees[i]);
	for(int i=0; i<N; i++) AllTreeGeneName.push_back(treeGeneName[i]);
}

void TreeLabeling(ofstream& orthoGroupOut)
{
	// Decide which labeling algorithm to use (NodeCentric or TreeCentric)
	vector<string> LabelResults;
	int totalSubstitutions;

	int N=AllTrees.size();
	maximumN=max(maximumN, N);

	if(N==0) return;

	if(N<5)
	{
		NodeCentric nc(AllTrees);
		totalSubstitutions=nc.totalSubstitutions;
		LabelResults=nc.optimalLabeling;
	}
	else
	{
		TreeCentric tc(AllTrees);
		totalSubstitutions=tc.totalSubstitutions;
		LabelResults=tc.optimalLabeling;
	}

	//cout<<"Optimal Labeling: "<<totalSubstitutions<<endl;
	//for(int i=0; i<N; i++)
	//	cout<<LabelResults[i]<<endl;

	// Tree Analysis
	//cout<<"Tree Analysis: "<<endl;

	TreeAnalysis ta(speciesTree, LabelResults, AllTreeGeneName);
	//ta.printAnalysis();
	//ta.printDetailedAnalysis();
	ta.printOrthoGroups(orthoGroupOut);
	ta.printGeneInfo();
	for(int i=0; i<ta.GeneBirth.size(); i++) AllGeneBirth.insert(ta.GeneBirth[i]);
	for(int i=0; i<ta.GeneDuplication.size(); i++) AllGeneDuplication.insert(ta.GeneDuplication[i]);
	for(int i=0; i<ta.GeneLoss.size(); i++) AllGeneLoss[ta.GeneLoss[i]]++;
	//cout<<"********************************************************"<<endl;
}

void printGeneInfo(char* filename)
{
	ofstream outfile(filename);
	outfile<<"Gene birth: ";
	for(set<string>::iterator it=AllGeneBirth.begin(); it!=AllGeneBirth.end(); it++)
		outfile<<*it<<"\t";
	outfile<<endl;

	outfile<<"Gene duplication: ";
	for(set<string>::iterator it=AllGeneDuplication.begin(); it!=AllGeneDuplication.end(); it++)
		outfile<<*it<<"\t";
	outfile<<endl;

	outfile<<"Gene loss: ";
	for(map<int,int>::iterator it=AllGeneLoss.begin(); it!=AllGeneLoss.end(); it++)
		outfile<<"Species"<<it->first<<"\t"<<it->second<<"\t";
	outfile<<endl;
	outfile.close();
}

int main(int argc, char** argv)
{
	if(argc!=6)
	{
		cout<<"Usage: MultiMSOAR2.0 <#species> <speciesTree> <GeneFamily> <-o GeneInfo> <-o OrthoGroups>"<<endl;
		exit(1);
	}

	ifstream infile(argv[2]);
	string tmpSpeciesTree;
	getline(infile, tmpSpeciesTree);
	infile.close();
	speciesTree="";
	string tmpL="";
	for(int i=0; i<tmpSpeciesTree.length(); i++)
	{
		if(tmpSpeciesTree[i]==',')
		{
			if(tmpL!="") { speciesTree+='1'; tmpL=""; }
		}
		else if(tmpSpeciesTree[i]==')')
		{
			if(tmpL!="")
			{
				speciesTree+='1';
			}
			speciesTree+='N';
			tmpL="";
			i++;
			while(i<tmpSpeciesTree.length() and tmpSpeciesTree[i]!=')' and tmpSpeciesTree[i]!=',')
				i++;
			i--;
		}
		else if(tmpSpeciesTree[i]==' ') continue;
		else
			tmpL+=tmpSpeciesTree[i];
	}
	

	S=atoi(argv[1]);
	//S=(speciesTree.size()+1)/2;

	int level=0, level2=1;
	for(int i=0; i<S; i++) for(int j=i+1; j<S; j++)
	{
		stringstream ss;
		ss<<"S"<<i<<"_S"<<j;
		ifstream file(ss.str().c_str());
		if( file.fail() )
		{
			cout<<"Cannot open file "<<ss.str()<<endl;
			exit(1);
		}

		string line;
		while(getline(file, line))
		{
			stringstream ss;
			ss<<line;
			string gene1, gene2;
			double score;
			ss>>gene1>>gene2>>score;

			species[gene1]=i;
			species[gene2]=j;

			adjacency[gene1].push_back(gene2);
			adjacency[gene2].push_back(gene1);

			edges[make_pair(gene1,gene2)]=score;
			edges[make_pair(gene2,gene1)]=score;
		}
	}
	
	map<int, set<string> > RealFamily;
	ifstream infile2(argv[3]);
	string oneFamily;
	int indexFamily=0;
	while(getline(infile2, oneFamily))
	{
		stringstream ss(oneFamily);
		string name;
		while(ss>>name) RealFamily[indexFamily].insert(name);
		indexFamily++;
	}
	infile2.close();


	// Generate groups using single linkage
	ofstream orthoGroupOut(argv[5]);
	for(map<int,set<string> >::iterator it=RealFamily.begin(); it!=RealFamily.end(); it++)
	{
		AllTrees.clear();
		AllTreeGeneName.clear();

		//if((it->second).size()<=1) continue;

		//cerr<<"Family "<<it->first<<"\t Layers N=";
		for(set<string>::iterator j=(it->second).begin(); j!=(it->second).end(); j++) if(visited.count(*j)==0)
		{
			group.clear();
			DFS(*j);
			//if(group.size()<=1) continue;
			Partition();
		}
		//cerr<<AllTrees.size()<<endl;
		TreeLabeling(orthoGroupOut);
	}
	orthoGroupOut.close();

	printGeneInfo(argv[4]);

	return 0;
}
