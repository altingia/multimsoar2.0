#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>

using namespace std;

class Tree
{
	public:
		map<long long,int> AccCost;
		map<long long,int> CurLabel;
		map<long long,long long> PreValue;
		Tree(){}
};


class TreeCentric
{
	public:

// Store the labeling results
int totalSubstitutions;
vector<string> optimalLabeling;

// N: the number of trees
// S: the number of internal nodes
int N, S;

// All the trees
vector<string> trees;

// Internal vector to store the labeling of the current tree
vector<int> w;


// Given a tree, the valid labeling of the tree and its corresponding
// cost is stored in validLabeling
map<long long,int> validLabeling;

void go2(string tree, long long label, int cost)
{
	if(tree.size()==1)
	{
		if(validLabeling.count(label)==0 or validLabeling[label]>cost)
			validLabeling[label]=cost;
		return;
	}
	int N_pos=tree.find('N');
	int left=tree[N_pos-2]-'0';
	int right=tree[N_pos-1]-'0';
	string prefix=tree.substr(0, N_pos-2);
	string suffix=tree.substr(N_pos+1);

	if(left==0 and right==0)
	{
		go2(prefix+string(1,'0')+suffix, (label<<1)|0, cost);
	}
	else if(left==1 and right==1)
	{
		go2(prefix+string(1,'1')+suffix, (label<<1)|1, cost);
	}
	else if((left==1 and right==0) or (left==0 and right==1))
	{
		go2(prefix+string(1,'1')+suffix, (label<<1)|1, cost+1);
		go2(prefix+string(1,'2')+suffix, (label<<1)|0, cost+1);
	}	
	else 
	{
		go2(prefix+string(1,'2')+suffix, (label<<1)|0, cost+((left|right)&1));
	}
}

void Valid_Internal_Labeling(string tree)
{
	validLabeling.clear();
	go2(tree,0,0);
	// cout<<"Valid Labeling: "<<validLabeling.size()<<endl;
}

// Print the bit value of an integer
string printValue(long long p)
{
	string s="";
	for(int i=0; i<S; i++)
	{
		s=string(1,'0'+(p&1))+s;
		p>>=1;
	}
	return s;
}

// Check 1-0*-1 constraint between trees (InterTree 10*1 constraint)
bool One_Oh_One_Constraint(long long p)
{
	vector<long long> vec;
	for(int i=0; i<S; i++)
	{
		vec.push_back(p&1LL);
		p>>=1LL;
	}
	reverse(vec.begin(), vec.end());

	vector<int> stack;
	int index=0;
	for(int i=0; i<trees[0].size(); i++)
	{
		if(trees[0][i]=='0' or trees[0][i]=='1')
		{
			int tmp=0;
			for(int j=0; j<N; j++)
				tmp+=trees[j][i]-'0';

			if(tmp>0) stack.push_back(1);
			else stack.push_back(tmp);
		}
		else
		{
			int right=stack.back();
			stack.pop_back();
			int left=stack.back();
			stack.pop_back();
			long long parent=vec[index++];
			if(parent==1)
			{
				if(left==2 or right==2) return false;
				else stack.push_back(1);
			}
			else
			{
				if(left==0 and right==0) stack.push_back(0);
				else stack.push_back(2);
			}
		}
	}
	return true;
}


// Check the 0-1 Constraint (if a node is 0, then at least one of its substree are all 0s)
bool Zero_One_Constraint(long long p)
{
	vector<int> vec;
	for(int i=0; i<S; i++)
	{
		vec.push_back(p&1);
		p>>=1;
	}
	reverse(vec.begin(), vec.end());

	vector<int> stack;
	int index=0;
	for(int i=0; i<trees[0].size(); i++)
	{
		if(trees[0][i]=='0' or trees[0][i]=='1')
		{
			int tmp=0;
			for(int j=0; j<N; j++)
				tmp+=trees[j][i]-'0';

			if(tmp>0) stack.push_back(1);
			else stack.push_back(tmp);
		}
		else
		{
			int right=stack.back();
			stack.pop_back();
			int left=stack.back();
			stack.pop_back();
			int parent=vec[index++];
			if(parent==0)
			{
				if(left>0 and right>0) return false;
				else if(left>0 or right>0) parent=1;
			}

			stack.push_back(parent);
		}
	}
	return true;
}

// Update the current tree
void UpdateCurrentTree(string tree, Tree* cur, Tree* pre)
{
	Valid_Internal_Labeling(tree);

	for(map<long long,int>::iterator i=(pre->AccCost).begin(); i!=(pre->AccCost).end(); i++)
		for(map<long long,int>::iterator j=validLabeling.begin(); j!=validLabeling.end(); j++)
		{
			long long preV=i->first;
			int preCost=i->second;
			long long curLabel=j->first;
			int curCost=j->second;

			long long curV=0;
			for(int k=0; k<S; k++)
			{
				int p=preV&1;
				int q=curLabel&1;
				curV=curV|(((long long)(p|q))<<k);
				preV>>=1;
				curLabel>>=1;
			}

			if((cur->AccCost).count(curV)==0 or (cur->AccCost[curV])>preCost+curCost)
			{
				(cur->AccCost)[curV]=preCost+curCost;
				(cur->CurLabel)[curV]=j->first;
				(cur->PreValue)[curV]=i->first;
			}
		}
}


TreeCentric(vector<string> input)
{
	trees=input;

	N=trees.size();
	S=trees[0].size()/2;

	vector<Tree*>  v(N+1);
	for(int i=0; i<N+1; i++) v[i]=new Tree();

	// Initialization
	(v[0]->AccCost)[0]=0;
	(v[0]->CurLabel)[0]=-1;


	for(int i=0; i<N; i++)
	{
		UpdateCurrentTree(trees[i], v[i+1], v[i]);
		//cout<<"Accumative Size: "<<(v[i+1]->AccCost).size()<<endl;
	}

	// Find the final optimal solution
	int totalSub=1<<30;
	long long finalV=-1;
	for(map<long long,int>::iterator i=(v[N]->AccCost).begin(); i!=(v[N]->AccCost).end(); i++) 
		if(Zero_One_Constraint(i->first) and One_Oh_One_Constraint(i->first))
		{
			if((i->second) < totalSub)
			{
				totalSub=i->second;
				finalV=i->first;
			}
		}
	// cout<<"\t"<<totalSub<<endl;

	optimalLabeling.resize(N);
	optimalLabeling[N-1]=printValue((v[N]->CurLabel)[finalV]);
	totalSubstitutions=totalSub;

	// Trace back the labeling of each tree backward
	for(int i=N-1; i>0; i--)
	{
		finalV=(v[i+1]->PreValue)[finalV];
		optimalLabeling[i-1]=printValue((v[i]->CurLabel)[finalV]);
		//cout<<"\t"<<(v[i]->AccCost)[finalV];
		//cout<<endl;
	}
	for(int i=0; i<N; i++)
	{
		int index=0;
		string stmp="";
		for(int j=0; j<trees[i].size(); j++)
			if(trees[i][j]!='N') stmp+=string(1,trees[i][j]);
			else stmp+=string(1, optimalLabeling[i][index++]);

		optimalLabeling[i]=stmp;
	}
}

};

