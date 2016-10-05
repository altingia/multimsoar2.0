#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>

#define MAXINT (1<<30)

using namespace std;

class Node
{
	public:
		map<long long,int> changes;
		map<long long,long long> leftV, rightV;
		Node* left;
		Node* right;
		Node(){left=NULL; right=NULL;}
		Node(long long v) { changes[v]=0; left=right=NULL; }
};


class NodeCentric
{
	public:

//Store the results;
int totalSubstitutions;
vector<string> optimalLabeling;

//N: the number of trees
int N;
long long checkZero;


//The input trees
vector<string> trees;

//Labeling Results
vector<long long> results;

map<long long, int> leftMap, rightMap;
map<long long, long long> leftMapNode, rightMapNode;
//map<long long, int> leftZero, rightZero;


vector<pair<int,int> > vp[50];

//Print the value
string printValue(long long n)
{
	string s="";
	for(int i=0; i<N; i++)
	{
		int p=n&3;
		s=string(1, '0'+p)+s;
		n>>=2LL;
	}
	return s;
}

bool moreThanTwoOnes(long long p)
{
	int cnt=0;
	for(int i=0; i<N; i++)
	{
		if(p&1) cnt++;
		p>>=2;
	}
	return cnt>1;
}


void PostOrderTraversal(Node* cur, long long value)
{
	//if(cur->left==NULL or cur->right==NULL) return;
	
	if(cur==NULL) return;
	PostOrderTraversal(cur->left, (cur->leftV)[value]);
	PostOrderTraversal(cur->right,(cur->rightV)[value]);

	results.push_back(value);
}

int getBit(long long p, int i)
{
	for(int j=0; j<i; j++) { p>>=2; }
	return p&3;
}

void go(Node* left, Node* right, Node* child)
{
	long long lv=((left->changes).begin())->first;
	long long rv=((right->changes).begin())->first;
	long long cv=((child->changes).begin())->first;

	for(int i=0; i<N; i++)
	{
		vp[i].clear();

		int lb=lv&3; lv>>=2;
		int rb=rv&3; rv>>=2;
		int cb=cv&3; cv>>=2;

		if( lb==0 and rb==0 )
		{
			vp[i].push_back(make_pair(0,0));
		}
		else
		{
			if(cb==0)
			{
				vp[i].push_back(make_pair(1,0));
				vp[i].push_back(make_pair(2,0));
			}
			else
			{
				vp[i].push_back(make_pair(1,1));
				vp[i].push_back(make_pair(2,1));
				vp[i].push_back(make_pair(2,2));
			}
		}
	}
}

map<long long, int> leftZero, rightZero;
map<long long, int> leftNonZeroMin, rightNonZeroMin;
map<long long, long long> leftNonZeroNode, rightNonZeroNode;


void allCombination(int pos, long long parentV, long long childV, int sub, Node* child)
{
	if(pos==N)
	{
		if((child->changes).count(childV))
		{
			if((parentV & checkZero)==0)
			{
				if(childV==0) 
					leftZero[parentV]=sub+(child->changes)[childV];
				else 
				{
					if(leftNonZeroMin.count(parentV)==0 or leftNonZeroMin[parentV] > sub+(child->changes)[childV])
					{
						leftNonZeroMin[parentV]=sub+(child->changes)[childV];
						leftNonZeroNode[parentV]=childV;
					}
				}
			}
			else
			{
				if( (childV & checkZero)==0 and childV!=0 ) return;

				if(leftMap.count(parentV)==0 or leftMap[parentV]>sub+(child->changes)[childV])
				{
					leftMap[parentV]=sub+(child->changes)[childV];
					leftMapNode[parentV]=childV;
				}
			}
		}
		else
			return;
	}
	else
	{
		for(int i=0; i<vp[pos].size(); i++)
		{
			long long p=vp[pos][i].first;
			long long q=vp[pos][i].second;

			allCombination(pos+1, parentV|(p<<(2*pos)), childV|(q<<(2*pos)), sub+((p&1)^(q&1)), child);
		}
	}
}

void allCombination2(int pos, long long parentV, long long childV, int sub, Node* child)
{
	if(pos==N)
	{
		if((child->changes).count(childV))
		{
			if((parentV & checkZero)==0)
			{
				if(childV==0)
					rightZero[parentV]=sub+(child->changes)[childV];
				else
				{
					if(rightNonZeroMin.count(parentV)==0 or rightNonZeroMin[parentV] > sub+(child->changes)[childV])
					{
						rightNonZeroMin[parentV] = sub+(child->changes)[childV];
						rightNonZeroNode[parentV] = childV;
					}
				}
			}
			else
			{
				if( (childV & checkZero)==0 and childV!=0 ) return;

				if(rightMap.count(parentV)==0 or rightMap[parentV]>sub+(child->changes)[childV])
				{
					rightMap[parentV]=sub+(child->changes)[childV];
					rightMapNode[parentV]=childV;
				}
			}
		}
		else
			return;
	}
	else
	{
		for(int i=0; i<vp[pos].size(); i++)
		{
			long long p=vp[pos][i].first;
			long long q=vp[pos][i].second;

			allCombination2(pos+1, parentV|(p<<(2*pos)), childV|(q<<(2*pos)), sub+((p&1)^(q&1)), child);
		}
	}
}

NodeCentric(vector<string> input)
{

	trees=input;

	N=trees.size();
	checkZero=0;
	for(int i=0; i<N; i++) checkZero=(checkZero<<2)|1LL;

	vector<Node*> stack;

	for(int i=0; i<trees[0].size(); i++)
	{
		if(trees[0][i]=='0' or trees[0][i]=='1')
		{
			long long tmp=0;
			for(int j=0; j<N; j++) tmp=(tmp<<2)|(trees[j][i]-'0');
			Node* newNode=new Node(tmp);
			stack.push_back(newNode);
		}
		else
		{
			Node* right=stack.back();
			stack.pop_back();
			Node* left=stack.back();
			stack.pop_back();

			Node* newNode=new Node();
			newNode->left=left;
			newNode->right=right;

			///////////////////
			leftMap.clear(); rightMap.clear();
			leftMapNode.clear(); rightMapNode.clear();
			leftZero.clear(); rightZero.clear();
			leftNonZeroMin.clear(); rightNonZeroMin.clear();
			leftNonZeroNode.clear(); rightNonZeroNode.clear();
			//leftZero=-1; rightZero=-1;
			//leftNonZeroMin=MAXINT; rightNonZeroMin=MAXINT;


			long long start=time(NULL);
			long long end;
			go(left, right, left);
//			cout<<"Go Time: "<<time(NULL)-start<<endl;
			allCombination(0, 0, 0, 0, left);
			end=time(NULL);
//			cout<<"Combination 1 done: "<<end-start<<endl;


			start=time(NULL);
			go(left, right, right);
//			cout<<"Go Time: "<<time(NULL)-start<<endl;
			allCombination2(0, 0, 0, 0, right);
			end=time(NULL);
//			cout<<"Combination 2 done: "<<end-start<<endl;


			for(map<long long, int>::iterator it=leftMap.begin(); it!=leftMap.end(); it++)
			{
				long long value=it->first;
				if(rightMap.count(value)==0) continue;
				(newNode->changes)[value]=leftMap[value]+rightMap[value];
				(newNode->leftV)[value]=leftMapNode[value];
				(newNode->rightV)[value]=rightMapNode[value];
			}
//			cout<<"Non zero done"<<endl;

			for(map<long long, int>::iterator it=leftZero.begin(); it!=leftZero.end(); it++)
			{
				long long value=it->first;
				if(rightZero.count(value) > 0)
				{
					if( (newNode->changes).count(value)==0 or (newNode->changes)[value]>leftZero[value]+rightZero[value] )
					{
						(newNode->changes)[value]=leftZero[value]+rightZero[value];
						(newNode->leftV)[value]=0LL;
						(newNode->rightV)[value]=0LL;
					}
				}
				if(rightNonZeroMin.count(value) > 0)
				{
					if( (newNode->changes).count(value)==0 or (newNode->changes)[value]>leftZero[value]+rightNonZeroMin[value] )
					{
						(newNode->changes)[value]=leftZero[value]+rightNonZeroMin[value];
						(newNode->leftV)[value]=0LL;
						(newNode->rightV)[value]=rightNonZeroNode[value];
					}
				}
			}
	
			for(map<long long, int>::iterator it=rightZero.begin(); it!=rightZero.end(); it++)
			{
				long long value=it->first;
				if(leftZero.count(value) > 0)
				{
					if( (newNode->changes).count(value)==0 or (newNode->changes)[value]>leftZero[value]+rightZero[value] )
					{
						(newNode->changes)[value]=leftZero[value]+rightZero[value];
						(newNode->leftV)[value]=0LL;
						(newNode->rightV)[value]=0LL;
					}
				}
				if(leftNonZeroMin.count(value) > 0)
				{
					if( (newNode->changes).count(value)==0 or (newNode->changes)[value]>leftNonZeroMin[value]+rightZero[value] )
					{
						(newNode->changes)[value]=leftNonZeroMin[value]+rightZero[value];
						(newNode->leftV)[value]=leftNonZeroNode[value];
						(newNode->rightV)[value]=0LL;
					}
				}
			}
//			cout<<"Zero done"<<endl;

			///////////////////

//			cout<<"Current Size: "<<(newNode->changes).size()<<endl;
			stack.push_back(newNode);
		}
	}

	//Find the optimal value
	int minSub=MAXINT;
	long long label=-1;

	for(map<long long, int>::iterator i=(stack[0]->changes).begin(); i!=(stack[0]->changes).end(); i++)
	{
		if(i->second < minSub)
		{
			minSub=i->second;
			label=i->first;
		}
	}
//	cout<<minSub;
//	cout<<"Minimum Substitutions: "<<minSub<<endl;
	totalSubstitutions=minSub;

	//Post-order Traversal
	PostOrderTraversal(stack[0], label);


	//Print labeling of each tree
	for(int i=N-1; i>=0; i--)
	{
		stringstream sss;
		for(int j=0; j<results.size(); j++)
			sss<< ((results[j]>>(2*i))&1);
		optimalLabeling.push_back(sss.str());
	}

}

};
