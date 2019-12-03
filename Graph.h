#ifndef GRAPH_H_
#define GRAPH_H_
#include "Motif.h"
#include <vector>
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
using namespace std;

class Graph{
public:
	int NodeSize;
	int EdgeSize;
	int LayerSize;
	//aggregate graph
    vector<dynamic_bitset<> > adjacencyList;
	int **adjacencyMatrix;
	int k;
	vector<dynamic_bitset<> >edgeLayers; //for each edge, store the set of layers that have this edge
	vector<pair<string, string>> edgeset; //for each edge, record original index
	vector<pair<int, int>> edgeIntSet; //for each edge, record the mapping id

public:
	Graph(int layerNum, int nodeNum, string networkName, double thre, unordered_map<string, int>& nodeMap): LayerSize(layerNum), NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
		 adjacencyMatrix = new int*[nodeNum];
		 for(int i=0;i<nodeNum;i++){
			 adjacencyMatrix[i] = new int[nodeNum]();
		 }
         k = ceil(layerNum*thre);
		 readNetworkFile(networkName, nodeMap);
	}

  ~Graph(){
		for(int i=0;i<NodeSize;i++) {
			 delete[] adjacencyMatrix[i];
		}
		delete[] adjacencyMatrix;
		adjacencyList.clear();
	}

  void findMotifs(vector<Motif*>& pattern, int caseN);
private:
  void readNetworkFile(string networkName, unordered_map<string, int>& nodeMap);
};

#endif
