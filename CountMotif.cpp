#include "Graph.h"
#include <iostream>
#include <chrono>
#include <map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <random>
using namespace std;

Graph* g;


//selectNode finds all nodes that are in the selected embeddings
int F2count(vector<Motif*>& pattern, unordered_set<string>& selectNode){
	int countF2 = 0;
	int size = pattern.size();

	//build overlap graph
	vector<dynamic_bitset<> > madjacencyList(size, dynamic_bitset<>(size));
	vector<dynamic_bitset<> > layerAdj(g->LayerSize, dynamic_bitset<>(size));
	dynamic_bitset<>::size_type it, itt;

	for(int i=0;i<size;i++){
		Motif* a=pattern[i];

        for (int j = 0; j < g->LayerSize; j++) {
			if (a->layers[j] == 1) {
				layerAdj[j][i] = 1;
			}
		}

		for(int j=i+1;j<size;j++){
			Motif* b = pattern[j];
			dynamic_bitset<> tmp = b->edges;
			tmp &= a->edges;
			if (tmp.any()) {
				madjacencyList[i][j] = 1;
				madjacencyList[j][i] = 1;
			}
		}
	}

	dynamic_bitset<> neighbors(size);
	dynamic_bitset<> networks(g->LayerSize);
	networks.set();
	dynamic_bitset<> embeddings(size);
	embeddings.set();
	vector<int> candidates;
    vector<string> strs;

    cout<<"--------------------Select motifs----------------"<<endl;

	while (embeddings.any()) {
		double curMinLoss = embeddings.count()+1;
	    //calculate embeddings's loss
	    it=embeddings.find_first();
	    while(it!=dynamic_bitset<>::npos){
	    	Motif* b=pattern[it];
	    	dynamic_bitset<> curNetworks = b->layers;

	    	neighbors = madjacencyList[it];
		    double loss = neighbors.count();
		    neighbors[it] = 1;
		    dynamic_bitset<> rmLayers = curNetworks;
		    rmLayers ^= networks;

		    if (rmLayers.any() && loss <= curMinLoss) {
		    	dynamic_bitset<> others(size);
			    itt = rmLayers.find_first();
			    while(itt != dynamic_bitset<>::npos) {
				    others |= layerAdj[itt];
				    itt = rmLayers.find_next(itt);
			    }

			    dynamic_bitset<> tmp = neighbors;
			    tmp &= others;
			    others ^= tmp;

			    itt = others.find_first();
			    while(itt!=dynamic_bitset<>::npos && loss <= curMinLoss){
				    Motif* c = pattern[itt];
				    dynamic_bitset<> tmp = c->layers;
				    tmp &= curNetworks;
				    if (tmp.count() < g->k) {
					    neighbors[itt] = 1;
					    loss += 1;
				    }
				    else {
					    loss += 0.1*(1 - tmp.count()*1.0/c->layers.count());
				    }
				    itt=others.find_next(itt);
			    }
		    }

            if (loss <= curMinLoss) {
			    if (curMinLoss - loss < 0.0000001) {
				    candidates.push_back(it);
			    }
			    else {
				    candidates.clear();
				    candidates.push_back(it);
				    curMinLoss = loss;
			    }
		    }
		    b->removes = neighbors;
		    it=embeddings.find_next(it);
	    }

	    //select the candidates with least loss
	    Motif* a;
	    if (candidates.size() == 1) {
	    	a = pattern[candidates[0]];
	    }
	    else {
	    	std::random_device r;
		    std::default_random_engine generator(r());
		    std::uniform_int_distribution<int> distribution(0,candidates.size()-1);
		    int index = distribution(generator);
		    a = pattern[candidates[index]];
	    }
	    candidates.clear();

	    it = a->edges.find_first();
	    while(it != dynamic_bitset<>::npos) {
	    	cout<<"("<<g->edgeset[it].first<<", "<<g->edgeset[it].second<<") ";
		    selectNode.insert(g->edgeset[it].first);
		    selectNode.insert(g->edgeset[it].second);
		    it = a->edges.find_next(it);
	    }
	    cout<<endl;

	    //cout<<"select "<<a->id<<" with loss "<<curMinLoss<<endl;
	    countF2++;
	    networks &= a->layers;

	    //delete a and its neighbors
	    neighbors = a->removes;
	    for (int i = 0; i < g->LayerSize; i++) {
	    	dynamic_bitset<> tmp = layerAdj[i];
		    tmp &= neighbors;
		    layerAdj[i] ^= tmp;
	    }

	    it = embeddings.find_first();
	    while(it!=dynamic_bitset<>::npos){
	    	Motif* b = pattern[it];
		    if (neighbors[it] == 1) {
			    embeddings[it] = 0;
			    delete pattern[it];
		    }
		    else {
			    b->layers &= networks;
			    dynamic_bitset<> tmp = madjacencyList[it];
			    tmp &= neighbors;
			    madjacencyList[it] ^= tmp;
		    }
		    it=embeddings.find_next(it);
	    }
    }

	neighbors.clear();
	pattern.clear();
	madjacencyList.clear();
	return countF2;
}

int main(int argc, char *argv[]){

	//parameters: which motif, numberOfLayers, threshold, network name (file path)
	int type = atoi(argv[1]);
	int numberLayers = atoi(argv[2]);
	double threshold = atof(argv[3]);
	string networkName = argv[4];

	//first read the original network, map node name to id (int)
	ifstream topology;
	topology.open(networkName);
	string str;
	vector<string> strs;
    unordered_map<string, int> nodeMap; //input string vs map id
	unordered_map<int, string> nameMap; //map id vs input string
	int index = 0;

	while (!topology.eof()) {
		getline(topology, str);
		if (str.size() > 0 && str[0] != '-') {
			strs.clear();
			boost::split(strs,str,boost::is_any_of("\t"));
			string s1 = strs[0];
			string s2 = strs[1];

			if (nodeMap.find(s1) == nodeMap.end()) {
				nameMap.insert(make_pair(index, s1));
				nodeMap.insert(make_pair(s1, index++));

			}

			if (nodeMap.find(s2) == nodeMap.end()) {
				nameMap.insert(make_pair(index, s2));
				nodeMap.insert(make_pair(s2, index++));
			}
		}
	}
	topology.close();

	int numberNodes = index;
	auto f2_start = chrono::steady_clock::now();
	g = new Graph(numberLayers, numberNodes, networkName, threshold, nodeMap);
	vector<Motif*> pattern;
	g->findMotifs(pattern,type);
	unordered_set<string> motifNodeSet;
	int numEmbeddings = F2count(pattern, motifNodeSet);
	delete g;
	auto f2_end = chrono::steady_clock::now();
	cout<<numEmbeddings<<"\t"<<chrono::duration_cast<chrono::milliseconds>(f2_end - f2_start).count()<<"ms"<<endl;

}
