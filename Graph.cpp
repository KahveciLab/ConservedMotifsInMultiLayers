#include "Graph.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <unordered_set>

void Graph::readNetworkFile(string networkName, unordered_map<string, int>& nodeMap){
	string str;
	string sourceNode, destNode, direction;
	vector<string> strs;

	ifstream in_stream;
	unordered_map<string, dynamic_bitset<>> edgeCount; //record the layers that each edge exist 
    in_stream.open(networkName);

	int layerIndex = 0;
	while(!in_stream.eof()) {
		getline(in_stream, str);
		if(str.size()>0){
			if (str[0] == '-') {
				layerIndex++;
			}
			else {
                strs.clear();
                boost::split(strs,str,boost::is_any_of("\t"));
                str = strs[0]+"\t"+strs[1];

				if (edgeCount.find(str) == edgeCount.end()) {
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	int edgeIndex = 0;
	for(auto& p: edgeCount) {
		dynamic_bitset<> tdb = p.second;
		if (tdb.count() >= k) {
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			sourceNode=strs[0];
			destNode=strs[1];
			int sourceId = nodeMap[sourceNode];
			int destId = nodeMap[destNode];

            if(sourceId != destId && adjacencyList[sourceId][destId]==0){
            	adjacencyMatrix[sourceId][destId] = edgeIndex;
                adjacencyList[sourceId][destId] = 1;
                edgeLayers.push_back(tdb);
                edgeset.push_back(make_pair(sourceNode, destNode));
                edgeIntSet.push_back(make_pair(sourceId, destId));
                edgeIndex++;
            }
		}
	}

	EdgeSize = edgeIndex;
	strs.clear();
    edgeCount.clear();
}

void Graph::findMotifs(vector<Motif*>& pattern, int caseN){
	dynamic_bitset<>::size_type it,itt,iit;
	dynamic_bitset<> ba(LayerSize);
	dynamic_bitset<> back(LayerSize);
	dynamic_bitset<> tmp(LayerSize);

	int index1,index2,index3, index4;
	index1=index2=index3=index4=-1;

	int count = 0;
	if(caseN == 1){ //Feed forward loop
		for(int i=0;i<NodeSize;i++){
			it=adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				index1 = adjacencyMatrix[i][it];
				ba = edgeLayers[index1];

				itt=adjacencyList[i].find_next(it);
				while(itt!=dynamic_bitset<>::npos){
					index2 = adjacencyMatrix[i][itt];
					tmp = edgeLayers[index2];
					tmp &= ba;
					if (tmp.count() >= k) {
						if(adjacencyList[it][itt]==1){
							index3=adjacencyMatrix[it][itt];
							back = edgeLayers[index3];
							back &= tmp;

							if (back.count() >= k) {
								dynamic_bitset<> edges(EdgeSize);
								edges[index1] = edges[index2] = edges[index3] = 1;
								Motif* t=new Motif(count, edges, back);
								pattern.push_back(t);
								count++;
							}
						}
						if(adjacencyList[itt][it]==1){
							index3=adjacencyMatrix[itt][it];
							back = edgeLayers[index3];
							back &= tmp;

							if (back.count() >= k) {
								dynamic_bitset<> edges(EdgeSize);
								edges[index1] = edges[index2] = edges[index3] = 1;
								Motif* t=new Motif(count, edges, back);
								pattern.push_back(t);
								count++;
							}
						}
					}
					itt=adjacencyList[i].find_next(itt);
				}
				it=adjacencyList[i].find_next(it);
			}
		}
	}
	else if(caseN==2){ //bifan
		for(int i=0;i<NodeSize-1;i++){
			if (adjacencyList[i].count() < 2) {
				continue;
			}
			for(int j=i+1;j<NodeSize;j++){
				if (adjacencyList[j].count() < 2) {
					continue;
				}
				dynamic_bitset<> tdb = adjacencyList[i];
				tdb &= adjacencyList[j];
				tdb[i] = 0;
				tdb[j] = 0;
				if(tdb.count()>=2){
					it = tdb.find_first();
					while(it!=dynamic_bitset<>::npos){
						index1 = adjacencyMatrix[i][it];
						index2 = adjacencyMatrix[j][it];
						tmp = edgeLayers[index1];
						tmp &= edgeLayers[index2];

						if (tmp.count() >= k) {
							itt=tdb.find_next(it);
							while(itt!=dynamic_bitset<>::npos){
								index3 = adjacencyMatrix[i][itt];
								index4 = adjacencyMatrix[j][itt];
								back = edgeLayers[index3];
								back &= edgeLayers[index4];

								if (back.count() >= k) {
									back &= tmp;
									if (back.count() >= k) {
                                        dynamic_bitset<> edges(EdgeSize);
                                        edges[index1] = edges[index2] = edges[index3] = edges[index4] = 1;
										Motif* t=new Motif(count, edges, back);
										pattern.push_back(t);
										count++;
									}
								}
								else {
									tdb[itt] = 0;
								}
								itt = tdb.find_next(itt);
							}
						}
						else {
							tdb[it] = 0;
						}
						it = tdb.find_next(it);
					}
				}
			}
		}
	}
	else if(caseN==3){ //bi-parallel
		for(int i=0;i<NodeSize;i++){
			if (adjacencyList[i].count() < 2) {
				continue;
			}

			it = adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				index1=adjacencyMatrix[i][it];
				tmp = edgeLayers[index1];
				itt = adjacencyList[i].find_next(it);
				while(itt != dynamic_bitset<>::npos){
					index2=adjacencyMatrix[i][itt];
					back = edgeLayers[index2];
					back &= tmp;

					if (back.count() >= k) {
						dynamic_bitset<> tdb=adjacencyList[it];
						tdb &= adjacencyList[itt];
						tdb[i]=0;
						tdb[it] = 0;
						tdb[itt] = 0;

						iit = tdb.find_first();
						while(iit!=dynamic_bitset<>::npos){
							index3=adjacencyMatrix[it][iit];
							index4=adjacencyMatrix[itt][iit];

							ba = edgeLayers[index3];
							ba &= edgeLayers[index4];
							ba &= back;

							if (ba.count() >= k) {
                                dynamic_bitset<> edges(EdgeSize);
                                edges[index1] = edges[index2] = edges[index3] = edges[index4] = 1;
								Motif* t=new Motif(count, edges, ba);
								pattern.push_back(t);
								count++;
							}
							iit = tdb.find_next(iit);
						}
					}

					itt = adjacencyList[i].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}
	}

	else if(caseN==4){ // cascade and delay
		for(int i=0;i < NodeSize - 2;i++){ //i is the smallest number among three nodes
			it = adjacencyList[i].find_next(i);
			while(it!=dynamic_bitset<>::npos){
				index1=adjacencyMatrix[i][it];
				tmp = edgeLayers[index1];
				itt = adjacencyList[it].find_next(i);
				while(itt != dynamic_bitset<>::npos){
					if(itt != it && adjacencyList[itt][i]==1){
						index2=adjacencyMatrix[it][itt];
						index3=adjacencyMatrix[itt][i];

						back = edgeLayers[index2];
						back &= edgeLayers[index3];
						back &= tmp;
						if (back.count() >= k) {
                            dynamic_bitset<> edges(EdgeSize);
                            edges[index1] = edges[index2] = edges[index3] = 1;
							Motif* t=new Motif(count, edges, back);
							pattern.push_back(t);
							count++;
						}
					}
					itt = adjacencyList[it].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}

	}
}
