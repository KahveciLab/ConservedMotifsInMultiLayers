#ifndef MOTIF_H_
#define MOTIF_H_
#include <boost/dynamic_bitset.hpp>
using namespace boost;

class Motif{
public:
     int id;
	   dynamic_bitset<> edges;
     dynamic_bitset<> layers;
     dynamic_bitset<> removes;

public:
	   Motif(int count, dynamic_bitset<> edge, dynamic_bitset<> networks): id(count){
       edges = edge;
       layers = networks;
     }

	   ~Motif(){}
};

#endif /* MOTIF_H_ */
