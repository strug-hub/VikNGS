#pragma once
#include "../Math/Math.h"
#include "../Enum/Depth.h"
#include <set>


class Group {
    VectorXi G;
    std::map<int, Depth> d;
    int ngroup;

    std::set<int> groupIDs;

    void updateCount(){
        groupIDs.clear();
        for(int i = 0; i< G.rows(); i++)
            groupIDs.insert(G[i]);
    }

public:

    Group(VectorXi groupID, std::map<int, Depth>& depthMap) : G(groupID), d(depthMap) {
        updateCount();
    }

    //toRemove: 0 = keep, 1 = remove
    inline void filterG(VectorXi& toRemove){
        this->G = extractRows(this->G, toRemove, 0);
        updateCount();
    }

    inline int operator[](int i) { return G[i]; }
    inline Depth depth(int groupID) { return d.at(groupID); }
    inline Depth depthAt(int i) { return d.at(G[i]); }
    inline VectorXi* getG() { return &G; }
    inline int countGroups() { return this->groupIDs.size(); }
    inline int size() { return this->G.rows(); }

};
