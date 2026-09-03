/*
    File: evolveTumour.hpp
    Model: particleCell -- random mutation model (RMM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE RMM STUDY]
*/

#ifndef EVOLVETUMOUR_HPP_INCLUDED
#define EVOLVETUMOUR_HPP_INCLUDED

#include "initTumour.hpp"

// one iteration of simulation
void oneIterOverdamp(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap, \
                     vector<double> &randNormNum0, vector<double> &randNormNumLife, vector<double> &randNormNumPthr, \
                     vector<double> &randNum01, vector<double> &randNum01b, \
                     vector<unordered_map<string,double>> &mitosisToWrite, \
                     vector<Node> &vecNodeBd, \
                     vector<Node> &vecNodeSC, \
                     vector<dVec> &vecBd2);

void updateCoord(vector<Node> & vecNode, vector<Group> & vecGroup, vector<double> &randNormNum0, \
                 vector<Node> & vecNodeSC);
void updateVeloc(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap, \
                 vector<double> &randNormNum0, \
                 vector<Node> &vecNodeBd, \
                 vector<Node> &vecNodeSC, \
             vector<dVec> &vecBd2);

void mitosis(vector<Group> &vecGroup, vector<Node> &vecNode, \
    vector<double> &randNormNumLife, vector<double> &randNormNumPthr, \
    vector<double> &randNum01, \
             vector<unordered_map<string,double>> &mitosisToWrite);
void apoptosis(vector<Group> &vecGroup, vector<Node> &vecNode);
void necrosis2(vector<Group> &vecGroup, vector<Node> &vecNode);
void recordMutations(unordered_map<int, unordered_map<string, int>> umap_mutations);
void intrinsicMutaProc(vector<Node> &vecNode, vector<double> &randNum01b);

// VoxMap algorithm
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode);

double dVecDist(dVec & p1, dVec & p2);

#endif
