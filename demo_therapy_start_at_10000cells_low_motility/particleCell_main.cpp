/*
    File: particleCell_main.cpp
    Model: particleCell -- random mutation model (RMM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE RMM STUDY]
*/

#include "initTumour.hpp"
#include "evolveTumour.hpp"

void writeParameterInfo(void);
int main(int argc, char** argv)
{
    // set up random number generator
    // int seed = chrono::system_clock::now().time_since_epoch().count();
    int seed = 123;
    mt19937 generator (seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);
    normal_distribution<double> normal0(0.0, sqrt(DT));
    normal_distribution<double> normalT(MITOSIS_T_DOUBLE_AVG, MITOSIS_T_DOUBLE_STD);
    normal_distribution<double> normalP(MITOSIS_CIP_THRP, MITOSIS_CIP_THRP_STD);

    // initialize voxMap
    const int dimX = DIM_VOXMAP;
    iVec dim = {dimX, dimX, dimX};
    if (SIM_DIM == "2D")
        iVec dim = {dimX, dimX, 1};
    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    VoxMap voxMap = VoxMap(dim, mapBeadIds, mapBeadIdsHalo);

    // initialize vectors for writing to file
    vector<unordered_map<string,double>> mitosisToWrite;

    // initialize vectors
    vector<Group> vecGroup; // collect all Groups (i.e., Subclones)
    vector<Node> vecNode;   // collect all Nodes (i.e., Cells)
    vector<Group> vecGroupBd;   // [NOT USED IN THE RMM STUDY]collect all Groups (i.e., different Boundaries: stroma, bronchiole, vasculature, ... Or simply container wall)
    vector<Node> vecNodeBd;     // [NOT USED IN THE RMM STUDY]collect all Nodes (i.e., components of the GroupBd)
    vector<Group> vecGroupSC;   // [NOT USED IN THE RMM STUDY]collect all explicit stroma Groups
    vector<Node> vecNodeSC;     // [NOT USED IN THE RMM STUDY]collect all explicit Stroma Cells

    vector<dVec> vecBd2; // [NOT USED IN THE RMM STUDY]collect locations (x, y, z) with boundary2 (i.e., high drag coefficient)

    // initialize groups
    initGroup(vecGroup, vecNode);

    // initialize node coordinates & velocties
    vector<double> randNum01, randNormNumLife;
    for (int i=0; i<vecNode.size()*3+1; i++)
    {
        randNum01.push_back(uniform01(generator));
        if (i<vecNode.size()+1)
        {
            if (MITOSIS_T_DOUBLE_DIST == "Gaussian")
                randNormNumLife.push_back(normalT(generator));
            if (MITOSIS_T_DOUBLE_DIST == "Erlang")
            {
                double erl, factor = 1;
                const double u_erl = MITOSIS_T_DOUBLE_STD*MITOSIS_T_DOUBLE_STD/MITOSIS_T_DOUBLE_AVG;
                const int    k_erl = (int) (MITOSIS_T_DOUBLE_AVG/u_erl);
                for (int kk = 0; kk < k_erl; kk++)
                    factor *= uniform01(generator);
                erl = -u_erl * log(factor);
                randNormNumLife.push_back(erl);
            }
        }
    }

    // default
    if (INIT_SPECIAL_PATTERN == "none")
        initDynamics(vecGroup, vecNode, randNum01, randNormNumLife);

    // [NOT USED IN THE RMM STUDY]
    if (INIT_SPECIAL_PATTERN == "three_line" || INIT_SPECIAL_PATTERN == "two_half" )
        initDynamicsSpecialPattern(vecGroup, vecNode, randNormNumLife, INIT_SPECIAL_PATTERN);
    if (INIT_SPECIAL_PATTERN == "from_saved" || INIT_SPECIAL_PATTERN == "from_saved_with_hole" || INIT_SPECIAL_PATTERN == "from_saved_keep_half"  )
        initGroupDynamicsFromSaved(vecGroup, vecNode, randNormNumLife, argc, argv);

    // (Note here: cleaned developmental code modules for initialising Boundary and Stroma)

    // SIMULATION starts
    double t_now = 0;
    long nIter = 0;
    vector<double> vecTumourSize;
    while (t_now <= T)
    {
        // one iteration of simulation
        vector<double> randNormNum0, randNormNumPthr, randNum01b;
        randNormNumLife.clear();
        randNum01.clear();
        for (int i = 0; i < vecNode.size()*7; i ++)
        {
            randNormNum0.push_back(normal0(generator));
            randNum01.push_back(uniform01(generator));
            randNum01b.push_back(uniform01(generator));
            if (i<vecNode.size()+1)
            {
                // randNum01.push_back(uniform01(generator));
                // if (MITOSIS_CIP_TYPE == "MechanicsVariable")
                //     randNormNumPthr.push_back(normalP(generator));
                if (MITOSIS_CIP_TYPE == "Mechanics")
                    randNormNumPthr.push_back(MITOSIS_CIP_THRP); // constant value used
                if (MITOSIS_T_DOUBLE_DIST == "Gaussian")
                    randNormNumLife.push_back(normalT(generator));
                if (MITOSIS_T_DOUBLE_DIST == "Erlang")
                {
                    double erl, factor = 1;
                    const double u_erl = MITOSIS_T_DOUBLE_STD*MITOSIS_T_DOUBLE_STD/MITOSIS_T_DOUBLE_AVG;
                    const int    k_erl = (int) (MITOSIS_T_DOUBLE_AVG/u_erl);
                    for (int kk = 0; kk < k_erl; kk++)
                        factor *= uniform01(generator);
                    erl = -u_erl * log(factor);
                    randNormNumLife.push_back(erl);
                }
            }
        }
        int nNode_old = vecNode.size();
        oneIterOverdamp(vecGroup, vecNode, voxMap, \
            randNormNum0, randNormNumLife, randNormNumPthr, \
            randNum01, randNum01b, mitosisToWrite, \
            vecNodeBd, vecNodeSC, vecBd2);
        int nNode_new = vecNode.size();

        // introduce chemotherapy
        if (THERAPY_ON)
        {
            if (TYPE_DRUG == "cytotoxic" || TYPE_DRUG == "cytostatic")
            {
                
                if (TYPE_DRUG_SCHEDULE == "continuous") // this schedule applies drug continuously
                {
                    bool flagTherapyConti = false;
                    // count all cells
                    if (vecNode.size() > DRUG_PER_NUM_CELL)
                        flagTherapyConti = true;

                    if (flagTherapyConti)
                    {
                        if (TYPE_DRUG_RESPONSE == "binary")
                            necrosis2(vecGroup, vecNode);

                    }
                }
                
            }
        }

        // print information
        if (nIter % FREQ_PRINT == 0)
        {
            stringstream PROC_ID_SS; PROC_ID_SS << PROC_ID;
            cout << "\n ----- pid = " << PROC_ID_SS.str() << " ----- t_now = " << t_now << " -----" << endl;

            double tumourSize;
            tumourSize = printGroupInfo(vecGroup, vecNode, t_now, seed);
            vecTumourSize.push_back(tumourSize);
        }

        // write coordinate & driver mutations to file
        if (nIter % FREQ_WRITE_DYNAMICS == 0 && FLAG_WRITE_DYNAMICS_THERAPY)
        {
            cout << "... writing cell dynamics under therapy ..." << endl;
            writeNodeDynamicsTherapy(vecNode, t_now, seed);
            cout << "done !" << endl;
        }

        t_now += DT;
        nIter ++;

    }


    return 0;
}
