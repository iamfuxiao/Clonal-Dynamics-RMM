/*
    File: initTumour.hpp
    Model: particleCell -- random mutation model (RMM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    [1] In the study Bailey_Bhargava_Fu_etal, key parameter explored in the Random Mutation Model (RMM)

    > BROWNIAN_DIFF_SCALE - which reflects a cell's intrinsic random motility term: values 0 or 1 (low) vs 6e5 (high)

    Expectation:
        higher BROWNIAN_DIFF_SCALE
            causes more fragmentation of subclones, including those harbouring therapy-resistance-confering mutations,
            and leads to more rapid dominance by resistance subclones

    [2] Other key parameters in RMM simulations explained

    > GROUP_SCALE_SEN - which with a value of 0 reflects that the cells harbouring therapy-resistance-confering mutations have zero
        sensitivity to the therapy
    > DRUG_EFFECT_RATE - which defines the death rate of cells sensitive to the therapy
    > NUM_MUT_PER_DIV or NUM_MUT_PER_DT - which define the rate of mutation accumulation in a cell, in RMM simulations
        with proliferation-dependent or proliferation-independent mutational processes implemented
    > NUM_MUT_DRIVER - which reflects the number of genes (out of the pool of 10000) that, when mutated, confer selective
        advantage in resistance to the therapy

    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE RMM STUDY]

*/

#ifndef INITTUMOUR_HPP_INCLUDED
#define INITTUMOUR_HPP_INCLUDED

#include <stdio.h>
// #include <math.h>
#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include <random>
#include <chrono>

#include <omp.h>
#include <unistd.h> // for getpid

using namespace std;

typedef struct
{
    double x, y, z;
} dVec;
typedef struct
{
    double x, y;
} dPair;
typedef struct
{
    int i, j, k;
} iVec;
typedef struct
{
    int i, j;
} iPair;

/* ~~~~~~~~~~~ Constants ~~~~~~~~~~ */
const double PI = 3.1415926535897;
const double K_BOLTZMANN = 1.38E-2; // (pN*nm/K)
const double TEMPERATURE = 310;     // (K)  body temperature
const string SIM_DIM = "2D";        // "2D" or "3D"
const pid_t PROC_ID = getpid();
const int DIM_VOXMAP = 16; // n x n x n
const int NUM_THREADS = 2;

/* ~~~~~~~~~~~ Optional Module Switchs ~~~~~~~~~~ */
const bool BROWNIAN_ON = true;    // brownian motion / stochastic force
const bool BOUNDARY_ON = false;   // [NOT USED IN THE RMM STUDY] boundary elements
const bool BOUNDARY2_ON = false;  // [NOT USED IN THE RMM STUDY] boundary modelled as a "field" that enhances local drag coefficient
const bool STROMA_ON = false;     // [NOT USED IN THE RMM STUDY] stromal cells
const bool THERAPY_ON = true;     // virtual chemotherapy
const bool PROPULSION_ON = false; // [NOT USED IN THE RMM STUDY] propulsive forces
const bool APOPTOSIS_ON = false;  // [NOT USED IN THE RMM STUDY] apoptosis

/* ~~~~~~~~~~~ Parameters ~~~~~~~~~~ */
const double DT = 0.01; // hour
const double T = 24 * 10;
const int FREQ_PRINT = (int)(24 / DT);
const int FREQ_WRITE = (int)(24 / DT);
const int FREQ_WRITE_DYNAMICS = (int)(24 / DT);
const int FREQ_WRITE_MUTATIONS = (int)(24 * 28 / DT);
const int FREQ_WRITE_MITOSIS = (int)(24 / DT);
const bool FLAG_WRITE_PDB = false;
const bool FLAG_WRITE_DYNAMICS = false;
const bool FLAG_WRITE_DYNAMICS_THERAPY = true; // this is used for writing an output in the RMM therapy simulations
const bool FLAG_WRITE_MITOSIS = false;
const bool FLAG_WRITE_MUTATIONS = false;
const double CELL_RADIUS = 15;                              // um
const double TUMOUR_RADIUS_MAX = 10000;                     // um; set to a large number, so size based stopper condition is not used.
const vector<int> VEC_TUMOUR_RADIUS_WRITE_MUTATIONS{12000}; // set to a large number, so detailed mutations are not saved in the RMM therapy simulations.

// Initial configuration
const string INIT_SPECIAL_PATTERN = "none"; // "none" - starting from single cell
const int INIT_NUM_GROUP = 1;
const int INIT_NUM_NODE_PER_GROUP = 1; // 1 if special pattern is none; 60 for three_line; 2500 for two_half
const double INIT_TUMOUR_RADIUS = cbrt(INIT_NUM_GROUP * INIT_NUM_NODE_PER_GROUP) * CELL_RADIUS;
const double INIT_TUMOUR_RADIUS_2D = sqrt(INIT_NUM_GROUP * INIT_NUM_NODE_PER_GROUP) * CELL_RADIUS;
const double INIT_FROM_SAVED_HOLE_RADIUS = 1000; // [NOT USED IN THE RMM STUDY] um   only if INIT_SPECIAL_PATTERN = "from_saved_with_hole"

// Brownian motion
const double GAMMA_VISCOUS = 1E2; // drag coefficient
// -- [ Key parameter studied in Bailey_Bhargava_Fu_etal RMM simulations]
const double BROWNIAN_DIFF_SCALE = 0;
const double BROWNIAN_DIFF_COEF_SQRT = sqrt(BROWNIAN_DIFF_SCALE * 2 * K_BOLTZMANN * TEMPERATURE * 3600 * 3600 / GAMMA_VISCOUS * 1E-9); // sqrt(um*um/h)

// Adhesion & Repulsion
const double ADHESION_LENGTH_ACT = 2 * CELL_RADIUS;   // um
const double ADHESION_LENGTH_EQU = 1.5 * CELL_RADIUS; // um
const double ADHESION_AREA_MIN = 1;                   // um^2

// ... Harmonic potential
const double ADHESION_STRENGTH = 1000;  // kg/h^2
const double REPULSION_STRENGTH = 1000; // kg/h^2
const double REPULSION_LENGTH_ACT = 2 * CELL_RADIUS;

// Proliferation & Mitosis
const double MITOSIS_T_DOUBLE_AVG = 10;        // h; if it's erlang, mean=ku
const double MITOSIS_T_DOUBLE_STD = 2;         // h; if it's erlang, var=ku^2
const string MITOSIS_T_DOUBLE_DIST = "Erlang"; // "Erlang", "Gaussian"
const string MITOSIS_CIP_TYPE = "Mechanics";   // "Mechanics"
const string MITOSIS_ENP_TYPE = "none";        // [NOT USED IN THE RMM STUDY] "Inherit", "Mechanics", "Geometry"
const int MITOSIS_CIP_THR_3D = 12;
const int MITOSIS_CIP_THR_2D = 6;                        // default: 6 for 2D, 12 for 3D; number of surrounding cells needs to be smaller than or equal to this number to permit mitosis
const double MITOSIS_CIP_THRP = -17 / CELL_RADIUS * 15.; //
const double MITOSIS_CIP_THRP_STD = 2;
const int PERIPHERY_MAX_NCB = 5; // [NOT USED IN THE RMM STUDY]

// Subclone evolution
const string TYPE_SUBCL = "emerge"; // "induce" - lineage tracing model (LTM), "emerge" - random mutation model (RMM), "none"
// .. Subclone induction (chemical exposure)
const int NUM_SUBCL = 10;                      // [NOT USED IN THE RMM STUDY]
const string TYPE_INDUCE_SUBCL = "stochastic"; // [NOT USED IN THE RMM STUDY]
const string TYPE_INDUCE_SUBCL_WHEN = "nhour"; // [NOT USED IN THE RMM STUDY]
const int NUM_CELL_INDUCE_SUBCL = 1000;        // [NOT USED IN THE RMM STUDY]
const int NUM_CELL_INDUCE_SUBCL_DET = 1000;    // [NOT USED IN THE RMM STUDY]
const int NUM_HOUR_INDUCE_SUBCL = 120;         // [NOT USED IN THE RMM STUDY]
const double PERC_CELL_INDUCE_SUBCL = 0.1;     // [NOT USED IN THE RMM STUDY]
// .. Subclone emergence (mutation)
const string TYPE_EMERGE_SUBCL = "random";                        // "random" (prolif. dependent); "random2" (prolif. independent)
const double RATE_MUT_PER_DIV = 5E-3;                             // [NOT USED IN THE RMM STUDY]
const double PROB_MUT_DRIVER = 5E-4;                              // probability that a mutation is a driver mutation (ONLY valid if it's TYPE_SUBCL is "emerge")
const int NUM_MUT_POOL = 10000;                                   // number of possible mutations in a candidate pool (SINCE 2019.03.18)
const int NUM_MUT_PER_DIV = 5;                                    // number of mutations per cell division (SINCE 2019.03.18)
const int NUM_MUT_DRIVER = (int)(PROB_MUT_DRIVER * NUM_MUT_POOL); // number of mutations being driver
const int NUM_MUT_DRIVER_DEL = (int)(0.5 * NUM_MUT_POOL);         // [NOT USED IN THE RMM STUDY] number of mutations being deleterious

const double RATE_MUT_PER_DT = NUM_MUT_PER_DIV / (MITOSIS_T_DOUBLE_AVG / DT); // used only if TYPE_EMERGE_SUBCL == "random2"

// [NOT USED IN THE RMM STUDY] Brownian rotation + active propulsion
const double BROWNIAN_DIFF_COEF_SQRT_PROP_ANGLE = 0.01; // sqrt(deg*deg/h)
const double PROPULSION_SIZE = 100;                     // pN

// Chemotherapy
const string TYPE_DRUG = "cytotoxic";             // "cytotoxic" (attack all sensitive cells)
const string TYPE_DRUG_SCHEDULE = "continuous";   // "continuous" (alway on)
const string TYPE_DRUG_RESPONSE = "binary";       // "binary" (any driver mutation confers complete drug resistance)
const double DRUG_DOSAGE = 1;                     // 1 is max dosage
const double DRUG_PERIOD = 24 * 7;                // hours
const int DRUG_PER_NUM_CELL = 10000;              // threshold size if used for "per_num_cell" ; minimum size if used for TYPE_DRUG_SCHEDULE == "continuous"
const int DRUG_PER_NUM_CELL_ADAPTIVE_STOP = 5000; // size to end therapy if TYPE_DRUG_SCHEDULE == "adaptive"

const int NUM_HOUR_INTRO_DRUG = DRUG_PERIOD; // hours post clonal induction
const double DRUG_EFFECT_RATE = 0.2;         // 1/h; // ONLY in use when TYPE_DRUG_SCHEDULE == "continuous"
// [NOT USED IN THE RMM STUDY] Mixing-modulating agent
const string TYPE_DRUG_2 = "promotility"; // "promotility" (enhances the BROWNIAN_DIFF_SCALE by a factor)
const double DRUG_2_IMPACT_DIFF = 6e5;    // 4e5 used for mix2a, 6e5 used for mix2b
const double DRUG_2_IMPACT_PTHR = -5;     // -5 used for mix2

// [NOT USED IN THE RMM STUDY] Apoptosis
const double APOPTOSIS_RATE = 1e-3; // 1/h;

// [OPTIONAL] Group-specific Scaling coefficients
// "MITOSIS", "MITOSIS_PTHR" (NOT WORKING), "MITOSIS2" (both beneficial & deleterious mutations; only implemented for RMM)
// "ADHESION", "FRICTION", "SENSITIVITY_DRUG", "SENSITIVITY_DRUG_WITH_COST", "APOPTOSIS", "RADIUS", "DIFFUSION"
const string TYPE_DRIVER_MUTATION = "SENSITIVITY_DRUG";
const double GROUP_SCALE_VIS = 1;  // [NOT USED IN THE RMM STUDY]scale GAMMA_VISCOUS
const double GROUP_SCALE_ADH = 1;  // [NOT USED IN THE RMM STUDY]scale ADHESION_STRENGTH
const double GROUP_SCALE_MIT = 1;  // scale MITOSIS_T_DOUBLE [or MITOSIS_CIP_THRP (not working!)]
const double GROUP_SCALE_MIT2 = 1; // scale MITOSIS_T_DOUBLE when there is a second type of mutation in RMM
const double GROUP_SCALE_SEN = 0;  // scale drug sensitivity 1 or 0
const double GROUP_SCALE_APO = 1;  // [NOT USED IN THE RMM STUDY]scale APOPTOSIS_RATE
const double GROUP_SCALE_RAD = 1;  // [NOT USED IN THE RMM STUDY]scale CELL_RADIUS
const double GROUP_SCALE_DIF = 1;  // [NOT USED IN THE RMM STUDY]scale BROWNIAN_DIFF_SCALE

// [NOT USED IN THE RMM STUDY]  Boundary effects
const string TYPE_BOUNDARY = "random";                       // "random", "vertical", "radial", "circum", "radcir", "ring"
const int INIT_NUM_GROUPBD = 200;                            // number of groupGd
const int INIT_NUM_NODE_PER_GROUPBD = 25;                    // number of nodes per groupBd
const double POS_BOUNDARY_LINE_X = 200;                      // location of a line-shaped boundary; only used if is single "vertical" boundary
const double INIT_RING_DIST2CEN = 2000;                      // distance of the ring boundary to center
const int INIT_RING_GAP_SIZE = 20;                           // gap size on the ring boundary
const double INIT_RING_BAR_GAP_RATIO = 1;                    // barrier to gap ratio
const double REPULSION_STRENGTH_BD = 2 * REPULSION_STRENGTH; // assumpting twice the cell-cell repulsion

// [NOT USED IN THE RMM STUDY]  Boundary effects 2
const string TYPE_BOUNDARY2 = "random";
const int INIT_NUM_GROUPBD2 = 100; // number of groupGd2
const int INIT_BD2_THICKNESS = 40;
const int INIT_BD2_LENGTH = 500;
const double BD2_SCALE_VIS_BY = 10; // scale the drag coefficient by this factor

// [NOT USED IN THE RMM STUDY] Stromal cells -- only implemented in lineage tracing model
const int INIT_SC_NCELL_SIDE = 100;  // number of cells for a square patch
const float INIT_SC_FRAC_KEEP = 0.5; // fraction of cells (in full hexagonal lattice) to be kept

/* ~~~~~~~~~~~ Classes ~~~~~~~~~~ */
/* -------------- Node (representing Cell) --------------- */
class Node
{
    bool isAlive;         // isAlive is a flag to indicate whether the cell is dead (apoptotic or necrotic)
    int id, groupId, gen; // gen refers to generation; founder cell has generation 0
    int lineageId;        // every clonal founder cell has a unique lineageId
    dVec coord, coord_prev, veloc;
    double propu_angle;
    double time_double, time_double_intrinsic;
    vector<int> vecNbNodeIds;
    vector<double> vecNbNodeDists;

    vector<int> vecMutIds; // record a list of mutations accumulated in the cell

    double pressure, radius;
    double pthr;          // assume pthr is cell-specific
    double strain_energy; // e.g., half the sum of potential energy in the volume exclusion spring
    vector<dPair> stress2d;
    vector<double> vecNbCommAreas;
    vector<double> vecNbPressures;
    unordered_map<string, int> umapCellDensity; // string indicates range; int indicates count

public:
    Node();
    Node(bool isAlive, int id, int lineageId, int groupId);
    Node(bool isAlive, int id, int lineageId, int groupId, int gen,
         dVec coord, dVec coord_prev, dVec veloc, double propu_angle,
         double time_double, double time_double_intrinsic,
         vector<int> vecNbNodeIds, vector<double> vecNbNodeDists,
         vector<int> vecMutIds,
         double pressure, double radius, double pthr,
         vector<dPair> stress2d,
         double strain_energy,
         vector<double> vecNbCommAreas, vector<double> vecNbPressures,
         unordered_map<string, int> umapCellDensity);
    ~Node();

    // setters
    void set_isAlive(bool isAlive) { this->isAlive = isAlive; }
    void set_id(int id) { this->id = id; }
    void set_lineageId(int lineageId) { this->lineageId = lineageId; }
    void set_groupId(int groupId) { this->groupId = groupId; }
    void set_gen(int gen) { this->gen = gen; }
    void set_coord(dVec coord) { this->coord = coord; }
    void set_coord_prev(dVec coord_prev) { this->coord_prev = coord_prev; }
    void set_veloc(dVec veloc) { this->veloc = veloc; }
    void set_propu_angle(double propu_angle) { this->propu_angle = propu_angle; }
    void set_time_double(double time_double) { this->time_double = time_double; }
    void set_time_double_intrinsic(double time_double_intrinsic) { this->time_double_intrinsic = time_double_intrinsic; }
    void set_vecNbNodeIds(vector<int> vecNbNodeIds) { this->vecNbNodeIds = vecNbNodeIds; }
    void set_vecNbNodeDists(vector<double> vecNbNodeDists) { this->vecNbNodeDists = vecNbNodeDists; }
    void set_vecMutIds(vector<int> vecMutIds) { this->vecMutIds = vecMutIds; }
    void set_pressure(double pressure) { this->pressure = pressure; }
    void set_pthr(double pthr) { this->pthr = pthr; }
    void set_strain_energy(double strain_energy) { this->strain_energy = strain_energy; }
    void set_radius(double radius) { this->radius = radius; }
    void set_stress2d(vector<dPair> stress2d) { this->stress2d = stress2d; }
    void set_vecNbCommAreas(vector<double> vecNbCommAreas) { this->vecNbCommAreas = vecNbCommAreas; }
    void set_vecNbPressures(vector<double> vecNbPressures) { this->vecNbPressures = vecNbPressures; }
    void set_umapCellDensity(unordered_map<string, int> umapCellDensity) { this->umapCellDensity = umapCellDensity; }

    // getters
    bool get_isAlive() const { return this->isAlive; }
    int get_id() const { return this->id; }
    int get_lineageId() const { return this->lineageId; }
    int get_groupId() const { return this->groupId; }
    int get_gen() const { return this->gen; }
    dVec get_coord() const { return this->coord; }
    dVec get_coord_prev() const { return this->coord_prev; }
    dVec get_veloc() const { return this->veloc; }
    double get_propu_angle() const { return this->propu_angle; }
    double get_time_double() const { return this->time_double; }
    double get_time_double_intrinsic() const { return this->time_double_intrinsic; }
    vector<int> get_vecNbNodeIds() const { return this->vecNbNodeIds; }
    vector<double> get_vecNbNodeDists() const { return this->vecNbNodeDists; }
    vector<int> get_vecMutIds() const { return this->vecMutIds; }
    double get_pressure() const { return this->pressure; }
    double get_pthr() const { return this->pthr; }
    double get_strain_energy() const { return this->strain_energy; }
    double get_radius() const { return this->radius; }
    vector<dPair> get_stress2d() const { return this->stress2d; }
    vector<double> get_vecNbCommAreas() const { return this->vecNbCommAreas; }
    vector<double> get_vecNbPressures() const { return this->vecNbPressures; }
    unordered_map<string, int> get_umapCellDensity() const { return this->umapCellDensity; }
};

/* -------------- Group (representing Subclone) ------------- */
class Group
{
    int id;
    vector<int> vecNodeIds;
    double scale_vis, scale_adh, scale_mit, scale_sen, scale_dif; // this is scaling coefficient on the global parameter value

public:
    Group();
    Group(int id, vector<int> vecNodeIds, double scale_vis, double scale_adh, double scale_mit, double scale_sen, double scale_dif);
    ~Group();

    // setters
    void set_id(int id) { this->id = id; }
    void set_vecNodeIds(vector<int> vecNodeIds) { this->vecNodeIds = vecNodeIds; }
    void set_scale_vis(double scale_vis) { this->scale_vis = scale_vis; }
    void set_scale_adh(double scale_adh) { this->scale_adh = scale_adh; }
    void set_scale_mit(double scale_mit) { this->scale_mit = scale_mit; }
    void set_scale_sen(double scale_sen) { this->scale_sen = scale_sen; }
    void set_scale_dif(double scale_dif) { this->scale_dif = scale_dif; }

    // getters
    int get_id() const { return this->id; }
    vector<int> get_vecNodeIds() const { return this->vecNodeIds; }
    double get_scale_vis() const { return this->scale_vis; }
    double get_scale_adh() const { return this->scale_adh; }
    double get_scale_mit() const { return this->scale_mit; }
    double get_scale_sen() const { return this->scale_sen; }
    double get_scale_dif() const { return this->scale_dif; }
};

/* -------------- VoxMap ----------------- */
class VoxMap
{
    iVec dim;
    unordered_map<int, vector<int>> mapBeadIds;     // contains Bead ids in VoxBox; key is collapsed from 3d indexing
    unordered_map<int, vector<int>> mapBeadIdsHalo; // contains Bead ids in halo of VoxBox;

public:
    VoxMap();
    VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo);
    ~VoxMap();

    // setters
    void set_dim(iVec dim) { this->dim = dim; }
    void set_mapBeadIds(unordered_map<int, vector<int>> mapBeadIds) { this->mapBeadIds = mapBeadIds; }
    void set_mapBeadIdsHalo(unordered_map<int, vector<int>> mapBeadIdsHalo) { this->mapBeadIdsHalo = mapBeadIdsHalo; }

    // getters
    iVec get_dim() const { return this->dim; }
    unordered_map<int, vector<int>> get_mapBeadIds() const { return this->mapBeadIds; }
    unordered_map<int, vector<int>> get_mapBeadIdsHalo() const { return this->mapBeadIdsHalo; }
};

/* ~~~~~~~~~~~ Functions ~~~~~~~~~~ */
void initGroup(vector<Group> &vecGroup, vector<Node> &vecNode);
void initGroupBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd);
void initDynamics(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNum01, vector<double> &randNormNumLife);
void initDynamicsSpecialPattern(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, string typePattern);
void createLineOfCoords(vector<dVec> &vecLineCoord, int ncell, dVec orientation, double stepsize);
void initGroupDynamicsFromSaved(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, int argc, char **argv);
void initDynamicsBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, vector<double> &randNum01);
void initDynamicsBdSpecialPattern(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, string typePattern);
void initStroma(vector<Group> &vecGroupSC, vector<Node> &vecNodeSC);
void createLineOfCoords(vector<dVec> &vecBd2, dVec orientation);
void initBd2(vector<dVec> &vecBd2, vector<double> &randNum01);

void currCloneCenterPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery, dVec &cloneCOM, bool flagCheckPeriphery);
void sortClonePeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery);
void fillGapPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery);
void calcCurvature(vector<Node> &vecNode, vector<int> &vecNodePeriphery, unordered_map<int, double> &mapNodeCurvature);

double printGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, int seed);
void writeNodeDynamicsTherapy(vector<Node> &vecNode, double t_now, int seed);

double dVecDist(dVec &p1, dVec &p2);

#endif
