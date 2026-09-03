/*
    File: initTumour.cpp
    Model: particleCell -- random mutation model (RMM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE RMM STUDY]
*/

#include "initTumour.hpp"

/* ~~~~~~~~~~~ Classes ~~~~~~~~~~ */
/* ------------ Node ------------- */
Node::Node() {}
Node::Node(bool isAlive, int id, int lineageId, int groupId)
{
    this->isAlive = isAlive;
    this->id = id;
    this->lineageId = lineageId;
    this->groupId = groupId;
}
Node::Node(bool isAlive, int id, int lineageId, int groupId, int gen,
           dVec coord, dVec coord_prev, dVec veloc, double propu_angle,
           double time_double, double time_double_intrinsic,
           vector<int> vecNbNodeIds, vector<double> vecNbNodeDists,
           vector<int> vecMutIds,
           double pressure, double radius, double pthr,
           vector<dPair> stress2d,
           double strain_energy,
           vector<double> vecNbCommAreas, vector<double> vecNbPressures,
           unordered_map<string, int> umapCellDensity)
{
    this->isAlive = isAlive;
    this->id = id;
    this->lineageId = lineageId;
    this->groupId = groupId;

    this->gen = gen;
    this->coord = coord;
    this->coord_prev = coord_prev;
    this->veloc = veloc;
    this->propu_angle = propu_angle;

    this->time_double = time_double;
    this->time_double_intrinsic = time_double_intrinsic;
    this->vecNbNodeIds = vecNbNodeIds;
    this->vecNbNodeDists = vecNbNodeDists;

    this->vecMutIds = vecMutIds;

    this->pressure = pressure;
    this->pthr = pthr;
    this->radius = radius;
    this->stress2d = stress2d;

    this->strain_energy = strain_energy;

    this->vecNbCommAreas = vecNbCommAreas;
    this->vecNbPressures = vecNbPressures;
    this->umapCellDensity = umapCellDensity;
}
Node::~Node() {}

/* ------------- Group -------------- */
Group::Group() {}
Group::Group(int id, vector<int> vecNodeIds, double scale_vis, double scale_adh, double scale_mit, double scale_sen, double scale_dif)
{
    this->id = id;
    this->vecNodeIds = vecNodeIds;
    this->scale_vis = scale_vis;
    this->scale_adh = scale_adh;
    this->scale_mit = scale_mit;
    this->scale_sen = scale_sen;
    this->scale_dif = scale_dif;
}
Group::~Group() {}

/* ------------- VoxMap ---------------- */
VoxMap::VoxMap() {}
VoxMap::VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo)
{
    this->dim = dim;
    this->mapBeadIds = mapBeadIds;
    this->mapBeadIdsHalo = mapBeadIdsHalo;
}
VoxMap::~VoxMap() {}

/* ~~~~~~~~~~~ Functions ~~~~~~~~~~ */
// function initGroup() is a one-time function
void initGroup(vector<Group> &vecGroup, vector<Node> &vecNode)
{
    bool isAlive = true;
    int currGroupId = vecGroup.size() + 1;
    int currNodeId = vecNode.size();
    int currLineageId = 0;
    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;

    for (int iG = 0; iG < INIT_NUM_GROUP; iG++)
    {
        vector<int> vecNodeIds;
        for (int iN = 0; iN < INIT_NUM_NODE_PER_GROUP; iN++)
        {
            Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
            vecNode.push_back(node0);
            vecNodeIds.push_back(currNodeId);
            currNodeId++;
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
        vecGroup.push_back(group0);

        cout << "... init Group ID = " << currGroupId << " ..." << endl;
        cout << "with scale_vis = " << vecGroup[currGroupId - 1].get_scale_vis() << endl;
        cout << "     scale_adh = " << vecGroup[currGroupId - 1].get_scale_adh() << endl;
        cout << "     scale_mit = " << vecGroup[currGroupId - 1].get_scale_mit() << endl;
        cout << "     scale_sen = " << vecGroup[currGroupId - 1].get_scale_sen() << endl;
        cout << "     scale_dif = " << vecGroup[currGroupId - 1].get_scale_dif() << endl;

        currGroupId++;
    }
}
// [NOT USED IN THE RMM STUDY] function initGroupBd() is a one-time function that creates initial configuration of boundary
void initGroupBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd)
{
    bool isAlive = true;
    int currGroupId = vecGroupBd.size() + 1;
    int currNodeId = vecNodeBd.size();
    int currLineageId = -1;
    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;

    // const int INIT_NUM_NODE_PER_GROUPBD = 100;

    // for (int iG = 0; iG < INIT_NUM_GROUP; iG ++)
    for (int iG = 0; iG < INIT_NUM_GROUPBD; iG++)
    {
        vector<int> vecNodeIds;
        for (int iN = 0; iN < INIT_NUM_NODE_PER_GROUPBD; iN++)
        {
            Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
            vecNodeBd.push_back(node0);
            vecNodeIds.push_back(currNodeId);
            currNodeId++;
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
        vecGroupBd.push_back(group0);

        cout << "... init GroupBd ID = " << currGroupId << " ..." << endl;
        currGroupId++;
    }
}

// function initDynamics() is a one-time function
void initDynamics(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNum01, vector<double> &randNormNumLife)
{
    double rInit = INIT_TUMOUR_RADIUS;
    double rInit2D = INIT_TUMOUR_RADIUS_2D;
    double r, theta, phi;
    int k = 0, kk = 0;

    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        // ref1: https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        // ref2: http://corysimon.github.io/articles/uniformdistn-on-sphere/
        // coordinate
        dVec coord;
        if (SIM_DIM == "3D")
        {
            r = rInit * cbrt(randNum01[k++]);
            theta = acos(2. * randNum01[k++] - 1);
            phi = 2. * PI * randNum01[k++];
            coord = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
        }
        if (SIM_DIM == "2D")
        {
            r = rInit2D * sqrt(randNum01[k++]);
            phi = 2. * PI * randNum01[k++];
            coord = {r * cos(phi), r * sin(phi), 0.};
        }
        it->set_coord(coord);
        it->set_coord_prev(coord);

        // velocity
        it->set_veloc({0., 0., 0.});

        // propulsion angle
        double propu_angle = rand() / (float)RAND_MAX * 2 * PI; // need to use random generator !
        it->set_propu_angle(propu_angle);
        // cout << propu_angle << endl;

        // time_double
        double time_double = randNormNumLife[kk++];
        it->set_time_double(time_double);
        it->set_time_double_intrinsic(time_double);

        // generation
        it->set_gen(0);

        // accumulate mutations
        vector<int> vecMutIds;
        it->set_vecMutIds(vecMutIds);

        // mechanics
        it->set_pressure(0.);
        it->set_pressure(MITOSIS_CIP_THRP);
        it->set_radius(CELL_RADIUS);
        it->set_strain_energy(0.);
    }
}
// [NOT USED IN THE RMM STUDY] function initDynamicsSpecialPattern() is a one-time function
void initDynamicsSpecialPattern(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, string typePattern)
{
    if (typePattern == "three_line")
    {
        // (1) create three cells to form a triangle
        // (2) elongate these cells along a line
        int nCellPerLine = vecNode.size() / 3;
        vector<dVec> vecLineCoord1, vecLineCoord2, vecLineCoord3;
        // ... verticle line
        vecLineCoord1.push_back({0, sqrt(3) / 3. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord1, nCellPerLine, {0, 1, 0}, 2 * CELL_RADIUS);
        // ... lower left line
        vecLineCoord2.push_back({-CELL_RADIUS, -sqrt(3) / 6. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord2, nCellPerLine, {-sqrt(3) / 2, -0.5, 0}, 2 * CELL_RADIUS);
        // ... lower right line
        vecLineCoord3.push_back({CELL_RADIUS, -sqrt(3) / 6. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord3, vecNode.size() - 2 * nCellPerLine, {sqrt(3) / 2, -0.5, 0}, 2 * CELL_RADIUS);

        int kk = 0;
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
        {
            dVec coord;
            if (it - vecNode.begin() < nCellPerLine)
                coord = vecLineCoord1[it - vecNode.begin()];
            else if (it - vecNode.begin() < 2 * nCellPerLine)
                coord = vecLineCoord2[(it - vecNode.begin()) % nCellPerLine];
            else
                coord = vecLineCoord3[(it - vecNode.begin()) % (2 * nCellPerLine)];

            it->set_coord(coord);
            it->set_coord_prev(coord);

            // velocity
            it->set_veloc({0., 0., 0.});

            // time_double
            double time_double = randNormNumLife[kk++];
            it->set_time_double(time_double);
            it->set_time_double_intrinsic(time_double);

            // generation
            it->set_gen(0);

            // mechanics
            it->set_pressure(0.);
            it->set_radius(CELL_RADIUS);
            it->set_strain_energy(0.);
        }
    }
    if (typePattern == "two_half")
    {
        // (1) define a square region of cells that fit into the square region
        // (2) assign the right half as Group 2
        const int ncell_side = (int)sqrt(INIT_NUM_NODE_PER_GROUP);
        const double sep = 2. * CELL_RADIUS;
        dVec coord_lower_left;
        coord_lower_left = {-ncell_side * 0.5 * sep, -ncell_side * 0.5 * sep * sqrt(3) * 0.5, 0};

        vector<int> vecNodeIdsChangeGroup;
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
        {
            dVec coord;
            int row = (it - vecNode.begin()) / ncell_side;
            int col = (it - vecNode.begin()) % ncell_side;

            if (row % 2 == 0)
                coord = {coord_lower_left.x + col * sep, coord_lower_left.y + row * sep * sqrt(3) * 0.5, 0};
            else if (row % 2 == 1)
                coord = {coord_lower_left.x + (col + 0.5) * sep, coord_lower_left.y + row * sep * sqrt(3) * 0.5, 0};

            it->set_coord(coord);
            it->set_coord_prev(coord);

            // velocity
            it->set_veloc({0., 0., 0.});

            // time_double
            double time_double = 10 * T; // make sure that growth is absent!
            it->set_time_double(time_double);
            it->set_time_double_intrinsic(time_double);

            // generation
            it->set_gen(0);

            // mechanics
            it->set_pressure(0.);
            it->set_radius(CELL_RADIUS);
            it->set_strain_energy(0.);

            // Move cells in the right half to Group 2
            if (coord.x > 0)
            {
                vecNodeIdsChangeGroup.push_back(it->get_id());
                // cout << it->get_id() << endl;
            }
            // else
            // cout << it->get_id() << ", " << coord.x << endl;
        }

        // ... step 1: remove nId from original Group
        // ... step 2: create a new Group
        // ... step 3: associate nId with new Group
        for (vector<int>::iterator it = vecNodeIdsChangeGroup.begin(); it != vecNodeIdsChangeGroup.end(); it++)
        {
            int nId = *it;
            int lId = (it - vecNodeIdsChangeGroup.begin()) + 1;
            int gId = vecNode[nId].get_groupId();
            // cout << nId << endl;

            // ... step 1 ...
            vector<int> vecNodeIds = vecGroup[gId - 1].get_vecNodeIds();
            vecNodeIds.erase(find(vecNodeIds.begin(), vecNodeIds.end(), nId));
            vecGroup[gId - 1].set_vecNodeIds(vecNodeIds);

            // ... step 2 ...
            double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;
            int currGroupId = vecGroup.size() + 1;
            vector<int> vecNodeIds_new;
            vecNodeIds_new.push_back(nId);

            if (currGroupId == 2)
            {
                Group group_new = Group(currGroupId, vecNodeIds_new, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
                vecGroup.push_back(group_new);
            }
            else
            {
                int joinGroupId = 2;
                vector<int> vecNodeIds2 = vecGroup[joinGroupId - 1].get_vecNodeIds();
                vecNodeIds2.push_back(nId);
                vecGroup[joinGroupId - 1].set_vecNodeIds(vecNodeIds2);

                currGroupId = joinGroupId;
            }

            // ... step 3 ...
            vecNode[nId].set_groupId(currGroupId);
            vecNode[nId].set_lineageId(lId);
        }
    }
}
// [NOT USED IN THE RMM STUDY] function createLineOfCoords() is a one-time function
// ... create a vector of (x,y,z) for BOUNDARY
void createLineOfCoords(vector<dVec> &vecLineCoord, int ncell, dVec orientation, double stepsize)
{
    // Note: ensure that vecLineCoord has 1 coord already!
    int cnt = 1;
    while (cnt < ncell)
    {
        dVec coord0 = vecLineCoord[cnt - 1];
        dVec coord1 = {coord0.x + orientation.x * stepsize,
                       coord0.y + orientation.y * stepsize,
                       coord0.z + orientation.z * stepsize};
        vecLineCoord.push_back(coord1);
        cnt++;
    }
}

// [NOT USED IN THE RMM STUDY] function createLineOfCoords2() is a one-time function
// ... create a vector of (x,y,z) for BOUNDARY2
void createLineOfCoords2(vector<dVec> &vecBd2, dVec orientation)
{
    // Note: ensure that vecLineCoord has 1 coord already!
    int cnt = 1, ncell = (int)INIT_BD2_LENGTH;
    int stepsize = 1;
    while (cnt < ncell)
    {
        dVec coord0 = vecBd2[vecBd2.size() - 1];
        dVec coord1 = {coord0.x + orientation.x * stepsize,
                       coord0.y + orientation.y * stepsize,
                       coord0.z + orientation.z * stepsize};
        vecBd2.push_back(coord1);
        cnt++;
    }
}

// [NOT USED IN THE RMM STUDY] function initGroupDynamicsFromSaved is a one-time function
void initGroupDynamicsFromSaved(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, int argc, char **argv)
{
    if (argc == 1 || argc == 2)
    {
        cout << "\nPlease input BOTH filename & time as arguments for reconstructing saved tumour! \n";
        exit(0);
    }
    if (argc == 3)
    {
        vecGroup.clear();
        vecNode.clear();

        // read from file and reconstruct tumour
        string filename = argv[1], tstr = argv[2];
        int tsaved = stoi(tstr);
        cout << ">> reading file " << filename << " ... " << endl;
        ifstream data(filename);
        string line;

        // create Groups
        int currGroupId = vecGroup.size() + 1;
        for (int iG = 0; iG < 20; iG++)
        {
            vector<int> vecNodeIds;
            double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;

            Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
            vecGroup.push_back(group0);

            currGroupId++;
        }

        // create Nodes
        if (data.is_open())
        {
            getline(data, line);
            while (getline(data, line))
            {
                stringstream ss(line);
                vector<string> items;
                string buf;
                while (ss >> buf)
                    items.push_back(buf);
                cout << line << endl;

                int t = stoi(items[0]);

                if (t == tsaved)
                {
                    int nId = stoi(items[1]);
                    int gId = stoi(items[3]);
                    int lId = stoi(items[2]);
                    double x = stod(items[5]);
                    double y = stod(items[6]);
                    double z = stod(items[7]);
                    cout << ">> creating Node ID = " << nId << endl;
                    cout << "... Group ID = " << gId << endl;
                    cout << "... Lineage ID = " << lId << endl;
                    cout << "... (x,y,z) = (" << x << ", " << y << ", " << z << ")" << endl;

                    // create Node & Group
                    bool isAlive = true;
                    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1.;
                    int currGroupId = gId;

                    Node node0 = Node(isAlive, nId, lId, currGroupId);
                    // coordinate
                    node0.set_coord({x, y, z});
                    node0.set_coord_prev({x, y, z});
                    // velocity
                    node0.set_veloc({0., 0., 0.});
                    // propulsion angle
                    double propu_angle = rand() / (float)RAND_MAX * 2 * PI; // need to use random generator !
                    node0.set_propu_angle(propu_angle);
                    // time_double
                    // double time_double = randNormNumLife[kk++];
                    double time_double = MITOSIS_T_DOUBLE_AVG;
                    node0.set_time_double(time_double);
                    node0.set_time_double_intrinsic(time_double);
                    // generation
                    node0.set_gen(0);
                    // accumulate mutations
                    vector<int> vecMutIds;
                    node0.set_vecMutIds(vecMutIds);
                    // mechanics
                    node0.set_pressure(0.);
                    node0.set_radius(CELL_RADIUS);
                    node0.set_strain_energy(0.);

                    vecNode.push_back(node0);
                    vector<int> vecNodeIds = vecGroup[currGroupId - 1].get_vecNodeIds();
                    vecNodeIds.push_back(nId);
                    vecGroup[currGroupId - 1].set_vecNodeIds(vecNodeIds);
                }

                if (t > tsaved)
                    break;
            }
        }

        // create a "hole" if INIT_SPECIAL_PATTERN = "from_saved_with_hole"
        // ... basically, set these cells in the hole as "dead" cells
        if (INIT_SPECIAL_PATTERN == "from_saved_with_hole")
        {
            cout << ">> creating hole in the center " << endl;
            for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
            {
                bool isAlive = it->get_isAlive();
                if (isAlive)
                {
                    // calculate distance to center
                    dVec coord = it->get_coord();
                    double dsq = coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;

                    if (dsq <= INIT_FROM_SAVED_HOLE_RADIUS * INIT_FROM_SAVED_HOLE_RADIUS) // condition for cell removal to occur
                    {
                        // ... setting isAlive to false
                        it->set_isAlive(false);

                        // ... remove from vecNodeIds of vecGroup
                        int nId = it->get_id();
                        int gId = it->get_groupId();

                        cout << "... removing cell ID = " << nId << " in group ID = " << gId << endl;

                        // cout << "checkpoint 0 .. " << endl;
                        vector<int> vecNodeIds = vecGroup[gId - 1].get_vecNodeIds();
                        // cout << "checkpoint 1 .. " << endl;
                        vecNodeIds.erase(find(vecNodeIds.begin(), vecNodeIds.end(), nId));
                        // cout << "checkpoint 2" << endl ;
                        vecGroup[gId - 1].set_vecNodeIds(vecNodeIds);

                        // ... set the groupId to 26 "Z"
                        it->set_groupId(26);
                    }
                }
            }
        }

        // create a "cut" tumour if INIT_SPECIAL_PATTERN = "from_saved_keep_half"
        if (INIT_SPECIAL_PATTERN == "from_saved_keep_half")
        {
            cout << ">> keeping only the lower half " << endl;
            for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
            {
                bool isAlive = it->get_isAlive();
                if (isAlive)
                {
                    dVec coord = it->get_coord();

                    if (coord.y >= 0) // condition for cell removal to occur
                    {
                        // ... setting isAlive to false
                        it->set_isAlive(false);

                        // ... remove from vecNodeIds of vecGroup
                        int nId = it->get_id();
                        int gId = it->get_groupId();

                        cout << "... removing cell ID = " << nId << " in group ID = " << gId << endl;

                        // cout << "checkpoint 0 .. " << endl;
                        vector<int> vecNodeIds = vecGroup[gId - 1].get_vecNodeIds();
                        // cout << "checkpoint 1 .. " << endl;
                        vecNodeIds.erase(find(vecNodeIds.begin(), vecNodeIds.end(), nId));
                        // cout << "checkpoint 2" << endl ;
                        vecGroup[gId - 1].set_vecNodeIds(vecNodeIds);

                        // ... set the groupId to 26 "Z"
                        it->set_groupId(26);
                    }
                }
            }
        }
    }
}

// [NOT USED IN THE RMM STUDY] function initDynamicsBd() is a one-time function that initialize coordinates of nodeBd
void initDynamicsBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, vector<double> &randNum01)
{
    // Assumption 1 : a straight line of nodes
    if (TYPE_BOUNDARY == "vertical")
    {
        // const double POS_BOUNDARY_LINE_X = 200;
        int nCellPerLine = vecNodeBd.size();
        double dy = TUMOUR_RADIUS_MAX / (double)nCellPerLine;
        // double y_min = -CELL_RADIUS*nCellPerLine;
        double y_min = -0.5 * dy * nCellPerLine;
        vector<dVec> vecLineCoord;
        vecLineCoord.push_back({POS_BOUNDARY_LINE_X, y_min, 0});
        // createLineOfCoords(vecLineCoord, nCellPerLine, {0,1,0}, 2*CELL_RADIUS);
        createLineOfCoords(vecLineCoord, nCellPerLine, {0, 1, 0}, dy);

        int kk = 0;
        for (vector<Node>::iterator it = vecNodeBd.begin(); it != vecNodeBd.end(); it++)
        {
            dVec coord;
            if (it - vecNodeBd.begin() < nCellPerLine)
                coord = vecLineCoord[it - vecNodeBd.begin()];

            it->set_coord(coord);
            it->set_coord_prev(coord);

            // velocity
            it->set_veloc({0., 0., 0.});

            // time_double  --   these NodeBds don't divide
            double time_double = T + DT;
            it->set_time_double(time_double);
            it->set_time_double_intrinsic(time_double);

            // generation
            it->set_gen(0);

            // mechanics
            it->set_pressure(0.);
            it->set_radius(CELL_RADIUS); // assuming the same as cancer cells
        }
    }
    if (TYPE_BOUNDARY == "random" || TYPE_BOUNDARY == "radial" || TYPE_BOUNDARY == "circum" || TYPE_BOUNDARY == "radcir")
    {
        // need to set coordinates for all groupBd
        int k = 0;
        for (vector<Group>::iterator it_g = vecGroupBd.begin(); it_g != vecGroupBd.end(); it_g++)
        {
            vector<int> vecNodeIds = it_g->get_vecNodeIds();
            int nCellThisLine = vecNodeIds.size();
            double x_start, y_start, dx_ori, dy_ori;

            // randomize the starting point
            // x_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
            // y_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
            x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            while (x_start * x_start + y_start * y_start > TUMOUR_RADIUS_MAX * TUMOUR_RADIUS_MAX)
            {
                // x_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                // y_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
                y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            }
            // dx_ori = 2.* ((double) rand() / (RAND_MAX)) - 1;
            // dy_ori = 2.* ((double) rand() / (RAND_MAX)) - 1;

            // This is random orientation
            dx_ori = 2. * randNum01[k++] - 1;
            dy_ori = 2. * randNum01[k++] - 1;
            while (dx_ori * dx_ori + dy_ori * dy_ori > 1)
            {
                dx_ori = 2. * randNum01[k++] - 1;
                dy_ori = 2. * randNum01[k++] - 1;
            }
            double ori_size = sqrt(dx_ori * dx_ori + dy_ori * dy_ori);
            dx_ori /= ori_size;
            dy_ori /= ori_size;
            //
            // This is radial orientation
            if (TYPE_BOUNDARY == "radial")
            {
                double theta = atan2(y_start, x_start);
                dx_ori = cos(theta);
                dy_ori = sin(theta);
            }
            //
            // This is circum orientation
            if (TYPE_BOUNDARY == "circum")
            {
                double theta = atan2(y_start, x_start);
                dy_ori = cos(theta);
                dx_ori = -sin(theta);
            }
            //
            // This is radial (left) - circum (right)
            if (TYPE_BOUNDARY == "radcir")
            {
                double theta = atan2(y_start, x_start);
                dx_ori = cos(theta); // left  : radial
                dy_ori = sin(theta);
                if (x_start > 0) // right : circum
                {
                    dy_ori = cos(theta);
                    dx_ori = -sin(theta);
                }
            }

            vector<dVec> vecLineCoord;
            vecLineCoord.push_back({x_start, y_start, 0});
            createLineOfCoords(vecLineCoord, nCellThisLine, {dx_ori, dy_ori, 0}, 1.5 * CELL_RADIUS);
            for (vector<int>::iterator it = vecNodeIds.begin(); it != vecNodeIds.end(); it++)
            {
                dVec coord;
                if (it - vecNodeIds.begin() < nCellThisLine)
                    coord = vecLineCoord[it - vecNodeIds.begin()];

                int nodeId = *it;
                vecNodeBd[nodeId].set_coord(coord);
                vecNodeBd[nodeId].set_coord_prev(coord);

                // velocity
                vecNodeBd[nodeId].set_veloc({0., 0., 0.});

                // time_double  --   these NodeBds don't divide
                double time_double = T + DT;
                vecNodeBd[nodeId].set_time_double(time_double);
                vecNodeBd[nodeId].set_time_double_intrinsic(time_double);

                // generation
                vecNodeBd[nodeId].set_gen(0);

                // mechanics
                vecNodeBd[nodeId].set_pressure(0.);
                vecNodeBd[nodeId].set_radius(CELL_RADIUS); // assuming the same as cancer cells
            }
        }
    }

    if (TYPE_BOUNDARY == "ring")
    {
        const double sep = CELL_RADIUS;
        const int nCellThisRing = (int)(2 * PI * INIT_RING_DIST2CEN / sep);
        const int gap_size = INIT_RING_GAP_SIZE;
        const int bar_size = (int)(INIT_RING_GAP_SIZE * INIT_RING_BAR_GAP_RATIO);
        const int cycle = gap_size + bar_size;
        vecGroupBd.clear();
        vecNodeBd.clear();
        bool isAlive = true;
        int currGroupId = vecGroupBd.size() + 1;
        int currNodeId = vecNodeBd.size();
        int currLineageId = -1;
        double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;

        vector<int> vecNodeIds;
        for (int i = 0; i < nCellThisRing; i++)
        {
            if (i % cycle < gap_size) // belongs to gap
                continue;
            else
            {
                double ang = i * sep / INIT_RING_DIST2CEN;
                double x = INIT_RING_DIST2CEN * cos(ang), y = INIT_RING_DIST2CEN * sin(ang), z = 0;
                Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
                vecNodeIds.push_back(currNodeId);
                currNodeId++;
                vecNodeBd.push_back(node0);

                cout << i << " : " << x << ", " << y << endl;

                int nodeId = node0.get_id();
                dVec coord = {x, y, z};
                vecNodeBd[nodeId].set_coord(coord);
                vecNodeBd[nodeId].set_coord_prev(coord);

                // velocity
                vecNodeBd[nodeId].set_veloc({0., 0., 0.});

                // time_double  --   these NodeBds don't divide
                double time_double = T + DT;
                vecNodeBd[nodeId].set_time_double(time_double);
                vecNodeBd[nodeId].set_time_double_intrinsic(time_double);

                // generation
                vecNodeBd[nodeId].set_gen(0);

                // mechanics
                vecNodeBd[nodeId].set_pressure(0.);
                vecNodeBd[nodeId].set_radius(CELL_RADIUS); // assuming the same as cancer cells
            }
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
        vecGroupBd.push_back(group0);

        cout << "... init GroupBd ID = " << currGroupId << " ..." << endl;
        currGroupId++;
    }
}
// [NOT USED IN THE RMM STUDY] ONLY valid if typePattern == "two_half"
void initDynamicsBdSpecialPattern(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, string typePattern)
{
    // create a square to enclose the cell region
    if (typePattern == "two_half")
    {
        // (1) define a square region of cells that fit into the square region
        // (2) assign the right half as Group 2
        const int ncell_side = (int)sqrt(INIT_NUM_NODE_PER_GROUP);
        const double sep = 2. * CELL_RADIUS, sep_bd = 0.5 * CELL_RADIUS;
        dVec coord_lower_left, coord_upper_right;
        coord_lower_left = {-ncell_side * 0.5 * sep - sep, -ncell_side * 0.5 * sep * sqrt(3) * 0.5 - sep, 0};
        // coord_upper_right = {-ncell_side*0.5*sep  +(ncell_side+1)*sep, -ncell_side*0.5*sep*sqrt(3)*0.5 +(ncell_side+1)*sep, 0};
        coord_upper_right = {-ncell_side * 0.5 * sep + (ncell_side + 1) * sep, -ncell_side * 0.5 * sep * sqrt(3) * 0.5 + (ncell_side * 0.5 * sqrt(3) + 1) * sep, 0};

        vector<dVec> vecLineCoordW, vecLineCoordE, vecLineCoordN, vecLineCoordS;
        vecLineCoordW.push_back(coord_lower_left);
        vecLineCoordS.push_back(coord_lower_left);
        vecLineCoordN.push_back(coord_upper_right);
        vecLineCoordE.push_back(coord_upper_right);

        const int nCellPerLine = (int)((coord_upper_right.x - coord_lower_left.x) / sep_bd) + 1;
        const double sep_bd_2 = (coord_upper_right.y - coord_lower_left.y) / nCellPerLine;

        createLineOfCoords(vecLineCoordW, nCellPerLine, {0, 1, 0}, sep_bd_2);
        createLineOfCoords(vecLineCoordE, nCellPerLine, {0, -1, 0}, sep_bd_2);

        createLineOfCoords(vecLineCoordN, nCellPerLine, {-1, 0, 0}, sep_bd);
        createLineOfCoords(vecLineCoordS, nCellPerLine, {1, 0, 0}, sep_bd);

        vector<vector<dVec>> vecVecCoord{vecLineCoordW, vecLineCoordE, vecLineCoordN, vecLineCoordS};

        // create boundary nodes
        bool isAlive = true;
        int currGroupId = vecGroupBd.size() + 1;
        int currNodeId = vecNodeBd.size();
        int currLineageId = -1;
        double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif;

        vector<int> vecNodeIds;
        for (int iN = 0; iN < 4 * nCellPerLine; iN++)
        {
            Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
            int line = iN / nCellPerLine;
            int posi = iN % nCellPerLine;
            dVec coord0 = vecVecCoord[line][posi];
            node0.set_coord(coord0);
            node0.set_veloc({0, 0, 0});

            vecNodeBd.push_back(node0);
            vecNodeIds.push_back(currNodeId);
            currNodeId++;
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
        vecGroupBd.push_back(group0);

        cout << "... init GroupBd ID = " << currGroupId << " ..." << endl;
    }
}

// [NOT USED IN THE RMM STUDY] ONLY valid if STROMA_ON == true
void initStroma(vector<Group> &vecGroupSC, vector<Node> &vecNodeSC)
{
    // (1) define a square region of cell coordinates, with certain fraction removed
    // (2) collect remaining coordinates & construct the vecNodeSC & vecGroupSC
    cout << "--- initialising stroma space ---" << endl;
    const int ncell_side = INIT_SC_NCELL_SIDE;
    const int ncell_full = ncell_side * ncell_side;
    const double sep = 2. * CELL_RADIUS;
    dVec coord_lower_left;
    coord_lower_left = {-ncell_side * 0.5 * sep, -ncell_side * 0.5 * sep * sqrt(3) * 0.5, 0};
    vector<dVec> vec_coord_remain;

    cout << ".. selecting from " << ncell_full << " cell coords .." << endl;

    // STEP -- (1)
    for (int inode = 0; inode < ncell_full; inode++)
    {
        dVec coord;

        if (rand() / (float)RAND_MAX < INIT_SC_FRAC_KEEP)
        {
            cout << ".. add cell coord .." << endl;
            int row = inode / ncell_side;
            int col = inode % ncell_side;

            if (row % 2 == 0)
                coord = {coord_lower_left.x + col * sep, coord_lower_left.y + row * sep * sqrt(3) * 0.5, 0};
            else if (row % 2 == 1)
                coord = {coord_lower_left.x + (col + 0.5) * sep, coord_lower_left.y + row * sep * sqrt(3) * 0.5, 0};

            // avoid coords too close to centre (overlapping with cancer cells)
            if (abs(coord.x) < sep * 0.75 && abs(coord.y) < sep * 0.75)
                continue;

            // push back to vector
            vec_coord_remain.push_back(coord);
        }
    }

    // STEP -- (2)
    int ncell_remain = vec_coord_remain.size();
    bool isAlive = true;
    int currGroupId = vecGroupSC.size() + 1;
    int currNodeId = vecNodeSC.size();
    int currLineageId = 0;
    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1., scale_dif = 1.;

    vector<int> vecNodeIds;
    for (int iN = 0; iN < ncell_remain; iN++)
    {
        Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
        // set additional properties
        // coordinate
        dVec coord = vec_coord_remain[iN];
        node0.set_coord(coord);
        node0.set_coord_prev(coord);

        cout << iN << ": (" << coord.x << ", " << coord.y << ", " << coord.z << ")" << endl;

        // velocity
        node0.set_veloc({0., 0., 0.});

        // time_double
        double time_double = 10 * T; // make sure that growth is absent!
        node0.set_time_double(time_double);
        node0.set_time_double_intrinsic(time_double);

        // generation
        node0.set_gen(0);

        // mechanics
        node0.set_pressure(0.);
        node0.set_radius(CELL_RADIUS);
        node0.set_strain_energy(0.);

        vecNodeSC.push_back(node0);
        vecNodeIds.push_back(currNodeId);
        currNodeId++;
    }
    Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen, scale_dif);
    vecGroupSC.push_back(group0);
    currGroupId++;

    cout << " --- done ! --- " << endl;
}

// [NOT USED IN THE RMM STUDY] ONLY valid if BOUNDARY2_ON == true
void initBd2(vector<dVec> &vecBd2, vector<double> &randNum01)
{
    if (TYPE_BOUNDARY2 == "random" || TYPE_BOUNDARY2 == "radial" || TYPE_BOUNDARY2 == "circum" || TYPE_BOUNDARY2 == "radcir")
    {
        int k = 0;
        for (int j = 0; j < INIT_NUM_GROUPBD2; j++)
        {
            double x_start, y_start, dx_ori, dy_ori;

            // randomise starting point
            x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;

            while (x_start * x_start + y_start * y_start > TUMOUR_RADIUS_MAX * TUMOUR_RADIUS_MAX)
            {
                // x_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                // y_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
                y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            }

            if (TYPE_BOUNDARY == "random")
            {
                // This is random orientation
                dx_ori = 2. * randNum01[k++] - 1;
                dy_ori = 2. * randNum01[k++] - 1;
                // while (dx_ori*dx_ori+dy_ori*dy_ori > 1)
                //{
                //     dx_ori = 2.* randNum01[k++] - 1;
                //     dy_ori = 2.* randNum01[k++] - 1;
                // }
                double ori_size = sqrt(dx_ori * dx_ori + dy_ori * dy_ori);
                dx_ori /= ori_size;
                dy_ori /= ori_size;
            }

            // push_back all pixels to vecBd2
            vecBd2.push_back({x_start, y_start, 0});
            createLineOfCoords2(vecBd2, {dx_ori, dy_ori, 0});
        }

        cout << "... end" << endl;
    }
}

// [NOT USED IN THE RMM STUDY] the following functions
void currCloneCenterPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery, dVec &cloneCOM, bool flagCheckPeriphery)
{
    double xall = 0, yall = 0, zall = 0;
    int nNode = vecNode.size();
    vector<int> vecNodePeriphery_tmp;
    for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
    {
        dVec coord = it_n->get_coord();
        // dVec coord = it_n->get_coord_prev();
        xall += coord.x;
        yall += coord.y;
        zall += coord.z;

        // [NOT GOOD] method 1: decide periphery by nb number
        // vector<int> vecNbNodeIds = it_n->get_vecNbNodeIds();
        // if (vecNbNodeIds.size() <= PERIPHERY_MAX_NCB)   // warning: this might have an error if some interior cells have few "neighbors"!
        //    vecNodePeriphery.push_back(it_n->get_id());
        // method 2: decide periphery by umapCellDensity
        unordered_map<string, int> umapCellDensity = it_n->get_umapCellDensity();
        if (umapCellDensity["d1"] <= PERIPHERY_MAX_NCB)
            vecNodePeriphery_tmp.push_back(it_n->get_id());
    }
    cloneCOM = {xall / nNode, yall / nNode, zall / nNode};

    if (flagCheckPeriphery) // this removes nodes that are wrongly included
    {
        // cout << "number of candidate nodes at colony periphery : " << vecNodePeriphery_tmp.size() << endl;

        // testing by check inter-cell distances within periphery
        // ... problem: including many internal cells
        if (false)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                bool flagIsNodePeriphery = false;
                double dist_test = 3 * CELL_RADIUS;

                dVec coord = vecNode[*it_i].get_coord();
                // dVec coord = vecNode[*it_i].get_coord_prev();
                for (vector<int>::iterator it_j = vecNodePeriphery_tmp.begin(); it_j != vecNodePeriphery_tmp.end(); it_j++)
                {
                    if (*it_i != *it_j)
                    {
                        dVec coord2 = vecNode[*it_j].get_coord();
                        // dVec coord2 = vecNode[*it_j].get_coord_prev();
                        double dist = dVecDist(coord, coord2);
                        if (dist <= dist_test)
                        {
                            // cout << endl;
                            // cout << *it_i << " with dist = " << dist/CELL_RADIUS << endl;
                            // cout << *it_i << ": " << coord.x/CELL_RADIUS  << ", " << coord.y/CELL_RADIUS  << endl;
                            // cout << *it_j << ": " << coord2.x/CELL_RADIUS << ", " << coord2.y/CELL_RADIUS << endl;
                            flagIsNodePeriphery = true;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // [BEST SO FAR: testing by moving cell radially
        // ... problem: miss some cells near protrusive structures
        if (true)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                // dVec coord = vecNode[*it_i].get_coord_prev();
                double dist = dVecDist(coord, cloneCOM), factorTest = (3 * CELL_RADIUS + dist) / dist;
                dVec vecTest = {(coord.x - cloneCOM.x) * factorTest,
                                (coord.y - cloneCOM.y) * factorTest,
                                (coord.z - cloneCOM.z) * factorTest};
                dVec pointTest = {cloneCOM.x + vecTest.x, cloneCOM.y + vecTest.y, cloneCOM.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // testing by moving cell towards empty edge (decided according to instantaneous force)
        // ... problem: wrong when cell has tensile forces
        if (false)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                dVec veloc = vecNode[*it_i].get_veloc();
                double speed = sqrt(veloc.x * veloc.x + veloc.y * veloc.y + veloc.z * veloc.z);
                double factorTest = 3 * CELL_RADIUS;
                dVec vecTest = {
                    veloc.x / speed * factorTest, veloc.y / speed * factorTest, veloc.z / speed * factorTest};
                dVec pointTest = {
                    coord.x + vecTest.x, coord.y + vecTest.y, coord.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // testing by moving cell towards empty edge (decided according to neighbors)
        // ...
        if (true)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                vector<int> vecNbNodeIds = vecNode[*it_i].get_vecNbNodeIds();
                dVec sumVec = {0, 0, 0};
                for (vector<int>::iterator it_j = vecNbNodeIds.begin(); it_j != vecNbNodeIds.end(); it_j++)
                {
                    dVec coord2 = vecNode[*it_j].get_coord();
                    dVec sumVec_tmp = sumVec;
                    sumVec = {
                        sumVec_tmp.x + coord.x - coord2.x,
                        sumVec_tmp.y + coord.y - coord2.y,
                        sumVec_tmp.z + coord.z - coord2.z};
                }
                double sumVecSize = sqrt(sumVec.x * sumVec.x + sumVec.y * sumVec.y + sumVec.z * sumVec.z);

                double factorTest = 3.5 * CELL_RADIUS;
                dVec vecTest = {
                    sumVec.x / sumVecSize * factorTest, sumVec.y / sumVecSize * factorTest, sumVec.z / sumVecSize * factorTest};
                dVec pointTest = {
                    coord.x + vecTest.x, coord.y + vecTest.y, coord.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery &&
                    find(vecNodePeriphery.begin(), vecNodePeriphery.end(), *it_i) == vecNodePeriphery.end())
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // cout << "*CONFIRMED* number of candidate nodes at colony periphery : " << vecNodePeriphery.size() << endl;
    }
}
void sortClonePeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery)
{
    // function: sort according to angular locations
    unordered_map<double, int> mapAngleNode;
    vector<double> vecAngle;
    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        dVec coord = vecNode[*it_i].get_coord();
        double ang = atan2(coord.y, coord.x);
        while (find(vecAngle.begin(), vecAngle.end(), ang) != vecAngle.end())
            ang += 0.01 * PI;
        vecAngle.push_back(ang);
        mapAngleNode[ang] = *it_i;
    }

    // sort the vector of angles in ascending order
    sort(vecAngle.begin(), vecAngle.end());

    // clear vecNodePeriphery as information is transfered to mapAngleNode
    vecNodePeriphery.clear();
    for (vector<double>::iterator it_d = vecAngle.begin(); it_d != vecAngle.end(); it_d++)
    {
        vecNodePeriphery.push_back(mapAngleNode[*it_d]);
    }
}
void fillGapPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery)
{
    // function:
    // ... by checking (1) if two consecutive peripheral cells are NOT in contact and (2) if they have common neighbor cells
    // ... add their common neighbor to the vector

    vector<int> vecNodePeriphery_temp;

    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        int nId1 = *it_i, nId2;
        vecNodePeriphery_temp.push_back(nId1);

        // node2
        if (it_i - vecNodePeriphery.end() == -1)
            nId2 = vecNodePeriphery[0];
        else
            nId2 = *(it_i + 1);

        // find if nId2 is in nId1 neighbor list
        vector<int> vecNbNodeIds1 = vecNode[nId1].get_vecNbNodeIds();
        if (find(vecNbNodeIds1.begin(), vecNbNodeIds1.end(), nId2) != vecNbNodeIds1.end()) // found !
            continue;
        else // not found
        {
            // try to find if they share common neighbors
            vector<int> vecNbNodeIds2 = vecNode[nId2].get_vecNbNodeIds();
            for (vector<int>::iterator it_j = vecNbNodeIds1.begin(); it_j != vecNbNodeIds1.end(); it_j++)
            {
                if (find(vecNbNodeIds2.begin(), vecNbNodeIds2.end(), *it_j) != vecNbNodeIds2.end() &&
                    find(vecNodePeriphery.begin(), vecNodePeriphery.end(), *it_j) == vecNodePeriphery.end() &&
                    find(vecNodePeriphery_temp.begin(), vecNodePeriphery_temp.end(), *it_j) == vecNodePeriphery_temp.end())
                {
                    vecNodePeriphery_temp.push_back(*it_j);
                    break; // only include 1 common neighbor
                }
            }
        }
    }

    // cout << ">>Number of peripheral cells<<" << endl;
    // cout << "before filling gap: " << vecNodePeriphery.size() << endl;
    // cout << "after  filling gap: " << vecNodePeriphery_temp.size() << endl;
    // cout << "   net add:         " << vecNodePeriphery_temp.size() - vecNodePeriphery.size() << endl;

    vecNodePeriphery.clear();
    vecNodePeriphery = vecNodePeriphery_temp;
}
void calcCurvature(vector<Node> &vecNode, vector<int> &vecNodePeriphery, unordered_map<int, double> &mapNodeCurvature)
{
    // reference link: http://www.ambrsoft.com/trigocalc/circle3d.htm
    // A(x^2+y^2) + Bx + Cy + D = 0
    // A = x1(y2-y3) - y1(x2-x3) + x2y3 - x3y2
    // B = (x1^2+y1^2)(y3-y2) + (x2^2+y2^2)(y1-y3) + (x3^2+y3^2)(y2-y1)
    // C = (x1^2+y1^2)(x2-x3) + (x2^2+y2^2)(x3-x1) + (x3^2+y3^2)(x1-x2)
    // D = (x1^2+y1^2)(x3y2-x2y3) + (x2^2+y2^2)(x1y3-x3y1) + (x3^2+y3^2)(x2y1-x1y2)
    // r = sqrt(   (B^2+C^2-4AD)  /  4A^2    )

    double x1, y1, x2, y2, x3, y3;
    double A, B, C, D, r, kappa;
    bool flagIsConvex = true;
    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        x2 = vecNode[*it_i].get_coord().x;
        y2 = vecNode[*it_i].get_coord().y;

        // x1, y1
        if (it_i - vecNodePeriphery.begin() == 0)
        {
            x1 = vecNode[vecNodePeriphery[vecNodePeriphery.size() - 1]].get_coord().x;
            y1 = vecNode[vecNodePeriphery[vecNodePeriphery.size() - 1]].get_coord().y;
        }
        else
        {
            x1 = vecNode[*(it_i - 1)].get_coord().x;
            y1 = vecNode[*(it_i - 1)].get_coord().y;
        }
        // x3, y3
        if (it_i - vecNodePeriphery.end() == -1)
        {
            x3 = vecNode[vecNodePeriphery[0]].get_coord().x;
            y3 = vecNode[vecNodePeriphery[0]].get_coord().y;
        }
        else
        {
            x3 = vecNode[*(it_i + 1)].get_coord().x;
            y3 = vecNode[*(it_i + 1)].get_coord().y;
        }

        // concavity, convexity
        dPair vec21 = {x1 - x2, y1 - y2};
        dPair vec23 = {x3 - x2, y3 - y2};
        double ang_21_23 = atan2(vec21.x * vec23.y - vec21.y * vec23.x, vec21.x * vec23.x + vec21.y * vec23.y);
        if (ang_21_23 <= 0)
            flagIsConvex = true;
        else
            flagIsConvex = false;

        // A, B, C, D
        A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
        B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
        C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
        D = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);

        // r, kappa
        r = sqrt((B * B + C * C - 4 * A * D) / (4 * A * A));
        kappa = 1. / r;

        if (flagIsConvex == true)
        {
            // cout << "nId = " << *it_i << "; kappa = " << kappa << endl;
            mapNodeCurvature[*it_i] = kappa;
        }
        if (flagIsConvex == false)
        {
            // cout << "nId = " << *it_i << "; kappa = " << -kappa << endl;
            mapNodeCurvature[*it_i] = -kappa;
        }
    }
}

// function printGroupInfo() is an iterative function
double printGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, int seed)
{
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    int colony_size = 0;
    for (vector<Group>::iterator it_c = vecGroup.begin(); it_c != vecGroup.end(); it_c++)
    {
        vector<int> vecBeadIds = it_c->get_vecNodeIds();
        int groupId = it_c->get_id();
        cout << "Group_id: " << groupId << "; Group_size: " << vecBeadIds.size() << " cells" << endl;
        colony_size += vecBeadIds.size();

        int beadIdM;
        dVec coordM, velocM;
        double speedSqM = -1.;
        double speedSqAvg = 0;
        double speedAvg = 0;

        double pressureAvg = 0;
        double pressureMax = -1000;
        double pressureMin = 1000;
        int beadIdMaxP, beadIdMinP;

        for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
        {
            int beadId = *it_n;
            // check if the cell isAlive
            bool isAlive = vecNode[beadId].get_isAlive();
            if (isAlive)
            {
                dVec veloc = vecNode[beadId].get_veloc();
                dVec coord = vecNode[beadId].get_coord();
                double speedSq = veloc.x * veloc.x + veloc.y * veloc.y + veloc.z * veloc.z;
                // speedSqAvg += speedSq;
                speedAvg += sqrt(speedSq);
                if (speedSq > speedSqM)
                {
                    speedSqM = speedSq;
                    beadIdM = beadId;
                    coordM = coord;
                    velocM = veloc;
                }

                double pressure = vecNode[beadId].get_pressure();
                pressureAvg += pressure;
                if (pressure < pressureMin)
                {
                    pressureMin = pressure;
                    beadIdMinP = beadId;
                }
                if (pressure > pressureMax)
                {
                    pressureMax = pressure;
                    beadIdMaxP = beadId;
                }

                // if (pressure < 0)
                //     cout << beadId << " pressure : " << pressure << endl;

                if (coord.x < x_min)
                    x_min = coord.x;
                if (coord.x > x_max)
                    x_max = coord.x;
                if (coord.y < y_min)
                    y_min = coord.y;
                if (coord.y > y_max)
                    y_max = coord.y;
                if (coord.z < z_min)
                    z_min = coord.z;
                if (coord.z > z_max)
                    z_max = coord.z;
            }
        }
    }

    cout << "--- Summary of Virtual Tumour ---" << endl;
    cout << "DIM_X (um) \t DIM_Y (um) \t DIM_Z (um) \t NUM_CELL" << endl;
    cout << x_max - x_min << "\t" << y_max - y_min << "\t" << z_max - z_min << "\t" << colony_size << endl;

    // write tumour size into file
    ofstream tumourFile;
    string tumourFileName;
    stringstream ssTime;
    ssTime << fixed << setprecision(1) << t_now;
    stringstream PROC_ID_SS;
    PROC_ID_SS << PROC_ID;

    stringstream seed_SS;
    seed_SS << seed;

    tumourFileName = "PID_" + PROC_ID_SS.str() + "_seed_" + seed_SS.str() + "_tumour_size.txt";
    tumourFile.open(tumourFileName, ios::app | ios::binary);
    if (t_now == 0)
        tumourFile << "time\txdim\tydim\tncell" << endl; // NEED TO WRITE NUMBER OF APOPTOTIC CELLS!!!!
    stringstream ssxdim, ssydim, ssnc;
    ssxdim << fixed << setprecision(3) << x_max - x_min;
    ssydim << fixed << setprecision(3) << y_max - y_min;
    ssnc << fixed << setprecision(1) << vecNode.size();
    tumourFile << ssTime.str() << "\t" << ssxdim.str() << "\t" << ssydim.str() << "\t" << ssnc.str() << endl;
    tumourFile.close();

    return 0.5 * (x_max - x_min + y_max - y_min);
}

// function writeNodeDynamicsTherapy() is an iterative function
void writeNodeDynamicsTherapy(vector<Node> &vecNode, double t_now, int seed)
{
    ofstream nodeFile;
    string nodeFileName;
    stringstream PROC_ID_SS;
    PROC_ID_SS << PROC_ID;
    stringstream seed_SS;
    seed_SS << seed;

    nodeFileName = "PID_" + PROC_ID_SS.str() + "_seed_" + seed_SS.str() + "_nodeDynamicsTherapy.csv";

    nodeFile.open(nodeFileName, ios::app | ios::binary);
    if (t_now == 0)
        nodeFile << "t(h), CellID, x(um), y, z, DriverMutations" << endl;

    for (vector<Node>::iterator it_s = vecNode.begin(); it_s != vecNode.end(); it_s++)
    {
        // check if the cell isAlive
        bool isAlive = it_s->get_isAlive();
        if (isAlive)
        {
            // write node dynamics
            int nId = it_s->get_id();
            dVec coord = it_s->get_coord();
            vector<int> vecMutIds = it_s->get_vecMutIds();
            string vecMutIdsStr = "", vecAllMutIdsStr = "";
            for (vector<int>::iterator it_m = vecMutIds.begin(); it_m != vecMutIds.end(); it_m++)
            {
                if (*it_m < NUM_MUT_DRIVER)
                    vecMutIdsStr += to_string(*it_m) + ";";
                vecAllMutIdsStr += to_string(*it_m) + ";";
            }

            stringstream sst, ssid, ssx, ssy, ssz;
            sst << t_now;
            ssid << nId;
            ssx << coord.x;
            ssy << coord.y;
            ssz << coord.z;

            nodeFile << sst.str() << ", " << ssid.str() << ", "
                     << ssx.str() << ", " << ssy.str() << ", " << ssz.str() << ", "
                     << vecMutIdsStr << endl;
        }
    }

    nodeFile.close();
}

double dVecDist(dVec &p1, dVec &p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}
