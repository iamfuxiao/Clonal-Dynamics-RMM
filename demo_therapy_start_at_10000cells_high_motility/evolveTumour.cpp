/*
    File: evolveTumour.cpp
    Model: particleCell -- random mutation model (RMM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE RMM STUDY]
*/

#include "evolveTumour.hpp"

// one iteration of simulation
void oneIterOverdamp(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap, \
                     vector<double> &randNormNum0, vector<double> &randNormNumLife, vector<double> &randNormNumPthr, \
                     vector<double> &randNum01, vector<double> &randNum01b, \
                     vector<unordered_map<string,double>> &mitosisToWrite, \
                     vector<Node> &vecNodeBd, vector<Node> &vecNodeSC, \
                 vector<dVec> &vecBd2)
{
    // [1] proliferation & mitosis
    mitosis(vecGroup, vecNode, randNormNumLife, randNormNumPthr, randNum01, mitosisToWrite);

    // [1.2] proliferation-INdependent mutation processes
    if (TYPE_SUBCL == "emerge" && TYPE_EMERGE_SUBCL == "random2")
        intrinsicMutaProc(vecNode, randNum01b);

    // [2] apoptosis [NOT USED IN THE RMM STUDY]

    // [3] assign each Node to VoxMap by doing division and modulo
    updateVoxMap(voxMap, vecNode);

    // [4] update coordinate
    updateCoord(vecNode, vecGroup, randNormNum0, \
                vecNodeSC);

    // [5] update velocity
    updateVeloc(vecGroup, vecNode, voxMap, randNormNum0, \
                vecNodeBd, vecNodeSC, vecBd2);

}

// function updateCoord() is an iterative function
void updateCoord(vector<Node> & vecNode, vector<Group> & vecGroup, vector<double> &randNormNum0, \
                 vector<Node> & vecNodeSC)
{
    int k = 0;

    //dVec cOC = {0,0,0}, cOC_tmp = cOC;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        int beadId = it-vecNode.begin();
        if (beadId != it->get_id())
            cout << beadId << ", " << it->get_id() << endl;

        // check if the cell isAlive
        bool isAlive = it->get_isAlive();
        if (isAlive)
        {
            int groupId = it->get_groupId();
            dVec coord = it->get_coord(), veloc = it->get_veloc();
            dVec coord_new;
            double brownian_x = 0, brownian_y = 0, brownian_z = 0;

            if (BROWNIAN_ON == true)
            {
                double scale_vis = 1.;
                if (TYPE_DRIVER_MUTATION == "FRICTION")
                    scale_vis = vecGroup[groupId-1].get_scale_vis();

                brownian_x = randNormNum0[k++]*BROWNIAN_DIFF_COEF_SQRT/sqrt(scale_vis);  // um
                brownian_y = randNormNum0[k++]*BROWNIAN_DIFF_COEF_SQRT/sqrt(scale_vis);
                if (SIM_DIM == "3D")
                    brownian_z = randNormNum0[k++]*BROWNIAN_DIFF_COEF_SQRT/sqrt(scale_vis);

            }

            coord_new = { coord.x + veloc.x*DT + brownian_x,
                          coord.y + veloc.y*DT + brownian_y,
                          coord.z + veloc.z*DT + brownian_z};

            it->set_coord(coord_new);
            it->set_coord_prev(coord);

        }

    }

    // [NOT USED IN THE RMM STUDY] also update stroma cell coordinates, if STROMA_ON == true
    if (STROMA_ON && TYPE_SUBCL == "induce")
    {
        for (vector<Node>::iterator it_sc = vecNodeSC.begin(); it_sc != vecNodeSC.end(); it_sc++)
        {
            dVec coord_sc = it_sc->get_coord(), veloc_sc = it_sc->get_veloc();
            dVec coord_sc_new;
            coord_sc_new = { coord_sc.x + veloc_sc.x*DT,
                              coord_sc.y + veloc_sc.y*DT,
                              coord_sc.z + veloc_sc.z*DT};
            it_sc->set_coord(coord_sc_new);
            it_sc->set_coord_prev(coord_sc);
        }
    }
}

// function updateVeloc() is an iterative function
void updateVeloc(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap, vector<double> &randNormNum0, \
                 vector<Node> &vecNodeBd, vector<Node> &vecNodeSC, vector<dVec> &vecBd2)
{
    // [0] clear veloc values in Nodes
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        it->set_veloc({0, 0, 0});

        // also clean pressure
        it->set_pressure(0);
        vector<dPair> stress2d;
        stress2d = {{0,0}, {0,0}};
        it->set_stress2d(stress2d);

        // also clean compressive 'strain energy'
        it->set_strain_energy(0.);

        // also clean the vecNbNodeIds in Node
        vector<int> vecNbNodeIds;
        it->set_vecNbNodeIds(vecNbNodeIds);
        vector<double> vecNbNodeDists;
        it->set_vecNbNodeDists(vecNbNodeDists);
        // also clean the vecNbCommAreas in Node
        vector<double> vecNbCommAreas;
        it->set_vecNbCommAreas(vecNbCommAreas);
        // also clean the vecNbPressures in Node
        vector<double> vecNbPressures;
        it->set_vecNbPressures(vecNbPressures);
        // also clean the umapCellDensity in Node
        unordered_map<string, int> umapCellDensity;
        umapCellDensity["d1"] = 0; // 1 cell diameter
        umapCellDensity["d2"] = 0; // 2 cell diameters
        umapCellDensity["d3"] = 0; // 3 cell diameters
        it->set_umapCellDensity(umapCellDensity);
    }
    // [NOT USED IN THE RMM STUDY] ... also for stroma cells
    for (vector<Node>::iterator it_sc = vecNodeSC.begin(); it_sc != vecNodeSC.end(); it_sc++)
    {
        it_sc->set_veloc({0, 0, 0});
    }

    double gamma_inv = 1./GAMMA_VISCOUS;

    // [2] Intercell adhesion (& repulsion) force
    string flagIntercell_method = "ha"; // "ha" harmonic; "ha2" harmonic with variable cell size; "lj" lennard-jones
    
    // ... harmonic potential
    if (flagIntercell_method == "ha")
    {
        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds();
        unordered_map<int, vector<int>> mapBeadIdsHalo = voxMap.get_mapBeadIdsHalo();

        //cout << "-------------------------------------------------" << endl;
        //cout << "total number of cells : " << vecNode.size() << endl;
        //cout << "boxId\tnCellInBox\tnCellInHalo" << endl;
        int totalCellInVoxMap = 0;
        # pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) \
        shared(mapBeadIds, mapBeadIdsHalo, vecNode, totalCellInVoxMap)
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            vector<int> vecBeadIds0, vecBeadIds1;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIds[boxId].begin(); it != mapBeadIds[boxId].end(); it ++)
                    vecBeadIds0.push_back(*it);
            }
            unordered_map<int, vector<int>>::const_iterator got2 = mapBeadIdsHalo.find(boxId);
            if ( got2 != mapBeadIdsHalo.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIdsHalo[boxId].begin(); it != mapBeadIdsHalo[boxId].end(); it ++)
                    vecBeadIds1.push_back(*it);
            }

            int nBeadInBox = vecBeadIds0.size();
            int nBeadInHalo = vecBeadIds1.size();

            #pragma omp critical
            totalCellInVoxMap += nBeadInBox;

            if (nBeadInBox > 0)
            {
                //#pragma omp critical
                //cout << boxId << "\t" << nBeadInBox << "\t" << nBeadInHalo << endl;
                // STEP 1 : adhesion between Beads within the Box
                if (nBeadInBox > 1)
                {
                    for (int l = 0; l < nBeadInBox-1; l++)
                    {
                        // [NOT USED IN THE RMM STUDY] update veloc in cell according to interaction with boundary
                        // ... vecNodeBd is non-empty only if BOUNDARY_ON == true
                        for (vector<Node>::iterator it_n_bd = vecNodeBd.begin(); it_n_bd != vecNodeBd.end(); it_n_bd++)
                        {
                            int bId1 = vecBeadIds0[l];
                            bool isAlive1 = vecNode[bId1].get_isAlive();
                            if (isAlive1)
                            {
                                dVec coord1 = vecNode[bId1].get_coord(), coord_bd = it_n_bd->get_coord();
                                double Lrep = REPULSION_LENGTH_ACT;
                                double L = dVecDist(coord1, coord_bd);

                                if (L < Lrep && L > 0)
                                {
                                    // for the purpose of monitoring stress
                                    double a12 = PI*(CELL_RADIUS*CELL_RADIUS - 0.25*L*L);
                                    if (a12 < ADHESION_AREA_MIN)
                                        a12 = ADHESION_AREA_MIN;
                                    double f_rep_size = 0;

                                    dVec veloc1 = vecNode[bId1].get_veloc(), veloc1_new;
                                    dVec v_rep = {0,0,0};

                                    double L_inv = 1./L;
                                    dVec unitVec12;
                                    unitVec12 = {(coord_bd.x-coord1.x)*L_inv,
                                                 (coord_bd.y-coord1.y)*L_inv,
                                                 (coord_bd.z-coord1.z)*L_inv};
                                    int gId1 = vecNode[bId1].get_groupId();

                                    double v_rep_size = REPULSION_STRENGTH_BD*(L-Lrep)*gamma_inv;   // um/h
                                    f_rep_size = v_rep_size*GAMMA_VISCOUS/3.6/3.6;      // pN
                                    v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                                    // DRIVER MUTATION -- FRICTION
                                    double scale_vis_1 = 1.;
                                    if (TYPE_DRIVER_MUTATION == "FRICTION")
                                    {
                                        scale_vis_1 = vecGroup[gId1-1].get_scale_vis();
                                    }

                                    veloc1_new = {  veloc1.x + v_rep.x/scale_vis_1,
                                                    veloc1.y + v_rep.y/scale_vis_1,
                                                    veloc1.z + v_rep.z/scale_vis_1 };
                                    vecNode[bId1].set_veloc(veloc1_new);

                                    // update normal pressure due to cell-boundary contact
                                    //double pressure_tensile = (f_adh_size+f_rep_size)/a12;  // Pa
                                    double pressure_tensile = f_rep_size/a12;  // Pa
                                    double pressure1 = vecNode[bId1].get_pressure();
                                    pressure1 += pressure_tensile;

                                    vecNode[bId1].set_pressure(pressure1);

                                    // update the 'strain energy' of the volume exlusion spring
                                    double strain_energy = 0.5*REPULSION_STRENGTH_BD*(L-Lrep)*(L-Lrep); // pN*um
                                    double strain_energy1 = vecNode[bId1].get_strain_energy();
                                    strain_energy1 += 0.5*strain_energy;
                                    vecNode[bId1].set_strain_energy(strain_energy1);
                                }
                            }
                        }

                        // [NOT USED IN THE RMM STUDY] update veloc in cell according to interaction with stroma
                        if (STROMA_ON && TYPE_SUBCL == "induce")
                        {
                            for (vector<Node>::iterator it_n_sc = vecNodeSC.begin(); it_n_sc != vecNodeSC.end(); it_n_sc++)
                            {
                                int bId1 = vecBeadIds0[l];
                                bool isAlive1 = vecNode[bId1].get_isAlive();
                                if (isAlive1)
                                {
                                    dVec coord1 = vecNode[bId1].get_coord(), coord_sc = it_n_sc->get_coord();
                                    double Lrep = REPULSION_LENGTH_ACT;
                                    double L = dVecDist(coord1, coord_sc);

                                    if (L < Lrep && L > 0)
                                    {
                                        // for the purpose of monitoring stress
                                        double a12 = PI*(CELL_RADIUS*CELL_RADIUS - 0.25*L*L);
                                        if (a12 < ADHESION_AREA_MIN)
                                            a12 = ADHESION_AREA_MIN;
                                        double f_rep_size = 0;

                                        dVec veloc1 = vecNode[bId1].get_veloc(), veloc1_new;
                                        dVec v_rep = {0,0,0};

                                        double L_inv = 1./L;
                                        dVec unitVec12;
                                        unitVec12 = {(coord_sc.x-coord1.x)*L_inv,
                                                     (coord_sc.y-coord1.y)*L_inv,
                                                     (coord_sc.z-coord1.z)*L_inv};
                                        int gId1 = vecNode[bId1].get_groupId();

                                        double v_rep_size = REPULSION_STRENGTH_BD*(L-Lrep)*gamma_inv;   // um/h
                                        f_rep_size = v_rep_size*GAMMA_VISCOUS/3.6/3.6;      // pN
                                        v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                                        // DRIVER MUTATION -- FRICTION
                                        double scale_vis_1 = 1.;
                                        if (TYPE_DRIVER_MUTATION == "FRICTION")
                                        {
                                            scale_vis_1 = vecGroup[gId1-1].get_scale_vis();
                                        }

                                        veloc1_new = {  veloc1.x + v_rep.x/scale_vis_1,
                                                        veloc1.y + v_rep.y/scale_vis_1,
                                                        veloc1.z + v_rep.z/scale_vis_1 };
                                        vecNode[bId1].set_veloc(veloc1_new);

                                        // --- also update the velocity of stroma cells
                                        dVec veloc2_sc = it_n_sc->get_veloc(), veloc2_sc_new;
                                        veloc2_sc_new = {  veloc2_sc.x - v_rep.x,
                                                           veloc2_sc.y - v_rep.y,
                                                           veloc2_sc.z - v_rep.z };
                                        it_n_sc->set_veloc(veloc2_sc_new);
                                        // --------

                                        // update normal pressure due to cell-boundary contact
                                        //double pressure_tensile = (f_adh_size+f_rep_size)/a12;  // Pa
                                        double pressure_tensile = f_rep_size/a12;  // Pa
                                        double pressure1 = vecNode[bId1].get_pressure();
                                        pressure1 += pressure_tensile;

                                        vecNode[bId1].set_pressure(pressure1);

                                        // update the 'strain energy' of the volume exlusion spring
                                        double strain_energy = 0.5*REPULSION_STRENGTH_BD*(L-Lrep)*(L-Lrep); // pN*um
                                        double strain_energy1 = vecNode[bId1].get_strain_energy();
                                        strain_energy1 += 0.5*strain_energy;
                                        vecNode[bId1].set_strain_energy(strain_energy1);
                                    }
                                }
                            }

                        }

                        // [NOT USED IN THE RMM STUDY] check if near BOUNDARY2
                        bool isBd2 = false;
                        if (BOUNDARY2_ON)
                        {
                            for (vector<dVec>::iterator it_dvec = vecBd2.begin(); it_dvec != vecBd2.end(); it_dvec ++)
                            {
                                int bId1 = vecBeadIds0[l];
                                bool isAlive1 = vecNode[bId1].get_isAlive();
                                if (isAlive1)
                                {
                                    dVec coord1 = vecNode[bId1].get_coord(), coord_bd = *it_dvec;
                                    double L = dVecDist(coord1, coord_bd);
                                    if (L <= 0.5*INIT_BD2_THICKNESS)
                                    {
                                        isBd2 = true;
                                        break;
                                    }
                                }
                            }
                        }

                        for (int m = l+1; m < nBeadInBox; m++)
                        {
                            int bId1 = vecBeadIds0[l], bId2 = vecBeadIds0[m];
                            // check if both cells are alive
                            bool isAlive1 = vecNode[bId1].get_isAlive(), isAlive2 = vecNode[bId2].get_isAlive();
                            if (isAlive1 && isAlive2)
                            {
                                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                                double Lact = ADHESION_LENGTH_ACT, Leq = ADHESION_LENGTH_EQU, Lrep = REPULSION_LENGTH_ACT;   // um
                                double L = dVecDist(coord1, coord2);   // um

                                // update umapCellDensity in both cells
                                if (L <= 6*CELL_RADIUS)
                                {
                                    unordered_map<string, int> umapCD1, umapCD2;
                                    umapCD1 = vecNode[bId1].get_umapCellDensity();
                                    umapCD2 = vecNode[bId2].get_umapCellDensity();
                                    umapCD1["d3"] += 1;
                                    umapCD2["d3"] += 1;

                                    if (L <= 4*CELL_RADIUS)
                                    {
                                        umapCD1["d2"] += 1;
                                        umapCD2["d2"] += 1;

                                        if (L <= 2*CELL_RADIUS)
                                        {
                                            umapCD1["d1"] += 1;
                                            umapCD2["d1"] += 1;
                                        }
                                    }

                                    vecNode[bId1].set_umapCellDensity(umapCD1);
                                    vecNode[bId2].set_umapCellDensity(umapCD2);
                                }

                                // update veloc in both cells
                                if (L < Lact && L > 0)
                                {
                                    // for the purpose of monitoring stress
                                    double a12 = PI*(CELL_RADIUS*CELL_RADIUS - 0.25*L*L);
                                    if (a12 < ADHESION_AREA_MIN)
                                        a12 = ADHESION_AREA_MIN;
                                    double f_adh_size = 0, f_rep_size = 0;
                                    double strain_energy = 0;

                                    dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
                                    dVec v_adh = {0,0,0}, v_rep = {0,0,0};

                                    double L_inv = 1./L;
                                    dVec unitVec12;
                                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};

                                    int gId1 = vecNode[bId1].get_groupId();
                                    int gId2 = vecNode[bId2].get_groupId();
                                    if (L > Leq)
                                    {
                                        //double v_adh_size = ADHESION_STRENGTH*(L-Leq)*gamma_inv*1E-6;   // um/h
                                        double v_adh_size = ADHESION_STRENGTH*(L-Leq)*gamma_inv;   // um/h
                                        f_adh_size = v_adh_size*GAMMA_VISCOUS/3.6/3.6;      // pN

                                        // DRIVER MUTATION -- ADHESION
                                        if (TYPE_DRIVER_MUTATION == "ADHESION")
                                        {
                                            double scale_adh_1 = vecGroup[gId1-1].get_scale_adh();
                                            double scale_adh_2 = vecGroup[gId2-1].get_scale_adh();
                                            v_adh_size *= min(scale_adh_1, scale_adh_2);
                                            //#pragma omp critical
                                            //{
                                            //    cout << "groupId1 = " << gId1 << "; scale_adh_1 = " << scale_adh_1 << endl;
                                            //    cout << "groupId2 = " << gId2 << "; scale_adh_2 = " << scale_adh_2 << endl;
                                            //    cout << min(scale_adh_1, scale_adh_2) << endl;
                                            //}
                                        }
                                        v_adh = {v_adh_size*unitVec12.x, v_adh_size*unitVec12.y, v_adh_size*unitVec12.z};

                                    }
                                    if (L < Lrep)
                                    {
                                        double v_rep_size = REPULSION_STRENGTH*(L-Lrep)*gamma_inv;   // um/h
                                        f_rep_size = v_rep_size*GAMMA_VISCOUS/3.6/3.6;      // pN
                                        v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                                        // update the 'strain energy' of the volume exlusion spring
                                        strain_energy = 0.5*REPULSION_STRENGTH*(L-Lrep)*(L-Lrep); // pN*um
                                    }

                                    // DRIVER MUTATION -- FRICTION
                                    double scale_vis_1 = 1., scale_vis_2 = 1.;
                                    if (TYPE_DRIVER_MUTATION == "FRICTION")
                                    {
                                        scale_vis_1 = vecGroup[gId1-1].get_scale_vis();
                                        scale_vis_2 = vecGroup[gId2-1].get_scale_vis();
                                    }

                                    // [NOT USED IN THE RMM STUDY]BOUNDARY2 scale drag coefficient
                                    if (BOUNDARY2_ON && isBd2)
                                    {
                                        scale_vis_1 *= BD2_SCALE_VIS_BY;
                                        scale_vis_2 *= BD2_SCALE_VIS_BY;
                                    }

                                    veloc1_new = {  veloc1.x + v_adh.x/scale_vis_1 + v_rep.x/scale_vis_1,
                                                    veloc1.y + v_adh.y/scale_vis_1 + v_rep.y/scale_vis_1,
                                                    veloc1.z + v_adh.z/scale_vis_1 + v_rep.z/scale_vis_1 };
                                    vecNode[bId1].set_veloc(veloc1_new);

                                    veloc2_new = {  veloc2.x - v_adh.x/scale_vis_2 - v_rep.x/scale_vis_2,
                                                    veloc2.y - v_adh.y/scale_vis_2 - v_rep.y/scale_vis_2,
                                                    veloc2.z - v_adh.z/scale_vis_2 - v_rep.z/scale_vis_2 };
                                    vecNode[bId2].set_veloc(veloc2_new);

                                    // push_back each other into vecNbNodeIds
                                    vector<int> vecNbNodeIds1 = vecNode[bId1].get_vecNbNodeIds();
                                    vecNbNodeIds1.push_back(bId2);
                                    vector<int> vecNbNodeIds2 = vecNode[bId2].get_vecNbNodeIds();
                                    vecNbNodeIds2.push_back(bId1);

                                    vecNode[bId1].set_vecNbNodeIds(vecNbNodeIds1);
                                    vecNode[bId2].set_vecNbNodeIds(vecNbNodeIds2);

                                    vector<double> vecNbNodeDists1 = vecNode[bId1].get_vecNbNodeDists();
                                    vecNbNodeDists1.push_back(L);
                                    vecNode[bId1].set_vecNbNodeDists(vecNbNodeDists1);
                                    vector<double> vecNbNodeDists2 = vecNode[bId2].get_vecNbNodeDists();
                                    vecNbNodeDists2.push_back(L);
                                    vecNode[bId2].set_vecNbNodeDists(vecNbNodeDists2);

                                    // update normal pressure due to cell-cell contact
                                    double pressure_tensile = f_rep_size/a12;  // Pa
                                    double pressure1 = vecNode[bId1].get_pressure();
                                    pressure1 += pressure_tensile;
                                    double pressure2 = vecNode[bId2].get_pressure();
                                    pressure2 += pressure_tensile;
                                    vector<double> vecNbPressures1 = vecNode[bId1].get_vecNbPressures();
                                    vector<double> vecNbPressures2 = vecNode[bId2].get_vecNbPressures();
                                    vecNbPressures1.push_back(pressure_tensile);
                                    vecNbPressures2.push_back(pressure_tensile);

                                    vecNode[bId1].set_pressure(pressure1);
                                    vecNode[bId2].set_pressure(pressure2);
                                    vecNode[bId1].set_vecNbPressures(vecNbPressures1);
                                    vecNode[bId2].set_vecNbPressures(vecNbPressures2);
                                    // --------


                                    // update strain energy in the volume exclusion part
                                    double strain_energy1 = vecNode[bId1].get_strain_energy();
                                    strain_energy1 += 0.5*strain_energy;
                                    vecNode[bId1].set_strain_energy(strain_energy1);
                                    double strain_energy2 = vecNode[bId2].get_strain_energy();
                                    strain_energy2 += 0.5*strain_energy;
                                    vecNode[bId2].set_strain_energy(strain_energy2);

                                    // --------

                                    dVec force_total = { (f_adh_size+f_rep_size)*unitVec12.x,
                                                         (f_adh_size+f_rep_size)*unitVec12.y,
                                                         (f_adh_size+f_rep_size)*unitVec12.z };
                                    dVec dist_vec = { L*unitVec12.x, L*unitVec12.y, L*unitVec12.z};
                                    vector<dPair> stress_2d = vecNode[bId1].get_stress2d();
                                    dPair stress_2d_xj = {stress_2d[0].x+dist_vec.x*force_total.x,
                                                           stress_2d[0].y+dist_vec.x*force_total.y};
                                    dPair stress_2d_yj = {stress_2d[1].x+dist_vec.y*force_total.x,
                                                           stress_2d[1].y+dist_vec.y*force_total.y};
                                    vector<dPair> stress_2d_new = {stress_2d_xj, stress_2d_yj};


                                    // --------
                                    // push_back each other into vecNbCommAreas
                                    vector<double> vecNbCommAreas1 = vecNode[bId1].get_vecNbCommAreas();
                                    vecNbCommAreas1.push_back(a12);
                                    vector<double> vecNbCommAreas2 = vecNode[bId2].get_vecNbCommAreas();
                                    vecNbCommAreas2.push_back(a12);

                                    vecNode[bId1].set_vecNbCommAreas(vecNbCommAreas1);
                                    vecNode[bId2].set_vecNbCommAreas(vecNbCommAreas2);
                                }
                            }

                        }
                    }
                }

                // STEP 2 : repulsion between one Bead within the Box and another Bead within Halo
                if (nBeadInBox > 0 && nBeadInHalo > 0)
                {
                    for (int l = 0; l < nBeadInBox; l++)
                    //for (int l = 0; l < nBeadInBox-1; l++) // old
                    {
                        for (int m = 0; m < nBeadInHalo; m++)
                        {
                            int bId1 = vecBeadIds0[l], bId2 = vecBeadIds1[m];
                            // check if both cells are alive
                            bool isAlive1 = vecNode[bId1].get_isAlive(), isAlive2 = vecNode[bId2].get_isAlive();
                            if (isAlive1 && isAlive2)
                            {
                                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                                double Lact = ADHESION_LENGTH_ACT, Leq = ADHESION_LENGTH_EQU, Lrep = REPULSION_LENGTH_ACT;   // um
                                double L = dVecDist(coord1, coord2);   // um

                                // update umapCellDensity in both cells
                                if (L <= 6*CELL_RADIUS)
                                {
                                    unordered_map<string, int> umapCD1;
                                    umapCD1 = vecNode[bId1].get_umapCellDensity();
                                    umapCD1["d3"] += 1;

                                    if (L <= 4*CELL_RADIUS)
                                    {
                                        umapCD1["d2"] += 1;

                                        if (L <= 2*CELL_RADIUS)
                                        {
                                            umapCD1["d1"] += 1;
                                            //#pragma omp critical
                                            //cout << bId1 << ", " << bId2 << endl;
                                        }
                                    }

                                    vecNode[bId1].set_umapCellDensity(umapCD1);
                                }

                                if (L < Lact && L > 0)
                                {
                                    // for the purpose of monitoring stress
                                    double a12 = PI*(CELL_RADIUS*CELL_RADIUS - 0.25*L*L);
                                    if (a12 < ADHESION_AREA_MIN)
                                        a12 = ADHESION_AREA_MIN;
                                    double f_adh_size = 0, f_rep_size = 0;
                                    double strain_energy = 0;


                                    dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
                                    dVec v_adh = {0,0,0}, v_rep = {0,0,0};

                                    double L_inv = 1./L;
                                    dVec unitVec12;
                                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};

                                    int gId1 = vecNode[bId1].get_groupId();
                                    int gId2 = vecNode[bId2].get_groupId();
                                    if (L > Leq)
                                    {
                                        double v_adh_size = ADHESION_STRENGTH*(L-Leq)*gamma_inv;   // um/h
                                        f_adh_size = v_adh_size*GAMMA_VISCOUS/3.6/3.6;      // pN
                                        v_adh = {v_adh_size*unitVec12.x, v_adh_size*unitVec12.y, v_adh_size*unitVec12.z};

                                    }
                                    if (L < Lrep)
                                    {
                                        //double v_rep_size = REPULSION_STRENGTH*(L-Leq)*gamma_inv;   // um/h
                                        double v_rep_size = REPULSION_STRENGTH*(L-Lrep)*gamma_inv;   // um/h
                                        f_rep_size = v_rep_size*GAMMA_VISCOUS/3.6/3.6;      // pN
                                        v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                                        // update the 'strain energy' of the volume exlusion spring
                                        strain_energy = 0.5*REPULSION_STRENGTH*(L-Lrep)*(L-Lrep); // pN*um

                                    }

                                    // DRIVER MUTATION -- FRICTION
                                    double scale_vis_1 = 1.;
                                    if (TYPE_DRIVER_MUTATION == "FRICTION")
                                    {
                                        scale_vis_1 = vecGroup[gId1-1].get_scale_vis();
                                    }

                                    // ONLY update bead in the VoxBox, not in the VoxHalo
                                    veloc1_new = {  veloc1.x + v_adh.x/scale_vis_1 + v_rep.x/scale_vis_1,
                                                    veloc1.y + v_adh.y/scale_vis_1 + v_rep.y/scale_vis_1,
                                                    veloc1.z + v_adh.z/scale_vis_1 + v_rep.z/scale_vis_1 };
                                    vecNode[bId1].set_veloc(veloc1_new);

//#pragma omp critical
                                    /*
                                    {
                                        veloc2_new = {  veloc2.x - v_adh.x - v_rep.x,
                                            veloc2.y - v_adh.y - v_rep.y,
                                            veloc2.z - v_adh.z - v_rep.z };
                                        vecNode[bId2].set_veloc(veloc2_new);
                                    }
                                     */

                                    // ONLY update bead in the VoxBox, not in the VoxHalo
                                    vector<int> vecNbNodeIds1 = vecNode[bId1].get_vecNbNodeIds();
                                    vecNbNodeIds1.push_back(bId2);
                                    vecNode[bId1].set_vecNbNodeIds(vecNbNodeIds1);
                                    vector<double> vecNbNodeDists1 = vecNode[bId1].get_vecNbNodeDists();
                                    vecNbNodeDists1.push_back(L);
                                    vecNode[bId1].set_vecNbNodeDists(vecNbNodeDists1);
//#pragma omp critical
                                    /*
                                    {
                                        vector<int> vecNbNodeIds2 = vecNode[bId2].get_vecNbNodeIds();
                                        vecNbNodeIds2.push_back(bId1);
                                        vecNode[bId2].set_vecNbNodeIds(vecNbNodeIds2);
                                        vector<double> vecNbNodeDists2 = vecNode[bId2].get_vecNbNodeDists();
                                        vecNbNodeDists2.push_back(L);
                                        vecNode[bId2].set_vecNbNodeDists(vecNbNodeDists2);
                                    }
                                    */

                                    // update normal pressure due to cell-cell contact
                                    //double pressure_tensile = (f_adh_size+f_rep_size)/a12;  // Pa
                                    double pressure_tensile = f_rep_size/a12;  // Pa
                                    double pressure1 = vecNode[bId1].get_pressure();
                                    pressure1 += pressure_tensile;
                                    vecNode[bId1].set_pressure(pressure1);

                                    vector<double> vecNbPressures1 = vecNode[bId1].get_vecNbPressures();
                                    vecNbPressures1.push_back(pressure_tensile);
                                    vecNode[bId1].set_vecNbPressures(vecNbPressures1);
                                    // --------

                                    // update strain energy in the volume exclusion part
                                    double strain_energy1 = vecNode[bId1].get_strain_energy();
                                    strain_energy1 += 0.5*strain_energy;
                                    vecNode[bId1].set_strain_energy(strain_energy1);
                                    // --------

                                    // push_back each other into vecNbCommAreas
                                    vector<double> vecNbCommAreas1 = vecNode[bId1].get_vecNbCommAreas();
                                    vecNbCommAreas1.push_back(a12);
                                    vecNode[bId1].set_vecNbCommAreas(vecNbCommAreas1);
                                }
                            }
                        }
                    }
                }

            }

        }

        if (totalCellInVoxMap != vecNode.size())
            cout << "WARNING: discrepancy between VoxMap and vecNode! " << totalCellInVoxMap << ", " << vecNode.size() << endl;

    }

}

// function mitosis() is an iterative function
void mitosis(vector<Group> &vecGroup, vector<Node> &vecNode, \
    vector<double> &randNormNumLife, vector<double> &randNormNumPthr, \
    vector<double> &randNum01, \
    vector<unordered_map<string,double>> &mitosisToWrite)
{
    vector<int> vecNodeIdsDivide;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        // !! check if the cell isAlive
        bool isAlive = it->get_isAlive();
        if (isAlive)
        {            // account for CIP: get vecNbNodeIds.size()
            bool flagCIP = false;
            if (MITOSIS_CIP_TYPE == "Mechanics")
            {
                // method 1: sum of all contacts
                int nId = it->get_id();
                int gId = vecNode[nId].get_groupId();
                double scale_mit = 1.;
                
                double stress = it->get_pressure();
                // scale the threshold in order to implement driver mutation with type "MITOSIS_PTHR"
                if (stress < scale_mit*MITOSIS_CIP_THRP)  // negative means compressive stress
                {
                    flagCIP = true;
                }

            }

            // decrease remaining time to mitosis
            double time_double = it->get_time_double();
            double t_minus = DT;
            
            time_double -= t_minus;
            it->set_time_double(time_double);

            // push_back node ids that will undergo division
            if (time_double <= 0 && flagCIP == false)
            {
                vecNodeIdsDivide.push_back(it-vecNode.begin());
            }

        }
    }

    // record the positions of Nodes when undergoing mitosis
    vector<int> vecNodePeriphery;
    dVec cloneCOM = {0,0,0};
    int nNode = vecNode.size();
    double cloneRad = sqrt(nNode)/2;   // measured in number of cells (R = sqrt(N)/2 * 2r)
    if (vecNodeIdsDivide.size())
    {
        bool flagCheckPeriphery = false;
        currCloneCenterPeriphery(vecNode, vecNodePeriphery, cloneCOM, flagCheckPeriphery);
    }

    // create a new cell which copies certain information from its parent cell
    int kk = 0, ll = 0, mm = 0;
    for (vector<int>::iterator it_i = vecNodeIdsDivide.begin(); it_i != vecNodeIdsDivide.end(); it_i++)
    {
        int nId = *it_i;
        int lId = vecNode[nId].get_lineageId();
        int gId = vecNode[nId].get_groupId();
        int gen = vecNode[nId].get_gen(), gen_new = gen+1;   // generation adds 1 for both cells
        int nId2 = vecNode.size();  
        vector<int> vecNbNodeIds = vecNode[nId].get_vecNbNodeIds();
        const int nNb = vecNbNodeIds.size();

        // create new Node
        bool isAlive = true;
        dVec coord_old = vecNode[nId].get_coord(), coord2 = coord_old, coord = coord_old;
        dVec veloc = vecNode[nId].get_veloc();
        Node n2 = Node(isAlive, nId2, lId, gId);
        vecNode.push_back(n2);

        // [this is to collect mitosis information] calculate distance of dividing cell to clone com and periphery
        double d2center = 0, d2periphery = cloneRad*1.5;
        d2center = dVecDist(coord_old, cloneCOM);
        
        // ... also collect the intrinsic and actual doubling time
        double time_double_intrinsic = vecNode[nId].get_time_double_intrinsic();
        double time_double_actual = time_double_intrinsic-vecNode[nId].get_time_double();  // note that the get_time_double() should return a negative value if the mitosis is delayed
        
        unordered_map<string, double> mapInfo;
        mapInfo["cloneSize"] = nNode;
        mapInfo["cloneRad"]  = cloneRad;
        mapInfo["lineageId"] = lId;
        mapInfo["nodeId"]    = nId;
        mapInfo["nodeId_child"] = nId2;
        mapInfo["d2CloneCenter"] = d2center/CELL_RADIUS/2;
        mapInfo["d2ClonePeriphery"] = d2periphery/CELL_RADIUS/2;
        mapInfo["time_double_intrinsic"] = time_double_intrinsic;
        mapInfo["time_double_actual"]    = time_double_actual;
        mapInfo["generation"] = gen;
        mitosisToWrite.push_back(mapInfo);

        // ... shift coordinate 
        double phi   = 2.* PI * randNum01[ll++], theta;

        coord  = { coord_old.x + CELL_RADIUS*0.25*cos(phi),
                   coord_old.y + CELL_RADIUS*0.25*sin(phi),
                   coord_old.z };
        coord2 = { coord_old.x - CELL_RADIUS*0.25*cos(phi),
                   coord_old.y - CELL_RADIUS*0.25*sin(phi),
                   coord_old.z };
        if (SIM_DIM == "3D")
        {
            theta = acos(2.*randNum01[ll++]-1);
            coord  = { coord_old.x + CELL_RADIUS*0.25*sin(theta)*cos(phi),
                        coord_old.y + CELL_RADIUS*0.25*sin(theta)*sin(phi),
                        coord_old.z + CELL_RADIUS*0.25*cos(theta)};
            coord2 = { coord_old.x - CELL_RADIUS*0.25*sin(theta)*cos(phi),
                       coord_old.y - CELL_RADIUS*0.25*sin(theta)*sin(phi),
                       coord_old.z - CELL_RADIUS*0.25*cos(theta)};
        }
        vecNode[nId].set_coord(coord);
        vecNode[nId2].set_coord(coord2);
        vecNode[nId].set_veloc(veloc);
        vecNode[nId2].set_veloc(veloc);

        // update time_double
        double scale_mit = 1.;
        double time_double_new = scale_mit*randNormNumLife[kk++];
        double time_double_new_2 = scale_mit*randNormNumLife[kk++];

        vecNode[nId].set_time_double( time_double_new );
        vecNode[nId].set_time_double_intrinsic( time_double_new );
        vecNode[nId2].set_time_double(time_double_new_2);
        vecNode[nId2].set_time_double_intrinsic(time_double_new_2);

        // update pthr
        vecNode[nId].set_pthr(randNormNumPthr[mm++]);
        vecNode[nId2].set_pthr(randNormNumPthr[mm++]);

        // update generation
        vecNode[nId].set_gen(gen_new);
        vecNode[nId2].set_gen(gen_new);

        // radius
        vecNode[nId2].set_radius(vecNode[nId].get_radius());

        // update Group vecNodeIds
        vector<int> vecNodeIds = vecGroup[gId-1].get_vecNodeIds();
        vecNodeIds.push_back(nId2);
        vecGroup[gId-1].set_vecNodeIds(vecNodeIds);

        // accumulate mutations (SINCE 2019.03.18)
        // NOTE : this is an NEWER implementation of "random" emerge model -- proliferation dependent
        if (TYPE_SUBCL == "emerge" && TYPE_EMERGE_SUBCL == "random") // emerge via mutation -- proliferation dependent
        {
            vector<int> vecMutIds = vecNode[nId].get_vecMutIds();
            for (int cnt_m = 0; cnt_m < NUM_MUT_PER_DIV; cnt_m ++)
            {
                int mut_id = static_cast<int>(std::floor(randNum01[ll++] * NUM_MUT_POOL));

                if (find(vecMutIds.begin(), vecMutIds.end(), mut_id) == vecMutIds.end())    
                    vecMutIds.push_back(mut_id);
            }
            vecNode[nId].set_vecMutIds(vecMutIds);
            vecNode[nId2].set_vecMutIds(vecMutIds);
        }
        //
        if (TYPE_SUBCL == "emerge" && TYPE_EMERGE_SUBCL == "random2") // emerge via mutation -- proliferation independent
        {
            vector<int> vecMutIds = vecNode[nId].get_vecMutIds();
            // not acquiring mutations here; but randomly over time in another function!
            vecNode[nId2].set_vecMutIds(vecMutIds);
        }
    }
}

// [NOT USED IN THE RMM STUDY] function apoptosis() models the spontaneous cell death
void apoptosis(vector<Group> &vecGroup, vector<Node> &vecNode)
{
    // Assumption 1: cells have a constant spontaneous apoptotic rate
    // ... death means setting isAlive to false
    //const double apop_prob = APOPTOSIS_RATE * DT;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        bool isAlive = it->get_isAlive();
        if (isAlive)
        {
            double apop_rand = random()/(float)RAND_MAX;
            double apop_prob = APOPTOSIS_RATE * DT;

            // deleterious mutations could result in enhanced apoptotic rate
            if (TYPE_SUBCL == "emerge" && TYPE_DRIVER_MUTATION == "APOPTOSIS")
            {
                vector<int> vecMutIds = it->get_vecMutIds();
                bool flagCellIsDriver = false;
                //cout << NUM_MUT_DRIVER << endl;
                for (vector<int>::iterator it_m = vecMutIds.begin(); it_m != vecMutIds.end(); it_m ++)
                {
                    if (*it_m < NUM_MUT_DRIVER)
                    {
                        flagCellIsDriver = true;
                        break;
                    }
                }

                if (flagCellIsDriver)
                {
                    double scale_apop = GROUP_SCALE_APO;
                    apop_prob *= scale_apop;
                    //cout << "cell id = " << nId << " harbors driver mutation " << endl;
                }
            }

            if (apop_rand < apop_prob)  // condition for cell death to occur
            {
                // ... setting isAlive to false
                it->set_isAlive(false);

                // ... remove from vecNodeIds of vecGroup
                int nId = it->get_id();
                int gId = it->get_groupId();
                vector<int> vecNodeIds = vecGroup[gId-1].get_vecNodeIds();
                vecNodeIds.erase(find(vecNodeIds.begin(), vecNodeIds.end(), nId));
                vecGroup[gId-1].set_vecNodeIds(vecNodeIds);

                // ... set the groupId to 26 "Z"
                it->set_groupId(26);
            }

        }
    }
}

// function necrosis2() models the cell's reponse to chemotherapy
// if TYPE_DRUG_SCHEDULE == "continuous"
// ... also checked if TYPE_DRUG_RESPONSE == "binary" at the higher level
void necrosis2(vector<Group> &vecGroup, vector<Node> &vecNode)
{
    // Assumption: cells have a constant necrosis rate under continuous therapy
    // ... death means setting isAlive to false
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        bool isAlive = it->get_isAlive();
        if (isAlive)
        {
            //double apop_rand = random()/(float)RAND_MAX;
            //double apop_prob = APOPTOSIS_RATE * DT;
            double necr_rand = random()/(float)RAND_MAX;
            double necr_prob = DRUG_EFFECT_RATE * DT;

            // drug-resistant mutations could result in full protection from therapy (zero apoptotic rate)
            if (TYPE_SUBCL == "emerge" && TYPE_DRUG_SCHEDULE == "continuous")
            {
                vector<int> vecMutIds = it->get_vecMutIds();
                bool flagCellIsDriver = false;
                //cout << NUM_MUT_DRIVER << endl;
                for (vector<int>::iterator it_m = vecMutIds.begin(); it_m != vecMutIds.end(); it_m ++)
                {
                    if (*it_m < NUM_MUT_DRIVER)
                    {
                        flagCellIsDriver = true;
                        break;
                    }
                }

                if (flagCellIsDriver)
                {
                    double scale_sen = GROUP_SCALE_SEN;
                    necr_prob *= scale_sen;
                }

                // [NOT USED IN THE RMM STUDY] if TYPE_DRUG == "cytostatic", further check cell's proliferative status
                if (TYPE_DRUG == "cytostatic")
                {
                    bool flagCIP = false;
                    int nId = it->get_id();
                    double stress = it->get_pressure();
                    if (stress < MITOSIS_CIP_THRP)  // negative means compressive stress
                        flagCIP = true;
                    if (flagCIP) // even sensitive cells are protected from cytostatic drug if they are restricted in proliferation
                        necr_prob = 0;
                }
            }

            if (necr_rand < necr_prob)  // condition for cell death to occur
            {
                // ... setting isAlive to false
                it->set_isAlive(false);

                // ... remove from vecNodeIds of vecGroup
                int nId = it->get_id();
                int gId = it->get_groupId();
                vector<int> vecNodeIds = vecGroup[gId-1].get_vecNodeIds();
                vecNodeIds.erase(find(vecNodeIds.begin(), vecNodeIds.end(), nId));
                vecGroup[gId-1].set_vecNodeIds(vecNodeIds);

                // ... set the groupId to 26 "Z"
                it->set_groupId(26);
            }

        }
    }
}

// ... introduce random mutations that are proliferation-independent "random2" emerge model
void intrinsicMutaProc(vector<Node> &vecNode, vector<double> &randNum01b)
{
    int ll = 0;
    for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
    {
        int nId = it_n -> get_id();
        bool isAlive = it_n ->get_isAlive();
        if (isAlive)
        {
            double prob = randNum01b[ll++];
            if (prob < RATE_MUT_PER_DT) // the cell will acquire a random mutation
            {
                vector<int> vecMutIds = vecNode[nId].get_vecMutIds();

                int mut_id = static_cast<int>(std::floor(randNum01b[ll++] * NUM_MUT_POOL));
                if (find(vecMutIds.begin(), vecMutIds.end(), mut_id) == vecMutIds.end())    // dont duplicate in writing
                    vecMutIds.push_back(mut_id);

                vecNode[nId].set_vecMutIds(vecMutIds);
            }
        }
    }
}

// function updateVoxMap() is an iterative function
// ... no need to check if cell isAlive
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode)
{
    /* indexing i = nz*(dim.i*dim.j) + ny*(dim.k) + nx */
    iVec dim = voxMap.get_dim();
    int nMax = dim.i*dim.j*dim.k;
    double rMax = TUMOUR_RADIUS_MAX;
    double aVoxX = rMax*2./dim.i;
    double aVoxY = rMax*2./dim.j;
    double aVoxZ = rMax*2./dim.k;

    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        int bId = it->get_id();
        dVec coord = it->get_coord();

        double xx, yy, zz;
        int nx, ny, nz;
        xx = (coord.x+rMax) / aVoxX;
        yy = (coord.y+rMax) / aVoxY;
        zz = (coord.z+rMax) / aVoxZ;

        if (xx < 0 || yy < 0 || zz < 0) // warn negative values
        {
            cout << "coord: " << coord.x << ", " << coord.y << ", " << coord.z << endl;
            cout << "xx = " << xx << "; yy = " << yy << "; zz = " << zz << endl;
            cout << "Please offset coordinates so that x, y, z values are positive !" << endl;
            exit(234);
        }

        nx = (int) floor( xx );
        ny = (int) floor( yy );
        nz = (int) floor( zz );
        int n = nz*(dim.i*dim.j) + ny*dim.i + nx;
        mapBeadIds[n].push_back(bId);

        // OFFSET in space to examine if this Node is within Halo space of neighboring 26 boxes
        int nx2, ny2, nz2;
        double aOffsetX = 4.*CELL_RADIUS/aVoxX, aOffsetY = 4.*CELL_RADIUS/aVoxY, aOffsetZ = 4.*CELL_RADIUS/aVoxZ;

        // TYPE 1 : 6 faces
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x
            if (nx2 == nx+1 && n+1 < nMax)
                mapBeadIdsHalo[n+1].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x
            if (nx2 == nx-1 && n-1 >= 0)
                mapBeadIdsHalo[n-1].push_back(bId);
            // ---------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y
            if (ny2 == ny+1 && n+dim.i < nMax)
                mapBeadIdsHalo[n+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y
            if (ny2 == ny-1 && n-dim.i >= 0)
                mapBeadIdsHalo[n-dim.i].push_back(bId);
            // ---------------------------------------------------
            nz2 = (int) floor( zz+aOffsetZ );         // plus y
            if (nz2 == nz+1 && n+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y
            if (nz2 == nz-1 && n-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i*dim.j].push_back(bId);
        }
        // TYPE 2 : 12 edges
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx+1 && ny2 == ny+1 && n+1+dim.i < nMax)
                mapBeadIdsHalo[n+1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y
            if (nx2 == nx+1 && ny2 == ny-1 && n+1-dim.i >= 0)
                mapBeadIdsHalo[n+1-dim.i].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx-1 && ny2 == ny+1 && n-1+dim.i < nMax)
                mapBeadIdsHalo[n-1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y
            if (nx2 == nx-1 && ny2 == ny-1 && n-1-dim.i >= 0)
                mapBeadIdsHalo[n-1-dim.i].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && nz2 == nz+1 && n+1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus z
            if (nx2 == nx+1 && nz2 == nz-1 && n+1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i*dim.j].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && nz2 == nz+1 && n-1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus z
            if (nx2 == nx-1 && nz2 == nz-1 && n-1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny+1 && nz2 == nz+1 && n+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus y, minus z
            if (ny2 == ny+1 && nz2 == nz-1 && n+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny-1 && nz2 == nz+1 && n-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y, minus z
            if (ny2 == ny-1 && nz2 == nz-1 && n-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i-dim.i*dim.j].push_back(bId);
        }
        // TYPE 3 : 8 corners
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz+1 && n+1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, plus y, minus z
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz-1 && n+1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz+1 && n+1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus y, minus z
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz-1 && n+1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz+1 && n-1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, plus y, minus z
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz-1 && n-1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz+1 && n-1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus y, minus z
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz-1 && n-1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i-dim.i*dim.j].push_back(bId);
        }
    }
    voxMap.set_mapBeadIds(mapBeadIds);
    voxMap.set_mapBeadIdsHalo(mapBeadIdsHalo);

    // printing
    bool flagPrint = false;
    if (flagPrint)
    {
        for (unordered_map<int, vector<int>>::iterator it = mapBeadIds.begin(); it != mapBeadIds.end(); it ++)
        {
            int boxId = it->first;
            vector<int> vecBeadIds = it->second;
            cout << "boxId = " << boxId << " contains beadIds: ";
            for (vector<int>::iterator it2 = vecBeadIds.begin(); it2 != vecBeadIds.end(); it2 ++)
                cout << *it2 << "\t";

            unordered_map<int, vector<int>>::const_iterator got = mapBeadIdsHalo.find(boxId);
            if ( got != mapBeadIdsHalo.end() )
            {
                cout << "with Halo beadIds: ";
                for (vector<int>::iterator it3 = mapBeadIdsHalo[boxId].begin(); it3 != mapBeadIdsHalo[boxId].end(); it3 ++)
                    cout << *it3 << "\t";
            }

            cout << "\n" << endl;
        }
    }
}
