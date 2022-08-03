//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DetermineBonds.h"
#include <GraphMol/RDKitBase.h>
//#include "EHTTools.h"
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

namespace RDKit {

// should this be object oriented?


void connectivityHuckel(RWMol &mol, int charge) {
    unsigned int numAtoms = mol.getNumAtoms();
    mol.getAtomWithIdx(0)->setFormalCharge(charge);
    //
}

void connectivityVdW(RWMol &mol, double covFactor) {
    unsigned int numAtoms = mol.getNumAtoms();
    double *distMat = MolOps::get3DDistanceMat(mol);
    
    double rcov[numAtoms];
    for (unsigned int i = 0; i < numAtoms; i++) {
        rcov[i] = PeriodicTable::getTable()->getRcovalent(mol.getAtomWithIdx(i)->getAtomicNum());
        rcov[i] *= covFactor;
    }

    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            if (distMat[i * numAtoms + j] <= (rcov[i] + rcov[j])) {
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
        }
    }
}

void determineConnectivity(RWMol &mol, bool useHuckel, int charge, double covFactor) { // accept reference
    unsigned int numAtoms = mol.getNumAtoms();
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            mol.removeBond(i, j);
            mol.getAtomWithIdx(i)->setNoImplicit(false);
            mol.getAtomWithIdx(j)->setNoImplicit(false);
        }
    }
    if (useHuckel) {
        connectivityHuckel(mol, charge);
    } else {
        connectivityVdW(mol, covFactor);
    }
}

//std::vector<int> possibleValences(Atom *atom) { // somehow return by reference
//    unsigned int atomNum = atom->getAtomicNum();
//    unsigned int numBonds = atom->getDegree(); //look at all options
//    INT_VECT valences = PeriodicTable::getTable()->getValenceList(atomNum);
//    auto smallest = valences.begin();
//    for (auto &valence : valences) {
//        if (valence > numBonds) {
//            break;
//        }
//        smallest++; // CHECK
//    }
//    std::vector<int> possible(smallest, valences.end());
//    return possible;
//}
//
//std::vector<int> getUnsaturated(RWMol &mol, std::vector<vector<int>> possible) { // somehow return by reference
//    std::vector<int> degree(mol.getNumAtoms(), 0);
//    for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
//        if ((possible[i].size() != 0) && (possible[i][possible[i].size() - 1] > mol.getAtomWithIdx(i)->getDegree())) {
//            degree[i] = possible[i][possible[i].size() - 1] - mol.getAtomWithIdx(i)->getDegree();
//        }
//    }
//    return {0};
//}
//
//std::vector<std::vector<int>> valenceCombinations(RWMol &mol, std::vector<std::vector<int>> possible) { // return by reference
//    unsigned int numAtoms = mol.getNumAtoms();
//
//    possible.resize(numAtoms);
//    unsigned int numCombos = 1;
//    for (unsigned int i = 0; i < numAtoms; i++) {
//        possible[i] = possibleValences(mol.getAtomWithIdx(i));
//        numCombos *= possible[i].size();
//    }
//
//    unsigned int orig = numCombos;
//
//    std::vector<std::vector<int>> combos(numCombos, std::vector<int>(numAtoms));
//    for (unsigned int i = 0; i < numAtoms; i++) {
//         unsigned int seg = numCombos / possible[i].size();
//         for (unsigned int p = 0; p < orig / numCombos; p++) {
//           for (unsigned int j = 0; j < possible[i].size(); j++) {
//             for (unsigned int k = 0; k < seg; k++) {
//                 combos[p*possible[i].size()*seg + j*seg + k][i] = possible[i][j];
//             }
//           }
//         }
//         numCombos /= possible[i].size();
//     }
//
//    return combos;
//}
//
//void determineBondOrder(RWMol &mol, bool ignoreChiral, bool allowChargedFragments, int charge) {
//    for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
//
//    }
//
//    std::vector<std::vector<int>> orders = valenceCombinations(mol);
//
//    for (auto &order : orders) {
//        // find all unsaturated atoms and how many bonds they can make
//        // if there are no unsaturated atoms, check everything's valid and return
//        // otherwise, find the best possibilities for another layer of pi bonds
//        // repeat
//    }
//}


} // namespace RDKit


