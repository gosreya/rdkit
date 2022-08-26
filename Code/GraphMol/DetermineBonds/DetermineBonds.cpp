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
#include <YAeHMOP/EHTTools.h>
#include <iostream>

namespace RDKit {

void connectivityHuckel(RWMol &mol, int charge) {
    auto numAtoms = mol.getNumAtoms();
    mol.getAtomWithIdx(0)->setFormalCharge(charge);
    EHTTools::EHTResults res;
    bool success = runMol(mol, res);
    double *mat = res.reducedOverlapPopulationMatrix.get();
    int matInd = 0;
    int range = 1;
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = 0; j < range; j++) {
            if (mat[matInd] >= 0.15 && i != j) {
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
            matInd++;
        }
        range++;
    }
}

void connectivityVdW(RWMol &mol, double covFactor) {
    auto numAtoms = mol.getNumAtoms();
    double *distMat = MolOps::get3DDistanceMat(mol);
    
    double rcov[numAtoms];
    for (unsigned int i = 0; i < numAtoms; i++) {
        rcov[i] = covFactor * PeriodicTable::getTable()->getRcovalent(mol.getAtomWithIdx(i)->getAtomicNum());
    }
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            if (distMat[i * numAtoms + j] <= (rcov[i] + rcov[j])) {
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
        }
    }
}

void determineConnectivity(RWMol &mol, bool useHuckel, int charge, double covFactor) {
    auto numAtoms = mol.getNumAtoms();
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            mol.removeBond(i, j);
            mol.getAtomWithIdx(i)->setNoImplicit(true);
            mol.getAtomWithIdx(j)->setNoImplicit(true);
        }
    }
    if (useHuckel) {
        connectivityHuckel(mol, charge);
    } else {
        connectivityVdW(mol, covFactor);
    }
}

} // namespace RDKit


