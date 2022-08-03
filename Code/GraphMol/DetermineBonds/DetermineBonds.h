//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#pragma once
#include <GraphMol/RDKitBase.h>

// check include gaurds 

namespace RDKit {

void determineConnectivity(RWMol &mol, bool useHuckel=false, int charge=0, double covFactor=1.3);

void determineBondOrder(RWMol &mol, bool ignoreChiral=false, bool allowChargedFragments=false, int charge=0);


} // namespace RDKit

