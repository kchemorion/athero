// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef ATHEROSCLEROSIS_H_
#define ATHEROSCLEROSIS_H_

#include "biodynamo.h"
#include "Agents.h"
#include "AgentMotion.h"
#include "Network.h"
#include "LDLrule.h"
#include "core/agent/cell.h"
#include "core/environment/environment.h"
#include "core/simulation.h"
#include "core/param/param.h"
#include "core/util/random.h"
#include <cmath>

namespace bdm {

  void SetupNetworks(Simulation* simulation);

}  // namespace bdm

#endif  // ATHEROSCLEROSIS_H_






