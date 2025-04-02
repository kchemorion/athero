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
#ifndef LDLRULE_H_
#define LDLRULE_H_

#include "biodynamo.h"
#include "Agents.h"
#include "Network.h"
#include "AgentMotion.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <string>


namespace bdm {

class ECMEnvironment;

// Structure to hold the info for LDL modification processes
struct TransformationInfo {
    enum class Type { oxLDL,  LDL };  //agLDL,
    Type type;
    Real3 position;
    double diameter;
    int color;
    
    TransformationInfo(Type t, const Real3& pos, double diam, int col)
        : type(t), position(pos), diameter(diam), color(col) {}
};

class LDLTransformationBehavior : public Behavior {
  private:
    int total_ldl_particles_;
    Simulation* simulation_;
    Random* random_;
    
    // Constant for aggregation ratio
   // static const int LDL_PER_AGGREGATE = 10;
    
    // Safety mechanism - only run on certain timesteps
  //  int timestep_counter_ = 0;
   // int run_frequency_ = 10;  // Run every 10 timesteps for safety

  public:
    LDLTransformationBehavior(int total_ldl_particles, 
                             Simulation* simulation, Random* random);
        
    Behavior* New() const override;
    Behavior* NewCopy() const override;
    
    void Run(Agent* agent) override; 
    
    // Helper methods
    void CreateAgent(const TransformationInfo& info);
 };
} //namespace bdm

#endif  // LDLRULE_H_
