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

#include "biodynamo.h"
#include "Agents.h"
#include "Network.h"
#include "LDLrule.h"
#include "AgentMotion.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <string>

//LDL modification 

namespace bdm {

LDLTransformationBehavior::LDLTransformationBehavior(int total_ldl_particles, 
                             Simulation* simulation, Random* random)
        : total_ldl_particles_(total_ldl_particles),
          simulation_(simulation),
          random_(random) {}
    
    Behavior* LDLTransformationBehavior::New() const { 
        return new LDLTransformationBehavior(total_ldl_particles_, simulation_, random_); 
    }
    
    Behavior* LDLTransformationBehavior::NewCopy() const { 
        return new LDLTransformationBehavior(*this); 
    }
    
    
    void LDLTransformationBehavior::CreateAgent(const TransformationInfo& info) {
       switch (info.type) {
           case TransformationInfo::Type::oxLDL: {
              auto* ox_ldl = new oxLDL(info.position);
              ox_ldl->SetDiameter(info.diameter);
              ox_ldl->SetCellColor(info.color);
              ox_ldl->AddBehavior(new BrownianMotion());
              simulation_->GetResourceManager()->AddAgent(ox_ldl);
              break;
        }
        /*
           case TransformationInfo::Type::agLDL: {
              auto* ag_ldl = new agLDL(info.position);
              ag_ldl->SetDiameter(info.diameter);
              ag_ldl->SetCellColor(info.color);
              // No Brownian motion for aggregated LDL since they stick to the proteoglycans in the intimal space
              simulation_->GetResourceManager()->AddAgent(ag_ldl);
              break;
             
        }
        */
           case TransformationInfo::Type::LDL: {
              auto* new_ldl = new LDL(info.position);
              new_ldl->SetDiameter(info.diameter);
              new_ldl->SetCellColor(info.color);
              new_ldl->AddBehavior(new BrownianMotion());
              simulation_->GetResourceManager()->AddAgent(new_ldl);
              break;
        }
    }
}
    void LDLTransformationBehavior::Run(Agent* agent) {
        // Get the ECM network module
        NetworkModule* ecm_network = ECMEnvironment::GetNetwork();
        if (!ecm_network) return;
        
        // Get oxidation and aggregation rates from ECM network
        double ox_ldl_rate = ecm_network->GetNodeState("oxLDL");
       // double ag_ldl_rate = ecm_network->GetNodeState("agLDL");
        
        // For normal LDL, use half the aggregation rate, since oxidized LDL aggregates at a higher rate
      //  double normal_ldl_agg_rate = ag_ldl_rate / 2.0;
        
        // Get all cells in the simulation
        auto* rm = simulation_->GetResourceManager();
        
        // Separate vectors for normal and oxidized LDL
        std::vector<LDL*> normal_ldl_particles;
        std::vector<oxLDL*> oxidized_ldl_particles;
        
       
        // Collect all regular LDL particles
        rm->ForEachAgent([&](Agent* agent) {
            if (!agent) return; // Ensure agent is not null
            // First check if this agent is a Cell
            if (Cell* cell = dynamic_cast<Cell*>(agent)) {
        // Then check if it's specifically an LDL cell
                if (auto* ldl = dynamic_cast<LDL*>(cell)) {
                    if (typeid(*ldl) == typeid(LDL)) {
                       normal_ldl_particles.push_back(ldl);
                    }
                }
                // Check for oxidized LDL
                else if (auto* ox_ldl = dynamic_cast<oxLDL*>(cell)) {
                oxidized_ldl_particles.push_back(ox_ldl);
                }
            }
        });
        
       // Count existing particles of each type for logging
       int initial_normal_ldl_count = normal_ldl_particles.size();
       int initial_oxidized_ldl_count = oxidized_ldl_particles.size();
    //   int initial_agg_ldl_count = 0;
    /*
       rm->ForEachAgent([&](Agent* agent) {
            if (Cell* cell = dynamic_cast<Cell*>(agent)) {
                if (dynamic_cast<agLDL*>(cell)) {
                initial_agg_ldl_count++;
                }
            }
            
        });
        */
      // Calculate how many to transform
      int to_oxidize = std::min(static_cast<int>(normal_ldl_particles.size() * ox_ldl_rate), 
                             static_cast<int>(normal_ldl_particles.size()));
    
      // Collections for safe manipulation
      std::vector<AgentUid> to_remove;
      std::vector<TransformationInfo> to_create;
      
      // OXIDATION PHASE: Process normal LDL oxidation
      int oxidized = 0;
      for (size_t i = 0; i < normal_ldl_particles.size() && oxidized < to_oxidize; i++) {
        LDL* ldl = normal_ldl_particles[i];
        
        // Store properties
        Real3 position = ldl->GetPosition();
        double diameter = ldl->GetDiameter();
        
        // Mark for removal
        to_remove.push_back(ldl->GetUid());
        
        // Mark for creation
        to_create.push_back(TransformationInfo(
            TransformationInfo::Type::oxLDL,
            position,
            diameter,
            4 // Color for OxidizedLDL
        ));
        
        oxidized++;
    }
    
    // Remove the oxidized particles from consideration for aggregation
    normal_ldl_particles.erase(normal_ldl_particles.begin(), 
                             normal_ldl_particles.begin() + oxidized);
    
    // NORMAL LDL AGGREGATION PHASE
    /*
    int normal_groups_possible = normal_ldl_particles.size() / LDL_PER_AGGREGATE;
    int normal_groups_to_aggregate = std::min(
        static_cast<int>(normal_groups_possible * normal_ldl_agg_rate),
        normal_groups_possible);
    int normal_ldl_to_aggregate = normal_groups_to_aggregate * LDL_PER_AGGREGATE;
    
    // OXIDIZED LDL AGGREGATION PHASE
    int oxidized_groups_possible = oxidized_ldl_particles.size() / LDL_PER_AGGREGATE;
    int oxidized_groups_to_aggregate = std::min(
        static_cast<int>(oxidized_groups_possible * ag_ldl_rate),
        oxidized_groups_possible);
    int oxidized_ldl_to_aggregate = oxidized_groups_to_aggregate * LDL_PER_AGGREGATE;
    
    // Track aggregation stats
    int normal_ldl_aggregated = 0;
    int oxidized_ldl_aggregated = 0;
    int agg_groups_created = 0;
    
    // Aggregate normal LDL
    for (size_t i = 0; i < normal_ldl_particles.size() && i < static_cast<size_t>(normal_ldl_to_aggregate); 
         i += LDL_PER_AGGREGATE) {
        
        // Make sure we have enough particles left
        if (i + LDL_PER_AGGREGATE > normal_ldl_particles.size()) break;
        
        // Use the position of the last LDL in the group
        Real3 position = normal_ldl_particles[i + LDL_PER_AGGREGATE - 1]->GetPosition();
        double diameter = normal_ldl_particles[i]->GetDiameter() * 3.0;
        
        // Mark all LDL particles in this group for removal
        for (size_t j = 0; j < LDL_PER_AGGREGATE; j++) {
            to_remove.push_back(normal_ldl_particles[i + j]->GetUid());
            normal_ldl_aggregated++;
        }
        
        // Mark for creation
        to_create.push_back(TransformationInfo(
            TransformationInfo::Type::agLDL,
            position,
            diameter,
            5 // Color for AggregatedLDL
        ));
        
        agg_groups_created++;
    }
    
    // Aggregate oxidized LDL
    for (size_t i = 0; i < oxidized_ldl_particles.size() && i < static_cast<size_t>(oxidized_ldl_to_aggregate); 
         i += LDL_PER_AGGREGATE) {
        
        // Make sure we have enough particles left
        if (i + LDL_PER_AGGREGATE > oxidized_ldl_particles.size()) break;
        
        // Use the position of the last LDL in the group
        Real3 position = oxidized_ldl_particles[i + LDL_PER_AGGREGATE - 1]->GetPosition();
        double diameter = oxidized_ldl_particles[i]->GetDiameter() * 3.0;
        
        // Mark all oxidized LDL particles in this group for removal
        for (size_t j = 0; j < LDL_PER_AGGREGATE; j++) {
            to_remove.push_back(oxidized_ldl_particles[i + j]->GetUid());
            oxidized_ldl_aggregated++;
        }
        
        // Mark for creation
        to_create.push_back(TransformationInfo(
            TransformationInfo::Type::agLDL,
            position,
            diameter,
            5 // Color for AggregatedLDL
        ));
        
        agg_groups_created++;
    }
    */
    
    // CREATE REPLACEMENT NORMAL LDL PHASE
    // Only replace normal LDL that was oxidized or aggregated, not oxidized LDL
    int normal_ldl_consumed = oxidized; //+ normal_ldl_aggregated;
    int to_create_count = normal_ldl_consumed;
    
    // Prepare replacement positions
    for (int i = 0; i < to_create_count; i++) {
        Real3 position = {
            random_->Uniform(0.0, 260.0),
            random_->Uniform(120.0, 240.0),
            random_->Uniform(0.0, 260.0)
        };
        
        // Mark for creation
        to_create.push_back(TransformationInfo(
            TransformationInfo::Type::LDL,
            position,
            0.5, // Standard diameter
            3    // Color for NormalLDL
        ));
    }
    
    // EXECUTION PHASE - first remove all marked agents
    for (const auto& uid : to_remove) {
        rm->RemoveAgent(uid);
    }
    
    // Then create all new agents
    for (const auto& info : to_create) {
        CreateAgent(info);
    }
        
       // Log transformation statistics
       std::cout << "LDL Transformation: " << std::endl
          << "  - Normal LDL: " << initial_normal_ldl_count 
          << " (Oxidized: " << oxidized 
         // << ", Aggregated: " << normal_ldl_aggregated 
          << ", Replaced: " << to_create_count << ")" << std::endl
          << "  - Oxidized LDL: " << initial_oxidized_ldl_count 
        //  << " (Aggregated: " << oxidized_ldl_aggregated 
          << ", Not replaced)" << std::endl
         // << "  - Aggregated groups created: " << agg_groups_created << std::endl
          << "  - Rates: oxLDL=" << ox_ldl_rate << std::endl;
         // << ", agLDL=" << ag_ldl_rate 
         // << " (normal LDL agg rate=" << normal_ldl_agg_rate << ")" << std::endl;
    }
    
   } //namespace

