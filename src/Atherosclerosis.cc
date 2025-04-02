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
#include "Atherosclerosis.h"
#include "Network.h"
#include "Agents.h"
#include "LDLrule.h"
#include "AgentMotion.h"


namespace bdm {

void SetupNetworks(Simulation* simulation) {
    auto* ctxt = simulation->GetExecutionContext();
    auto* random = simulation->GetRandom();
    size_t spacing = 20; // Spacing between cells
    
    // Define nodes across all cell types
    std::vector<std::string> node_names = {"SuperO", "LDL_norm", "HDL_norm", "oxLDL", "IFNg", "IL12", "IL6", "NO", "MCP1", "KLF2", "Nfkb", "TGFb", "IP10", "MIF", "VCAM1", "MMP9", "LRP1", "agLDL", "ERK5", "PDGFb", "TNFa", "IL1b"};
    //int num_nodes = node_names.size();
    int h = 1;  // Sensitivity parameter
    
    // Initial state of the nodes, 0.0 unless specified
    std::unordered_map<std::string, double> initial_values = {  
        {"Nfkb", 1.0},    //set ERK5 to 1 to explore the mechanobiology of laminar flow, set Nfkb to 1 to explore the mechanobiology of oscilatory shear stress
        {"HDL_norm", 1.0}, //set node values for HDL & LDL based on their concentration
        {"LDL_norm", 1.0},
    };
    
    // Create initial state with these values
    std::vector<double> initial_state = CreateInitialState(node_names, initial_values);
    
    // Define endothelial cell regulatory connections
    std::vector<std::tuple<std::string, std::string, double>> endothelial_activations = {
        {"LDL_norm", "SuperO", 1.0},
        {"oxLDL", "SuperO", 5.0},
        {"IFNg", "SuperO", 0.01},
        {"IL12", "SuperO", 0.01},
        {"IL6", "SuperO", 0.01},
        {"HDL_norm", "NO", 0.01},
        {"NO", "NO", 0.01},
        {"KLF2", "NO", 0.01},
        {"TGFb", "MCP1", 0.01},
        {"oxLDL", "MCP1", 1.0},
        {"IFNg", "MCP1", 0.01},
        {"IL12", "MCP1", 0.01},
        {"Nfkb", "MCP1", 0.01},
        {"oxLDL", "IP10", 1.0},
        {"Nfkb", "IP10", 0.01},
        {"IFNg", "IP10", 0.01},
        {"IL12", "IP10", 0.01},
        {"ERK5", "KLF2", 0.01},
        {"Nfkb", "MIF", 0.01},
        {"IL6", "VCAM1", 0.01},
        {"Nfkb", "IL6", 0.01},
        {"oxLDL", "IL6", 1.0},
        {"IL1b", "IL6", 0.01},
        {"TNFa", "IL6", 0.01},
        {"oxLDL", "IL1b", 1.0},
        {"LDL_norm", "IL1b", 0.01},
        {"TNFa", "IL1b", 0.01},
        {"IL1b", "TNFa", 0.01},
        {"oxLDL", "TNFa", 1.0},
        {"IFNg", "TNFa", 0.01},
    };
    
    std::vector<std::tuple<std::string, std::string, double>> endothelial_inhibitions = {
        {"HDL_norm", "SuperO", 0.01},
        {"NO", "SuperO", 0.01},
        {"LDL_norm", "NO", 1.0},
        {"oxLDL", "NO", 5.0},
        {"Nfkb", "NO", 0.01},
        {"IFNg", "NO", 0.01},
        {"IL12", "NO", 0.01},
        {"IL6", "NO", 0.01},
        {"SuperO", "NO", 0.01},
        {"NO", "MCP1", 0.01},
        {"KLF2", "MCP1", 0.01},
        {"NO", "IP10", 0.01},
        {"KLF2", "IP10", 0.01},
        {"oxLDL", "KLF2", 0.01},
        {"KLF2", "MIF", 0.01},
        {"KLF2", "VCAM1", 0.01},
        {"KLF2", "IL6", 0.01},
        {"HDL", "IL6", 0.01},
        {"NO", "IL6", 0.01},
        {"TGFb", "IL6", 0.01},
        {"TGFb", "IL1b", 0.01},
        {"NO", "IL1b", 0.01},
        {"TGFb", "TNFa", 0.01},
        {"HDL", "TNFa", 0.01},
        {"NO", "TNFa", 0.01},
    };
    
    auto endothelial_Mact = CreateMatrixFromConnections(node_names, endothelial_activations);
    auto endothelial_Minh = CreateMatrixFromConnections(node_names, endothelial_inhibitions);
    
    // Define smooth muscle cell regulatory connections
    std::vector<std::tuple<std::string, std::string, double>> smc_activations = {
        {"oxLDL", "IL12", 1.0},
        {"SuperO", "IL12", 0.01},
        {"agLDL", "LRP1", 0.01},
    };
    
    std::vector<std::tuple<std::string, std::string, double>> smc_inhibitions = {
        {"LRP1", "TGFb", 0.01}
    };
    
    auto smc_Mact = CreateMatrixFromConnections(node_names, smc_activations);
    auto smc_Minh = CreateMatrixFromConnections(node_names, smc_inhibitions);
    
    // Define extracellular matrix regulatory connections ///////////////////////////////////
    std::vector<std::tuple<std::string, std::string, double>> ecm_activations = {
        {"SuperO", "oxLDL", 0.01},
        {"RNS", "oxLDL", 0.01},
        {"IL6", "oxLDL", 0.01},
        {"MCP1", "oxLDL", 0.01},
        {"agLDL", "oxLDL", 0.01},
        {"TNFa", "oxLDL", 0.01},
        {"IL1b", "oxLDL", 0.01},
        {"Nfkb", "agLDL", 0.01},
        {"MMP9", "agLDL", 0.01},
        {"IL6", "agLDL", 0.01},
        {"MCP1", "agLDL", 0.01},
        {"oxLDL", "agLDL", 0.01},
        {"TNFa", "agLDL", 0.01},
        {"IL1b", "agLDL", 0.01},
    };
    
    std::vector<std::tuple<std::string, std::string, double>> ecm_inhibitions = {
        {"NO", "oxLDL", 0.01},
        {"HDL", "oxLDL", 0.01},
        {"KLF2", "oxLDL", 0.01}, ////////////////////////
        {"NO", "agLDL", 0.01},
        {"HDL", "agLDL", 0.01},
    };
    
    auto ecm_Mact = CreateMatrixFromConnections(node_names, ecm_activations); ////////////////
    auto ecm_Minh = CreateMatrixFromConnections(node_names, ecm_inhibitions);
    
    // Initialize the ECM environment////////////////////////
    ECMEnvironment::Initialize(node_names, initial_state, ecm_Mact, ecm_Minh, h);
    
    // Set any constant nodes in the ECM network
    if (auto* ecm_network = ECMEnvironment::GetNetwork()) {
        ecm_network->SetConstantNode("Nfkb", 1.0);
        ecm_network->SetConstantNode("LDL_norm", 1.0);
        ecm_network->SetConstantNode("HDL_norm", 1.0);
    }
    real_t x_coordEC, y_coordEC, z_coordEC; // endothelial cell coordinates
    
    // Create endothelial cells in the 3D grid
    for (double x = 0.5; x < 13.5; ++x) {
        for (double z = 0.5; z < 13.5; ++z) {
            x_coordEC = x * spacing; // x-coordinate
            z_coordEC = z * spacing; // z-coordinate
            y_coordEC = 12.5 * spacing; // y-coordinate
            // Create a new cell
            EndothelialCell* cell = new EndothelialCell({x_coordEC, y_coordEC, z_coordEC});
            cell->SetDiameter(20);
            cell->SetCellColor(0);
            
            cell->InitializeNetwork("endothelial", node_names, initial_state,           
                                   endothelial_Mact, endothelial_Minh, h);
            cell->GetNetwork().SetConstantNode("Nfkb", 1.0);  
            cell->GetNetwork().SetConstantNode("LDL_norm", 1.0); 
            cell->GetNetwork().SetConstantNode("HDL_norm", 1.0);                         
            cell->AddBehavior(new NetworkUpdateBehavior());           
            // Store reference to the first cell for analysis 
            if (endothelial_cell_ptr == nullptr) {
                // Need to cast since we're storing as INetworkEnabled*
                endothelial_cell_ptr = static_cast<INetworkEnabled*>(cell);
            }
            // Add cells to simulation
            ctxt->AddAgent(cell);
        }
    }
    
    real_t x_coordSMC, y_coordSMC, z_coordSMC; // smooth muscle cell coordinates
    
    // Create smooth muscle cells in the 3D grid
    for (double x = 0.5; x < 13.5; ++x) {
        for (double y = 0.5; y < 6.5; ++y) {
            for (double z = 0.5; z < 13.5; ++z) {
                x_coordSMC = x * spacing; // x-coordinate
                y_coordSMC = y * spacing; // y-coordinate
                z_coordSMC = z * spacing; // z-coordinate
                SmoothMuscleCell* cell = new SmoothMuscleCell({x_coordSMC, y_coordSMC, z_coordSMC});
                cell->SetDiameter(20);
                cell->SetCellColor(1);
                cell->InitializeNetwork("smc", node_names, initial_state,           
                                   endothelial_Mact, endothelial_Minh, h);
            //cell->GetNetwork().SetConstantNode("ERK5", 1.0);  
            cell->GetNetwork().SetConstantNode("LDL_norm", 1.0); 
            cell->GetNetwork().SetConstantNode("HDL_norm", 1.0);                          
            cell->AddBehavior(new NetworkUpdateBehavior());            
            // Store reference to the first cell for analysis 
            if (smc_cell_ptr == nullptr) {
                // Need to cast since we're storing as INetworkEnabled*
                smc_cell_ptr = static_cast<INetworkEnabled*>(cell);
            }
                ctxt->AddAgent(cell);
            }
        }
    }
    
 real_t x_coord_monocyte, y_coord_monocyte, z_coord_monocyte; // random coordinates for the single monocyte initially in the intimal space
 
 // Initialize a single monocyte randomly positioned within the boundaries of the intimal space
 

   x_coord_monocyte = random->Uniform(0.0, 260.0);  //boundaries of the intima
   y_coord_monocyte = random->Uniform(120, 240.0);
   z_coord_monocyte = random->Uniform(0.0, 260.0);
   Monocyte* cell = new Monocyte({x_coord_monocyte, y_coord_monocyte, z_coord_monocyte});
    // set cell parameters
   cell->SetDiameter(12);
   cell->SetCellColor(2);
   cell->AddBehavior(new BrownianMotion());
   ctxt->AddAgent(cell);  // put the created cell in our cell structure
  
 //create LDL particles in the intima
  
  constexpr double Avogadro_number = 6.022e23; // the Avogadro number is used to translate concentrations into the number of particles in the intima
  constexpr double intimal_volume = 8.112e6;       //the volume of the intima in micrometers cubed
  constexpr double mgdL_to_gmicrom = 1e-17;        // converts mg/dL to gram per micrometer cubed
  int LDL_concentration = 200;
  const double LDL_molecular_weight = 3e6; //grams per mole, approximated value
  double LDL_mgdl_to_particles = (LDL_concentration * mgdL_to_gmicrom * Avogadro_number * intimal_volume) / LDL_molecular_weight;
  int LDL_particles = static_cast<int>(std::round(LDL_mgdl_to_particles / 1e6));
  
  std::cout << "number of LDL particles is " << LDL_particles << std::endl;

  real_t x_coord_LDL, y_coord_LDL, z_coord_LDL; // random coordinates for the LDL particles in the intimal space
  
  //Initialize LDL particles randomly positioned within the boundaries of the intimal space
  
  for (int i=0; i<LDL_particles; ++i) {
    x_coord_LDL = random->Uniform(80.0, 200.0);  //boundaries of the intima
    y_coord_LDL = random->Uniform(120.0, 200.0);
    z_coord_LDL = random->Uniform(80.0, 200.0);
    LDL* cell = new LDL({x_coord_LDL, y_coord_LDL, z_coord_LDL});
    // set particle parameters
    cell->SetDiameter(0.5);   //0.025
    cell->SetCellColor(3);
    cell->AddBehavior(new BrownianMotion());
    ctxt->AddAgent(cell);  // put the created molecules in our cell structure
  }
  
    //Create special controller agent for ECM processes and LDL modification (oxidation & aggregation)///////////////////////////////////////////////
    Cell* ecm_controller = new Cell({0, 0, 0});
    ecm_controller->SetDiameter(0.1);
    
    // Add ECM update behavior
    ecm_controller->AddBehavior(new ECMUpdateBehavior());
    
    // Add LDL transformation behavior
    ecm_controller->AddBehavior(new LDLTransformationBehavior(
        LDL_particles,   // Total regular LDL particles to maintain
        simulation,      // Simulation
        random           // Random number generator
    ));
    
    // Add controller to simulation/////////////////////////////////
    ctxt->AddAgent(ecm_controller);
    
  //create HDL particles in the intima
  
  int HDL_concentration = 75;
  const double HDL_molecular_weight = 3e5;
  
  double HDL_mgdl_to_particles = (HDL_concentration * mgdL_to_gmicrom * Avogadro_number * intimal_volume) / HDL_molecular_weight;
  int HDL_particles = static_cast<int>(std::round(HDL_mgdl_to_particles / 1e7));
  
  std::cout << "number of HDL particles is " << HDL_particles << std::endl;

  real_t x_coord_HDL, y_coord_HDL, z_coord_HDL; // random coordinates for the HDL particles in the intimal space
  
  //Initialize HDL particles randomly positioned within the boundaries of the intimal space
  
  for (int i=0; i<HDL_particles; ++i) {
    x_coord_HDL = random->Uniform(0.0, 260.0);  //boundaries of the intima
    y_coord_HDL = random->Uniform(120.0, 240.0);
    z_coord_HDL = random->Uniform(0.0, 260.0);
    HDL* cell = new HDL({x_coord_HDL, y_coord_HDL, z_coord_HDL});
  // set particle parameters
    cell->SetDiameter(0.5);  //0.01
    cell->SetCellColor(4);
    cell->AddBehavior(new BrownianMotion());
    ctxt->AddAgent(cell);  // put the created lipid molecules in our cell structure
  }

 //Initialize auto-antibodies in the intima
  
  const double AAB_molecular_weight {1.5e5}; // molecular weight of auto-antibodies (grams per mole)
  const int AAB_concentration {500};
  double mcg_mcl_AAB = AAB_concentration * 3.84 / 100; //converts the AAB concentration in mU to microgram/microliter
  double AAB_molecules = 1e5 * (mcg_mcl_AAB * mgdL_to_gmicrom * Avogadro_number * intimal_volume) / AAB_molecular_weight;
  int AAB_particles = static_cast<int>(std::round(AAB_molecules / 1e12));
  
  std::cout << "number of AAB particles is " << AAB_particles << std::endl;

  real_t x_coord_AAB, y_coord_AAB, z_coord_AAB; // random coordinates for the AAB particles in the intimal space
  
  //Create AAB particles randomly positioned within the boundaries of the intimal space
  
  for (int i=0; i<AAB_particles; ++i) {
    x_coord_AAB = random->Uniform(20.0, 240.0);  //boundaries of the intima
    y_coord_AAB = random->Uniform(140.0, 220.0);
    z_coord_AAB = random->Uniform(20.0, 240.0);
    AAB* cell = new AAB({x_coord_AAB, y_coord_AAB, z_coord_AAB});
  // set particle parameters
    cell->SetDiameter(0.5);  //0.015
    cell->SetCellColor(5);
    cell->AddBehavior(new BrownianMotion());
    ctxt->AddAgent(cell);  // put the created lipid molecules in our cell structure
  }
  
  


}
//define the boundary conditions of the 3D environment (a cube with periodic boundaries)

inline int Simulate(int argc, const char** argv) {
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kTorus;
    param->min_bound = 0; 
    param->max_bound = 260; 
    param->simulation_time_step = 0.01;  
  };

  Simulation simulation(argc, argv, set_param);
  // Setup networks and create cells
  SetupNetworks(&simulation);
  
 // try {
    // Run the simulations
    simulation.GetScheduler()->Simulate(100);
    
    // Ensure memory cleanup
    //PositionTracker::Cleanup();
  
    std::cout << "Simulation completed successfully!" << std::endl;
    
     if (endothelial_cell_ptr) {
      std::cout << "Endothelial NO: " << endothelial_cell_ptr->GetNetwork().GetNodeState("NO") << std::endl; ////////////////////
      
      std::cout << "Endothelial Nfkb: " << endothelial_cell_ptr->GetNetwork().GetNodeState("Nfkb") << std::endl;
      std::cout << "Endothelial MIF: " << endothelial_cell_ptr->GetNetwork().GetNodeState("MIF") << std::endl;
      std::cout << "Endothelial VCAM1: " << endothelial_cell_ptr->GetNetwork().GetNodeState("VCAM1") << std::endl;
      std::cout << "Endothelial IL6: " << endothelial_cell_ptr->GetNetwork().GetNodeState("IL6") << std::endl;
      std::cout << "Endothelial MCP1: " << endothelial_cell_ptr->GetNetwork().GetNodeState("MCP1") << std::endl;
      std::cout << "Endothelial superoxide: " << endothelial_cell_ptr->GetNetwork().GetNodeState("SuperO") << std::endl;
      std::cout << "Endothelial KLF2: " << endothelial_cell_ptr->GetNetwork().GetNodeState("KLF2") << std::endl;
      std::cout << "Endothelial IP10: " << endothelial_cell_ptr->GetNetwork().GetNodeState("IP10") << std::endl;
      std::cout << "Endothelial ERK5: " << endothelial_cell_ptr->GetNetwork().GetNodeState("ERK5") << std::endl;
      std::cout << "EC HDL: " << endothelial_cell_ptr->GetNetwork().GetNodeState("HDL_norm") << std::endl;
      std::cout << "EC LDL: " << endothelial_cell_ptr->GetNetwork().GetNodeState("LDL_norm") << std::endl;
      std::cout << "EC oxLDL: " << endothelial_cell_ptr->GetNetwork().GetNodeState("oxLDL") << std::endl;
    } else {
      std::cout << "No endothelial cell was tracked" << std::endl;
    }
    if (smc_cell_ptr) {
      std::cout << "SMC NO: " << smc_cell_ptr->GetNetwork().GetNodeState("NO") << std::endl;
      std::cout << "SMC Nfkb: " << smc_cell_ptr->GetNetwork().GetNodeState("Nfkb") << std::endl;
      std::cout << "SMC MCP1: " << smc_cell_ptr->GetNetwork().GetNodeState("MIF") << std::endl;
      std::cout << "SMC VCAM1: " << smc_cell_ptr->GetNetwork().GetNodeState("VCAM1") << std::endl;
      std::cout << "SMC IL6: " << smc_cell_ptr->GetNetwork().GetNodeState("IL6") << std::endl;
      std::cout << "SMC MCP1: " << smc_cell_ptr->GetNetwork().GetNodeState("MCP1") << std::endl;
      std::cout << "SMC superoxide: " << smc_cell_ptr->GetNetwork().GetNodeState("SuperO") << std::endl;
      std::cout << "SMC KLF2: " << smc_cell_ptr->GetNetwork().GetNodeState("KLF2") << std::endl;
      std::cout << "SMC IP10: " << smc_cell_ptr->GetNetwork().GetNodeState("IP10") << std::endl;
      std::cout << "SMC ERK5: " << smc_cell_ptr->GetNetwork().GetNodeState("ERK5") << std::endl;
      std::cout << "SMC HDL: " << smc_cell_ptr->GetNetwork().GetNodeState("HDL_norm") << std::endl;
      std::cout << "SMC LDL: " << smc_cell_ptr->GetNetwork().GetNodeState("LDL_norm") << std::endl;
      std::cout << "SMC oxLDL: " << smc_cell_ptr->GetNetwork().GetNodeState("oxLDL") << std::endl;
    } else {
      std::cout << "No SMC cell was tracked" << std::endl;
    }
    
    return 0;
    
 /* } catch (const std::exception& e) {
    // Attempt cleanup even in case of exception
    try {
      PositionTracker::Cleanup();
    } catch (...) {
      std::cerr << "Error during emergency cleanup" << std::endl;
    }
    
    std::cerr << "Simulation failed with error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    // Handle unknown exceptions
    try {
      PositionTracker::Cleanup();
    } catch (...) {
      std::cerr << "Error during emergency cleanup" << std::endl;
    }
    
    std::cerr << "Simulation failed with unknown error" << std::endl;
    return 1;
  } */
//}
} 
}// namespace

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }




