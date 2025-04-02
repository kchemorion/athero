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

//This file defines the agents in the atherosclerosis simulation

#ifndef BDM_AGENTS_H_
#define BDM_AGENTS_H_


#include "biodynamo.h"
#include "Network.h"
#include "LDLrule.h"
#include "AgentMotion.h"
#include "core/agent/cell.h"



namespace bdm {


// Define the endothelial cell
class EndothelialCell : public Cell, public INetworkEnabled {    
  BDM_AGENT_HEADER(EndothelialCell, Cell, 1);
 private:
  NetworkModule network_;           
 public:
  EndothelialCell() {}
  explicit EndothelialCell(const Real3& position) : Base(position) {}
  virtual ~EndothelialCell() {}
  
 
  void InitializeNetwork(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& mact,
        const std::vector<std::vector<double>>& minh,
        int h) {
        network_ = NetworkModule(cell_type, node_names, initial_state, mact, minh, h);
    }
    
    // Implement the interface method 
    NetworkModule& GetNetwork() override { return network_; }
    
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  
};


// Define the vascular smooth muscle cell

class SmoothMuscleCell : public Cell, public INetworkEnabled {
  BDM_AGENT_HEADER(SmoothMuscleCell, Cell, 1);
  
 private:
  NetworkModule network_;           
 public:
  int cell_color_;
  bool elastin_is_degraded;
  //active, switch, lipid counter
  
  SmoothMuscleCell() {}
  explicit SmoothMuscleCell(const Real3& position) : Base(position) {}
  virtual ~SmoothMuscleCell() {}
  
  void InitializeNetwork(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& mact,
        const std::vector<std::vector<double>>& minh,
        int h) {
        network_ = NetworkModule(cell_type, node_names, initial_state, mact, minh, h);
    }
    
    // Implement the interface method 
    NetworkModule& GetNetwork() override { return network_; }
    
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
};

// Define the monocyte
class Monocyte : public Cell {
  BDM_AGENT_HEADER(Monocyte, Cell, 1);

 public:
  Monocyte() {}
  explicit Monocyte(const Real3& position) : Base(position) {}
  virtual ~Monocyte() {}
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  //birth
};


//define the macrophage
class Macrophage : public Cell {
  BDM_AGENT_HEADER(Macrophage, Cell, 1);
  
 public:
  Macrophage() {}
  explicit Macrophage(const Real3& position) : Base(position) {}
  virtual ~Macrophage() {}
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  bool is_active ();
  //age, lifespan, ingested lipid counter
  
};


//define the T cell
class Tcell : public Cell {
  BDM_AGENT_HEADER(Tcell, Cell, 1);

 public:
  Tcell() {}
  explicit Tcell(const Real3& position) : Base(position) {}
  virtual ~Tcell() {}
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  //age
  
};

//define the foam cell
class FoamCell : public Cell {
  BDM_AGENT_HEADER(FoamCell, Cell, 1);

 public:
  FoamCell() {}
  explicit FoamCell(const Real3& position) : Base(position) {}
  virtual ~FoamCell() {}
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  bool is_dead {false};
  //age, lipid counter, bursts counter
  
};


//define oxidized LDL
class oxLDL : public Cell {   
  BDM_AGENT_HEADER(oxLDL, Cell, 1);

 public:
  oxLDL() {}
  explicit oxLDL(const Real3& position) : Base(position) {}
  virtual ~oxLDL() {}
  const double oxLDL_molecular_weight {3e6}; // molecular weight of oxLDL (g/mol)
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; } 
  
  //oxidation rule (previously every 2 ticks)
  //aggregation rule (previously every 4 ticks)
  
};

//define aggregated LDL
class agLDL : public Cell {   
  BDM_AGENT_HEADER(agLDL, Cell, 1);

 public:
  agLDL() {}
  explicit agLDL(const Real3& position) : Base(position) {}
  virtual ~agLDL() {}
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }  
};

//define low density lipoprotein (LDL)
class LDL : public Cell {   
  BDM_AGENT_HEADER(LDL, Cell, 1);

 public:
  LDL() {}
  explicit LDL(const Real3& position) : Base(position) {}
  virtual ~LDL() {}
  const double LDL_molecular_weight {3e6}; // molecular weight of LDL (g/mol)
  int LDL_concentration {}; // concentration of LDL in mg/dL, which can vary based on simulation scenario
  double LDL_norm = 1/105 * LDL_concentration - 95/105;  //normalizes the value of LDL between 0 & 1 based on its concentration to use LDL_norm as an input to nodes in the regulatory network
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
};

//define high density lipoprotein (HDL)
class HDL : public Cell {   
  BDM_AGENT_HEADER(HDL, Cell, 1);

 public:
  HDL() {}
  explicit HDL(const Real3& position) : Base(position) {}
  virtual ~HDL() {}
  
  const double HDL_molecular_weight {3e5}; // molecular weight of HDL 
  int HDL_concentration {}; // concentration of HDL in mg/dL, which can vary based on simulation scenario
  double HDL_norm = 1/50 * HDL_concentration - 0.5;  //normalizes the value of HDL between 0 & 1 based on its concentration to use HDL_norm as an input to nodes in the regulatory network
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  
  
};

//define auto-antibodies (AAB)
class AAB : public Cell {   
  BDM_AGENT_HEADER(AAB, Cell, 1);

 public:
  AAB() {}
  explicit AAB(const Real3& position) : Base(position) {}
  virtual ~AAB() {}
  
  const int AAB_molecular_weight {150000}; // molecular weight of auto-antibodies (grams per mole)
  const int AAB_concentration {500}; // concentration of AAB in mU, which is the same in all simulation cases
  int cell_color_;
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  
  //lipid counter
  
};

// ECM class to track extracellular matrix processes//////////////////////////////////////
class ECMEnvironment {
private:
    static NetworkModule* ecm_network_;

public:
    // Initialize the ECM network
    static void Initialize(
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& mact,
        const std::vector<std::vector<double>>& minh,
        int h) {
        
        if (!ecm_network_) {
            ecm_network_ = new NetworkModule("ecm", node_names, initial_state, mact, minh, h);
        }
    }
    
    // Get the network module
    static NetworkModule* GetNetwork() {
        return ecm_network_;
    }
    
    // Update the network for one timestep
    static void Update() {
        if (ecm_network_) {
            ecm_network_->Step();
        }
    }
    
    // Clean up
    static void Cleanup() {
        delete ecm_network_;
        ecm_network_ = nullptr;
    }
};

} //namespace

 
#endif  // BDM_AGENTS_H_
