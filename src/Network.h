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
#ifndef NETWORK_H_
#define NETWORK_H_

#include "biodynamo.h"
#include "Agents.h"
#include "LDLrule.h"
#include "AgentMotion.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <string>

namespace bdm {

class Simulation;

// Define constant
constexpr double DELTA = 1.0; // Replace with your actual value

/**
 * Helper function to calculate derivatives based on cell type-specific
 * activation and inhibition matrices
 */
std::vector<double> CalculateDerivatives(
    const std::vector<double>& X_i,
    const std::vector<std::vector<double>>& Mact,
    const std::vector<std::vector<double>>& Minh,
    int h,
    int NumOfNodes);

/**
 * Network module to be used within Cell classes
 */
class NetworkModule {
  private:
    // Cell type identifier
    std::string cell_type_;
    
    //constant nodes
    std::unordered_map<std::string, double> constant_nodes_;
    
    // Node information
    std::vector<std::string> node_names_;
    std::unordered_map<std::string, int> node_name_to_index_;
    
    // Current state of all nodes
    std::vector<double> state_;
    
    // Cell type-specific activation and inhibition matrices
    std::unordered_map<std::string, std::vector<std::vector<double>>> cell_type_Mact_;
    std::unordered_map<std::string, std::vector<std::vector<double>>> cell_type_Minh_;
    
    // Parameter h (sensitivity parameter)
    int h_;
    
    // Number of nodes in the network
    int num_nodes_;
    
  public:
    NetworkModule();
    
    // Constructor with cell type and parameters
    NetworkModule(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        int h);
  
    
    // Constructor that takes pre-built cell type matrices
    NetworkModule(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& cell_Mact,
        const std::vector<std::vector<double>>& cell_Minh,
        int h);
        
    // method to keep the value of a node constant throughout the simulation
    void SetConstantNode(const std::string& node_name, double value);
        
    
    // Register a new cell type's regulatory matrices
    void RegisterCellType(
        const std::string& cell_type,
        const std::vector<std::vector<double>>& type_Mact,
        const std::vector<std::vector<double>>& type_Minh);
    
    // Get current cell type
    const std::string& GetCellType() const;
    
    // Run a single integration step using Runge-Kutta 4th order method
    void Step();
    
    // Run full integration for a specified number of steps
    void RunFullIntegration(int num_steps);
    
    // Get current state vector
    const std::vector<double>& GetStateVector() const;
    
    // Get state of a specific node by name
    double GetNodeState(const std::string& node_name) const;
    
    // Set state of a specific node by name
    void SetNodeState(const std::string& node_name, double value);
    
    // Get node names
    const std::vector<std::string>& GetNodeNames() const;
    
    // Get parameter h
    int GetH() const;
};

//  NetworkCell class with public network member
class NetworkCell : public Cell {
  public:
    // Make NetworkModule public so it can be accessed directly
    NetworkModule network;
    
    NetworkCell();
    
    // Constructor with network parameters
    NetworkCell(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& cell_Mact,
        const std::vector<std::vector<double>>& cell_Minh,
        int h);
        
    //function to set constant nodes  
    void SetConstantNode(const std::string& node_name, double value) {
        network.SetConstantNode(node_name, value);
    }
    
    // Required BioDynaMo methods for Cell
    Cell* New() const override;
    Cell* NewCopy() const override;
    
    
    void UpdateNetwork();
    
    
};

class INetworkEnabled {
  public:
    virtual ~INetworkEnabled() = default;
    virtual NetworkModule& GetNetwork() = 0;
};

// BioDynaMo behavior for explicit network updates 
class NetworkUpdateBehavior : public Behavior {           
  public:
    NetworkUpdateBehavior();
    
    // Required by BioDynaMo
    Behavior* New() const override;
    Behavior* NewCopy() const override;
    
    void Run(Agent* agent) override;
};

// Behavior to update the ECM network//////////////////////
class ECMUpdateBehavior : public Behavior {
  public:
    ECMUpdateBehavior();
    
    Behavior* New() const override; 
    Behavior* NewCopy() const override; 
    
    void Run(Agent* agent) override; 
};

extern INetworkEnabled* endothelial_cell_ptr;
extern INetworkEnabled* smc_cell_ptr;

// Helper function to create initial state from a name->value map
std::vector<double> CreateInitialState(
    const std::vector<std::string>& node_names,
    const std::unordered_map<std::string, double>& initial_values);

// Create connection matrix from a list of interactions
std::vector<std::vector<double>> CreateMatrixFromConnections(
    const std::vector<std::string>& node_names,
    const std::vector<std::tuple<std::string, std::string, double>>& connections);


} // namespace bdm

#endif  // NETWORK_H_
