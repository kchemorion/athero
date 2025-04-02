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

namespace bdm {

// Define constant
//constexpr double DELTA = 1.0; // Replace with your actual value

/**
 * Helper function to calculate derivatives based on cell type-specific
 * activation and inhibition matrices
 */
std::vector<double> CalculateDerivatives(
    const std::vector<double>& X_i,
    const std::vector<std::vector<double>>& Mact,
    const std::vector<std::vector<double>>& Minh,
    int h,
    int NumOfNodes) {
    
    // Initialize vectors for calculations
    std::vector<double> w(NumOfNodes, 0.0); // Total input of a node
    std::vector<double> f(NumOfNodes, 0.0); // The derivatives
    std::vector<double> Ract(NumOfNodes);   // Activating links for a node
    std::vector<double> Rinh(NumOfNodes);   // Inhibiting links for a node
    std::vector<int> gamma(NumOfNodes, 1); // decay term initialized to 1 for all nodes
    
    // Calculate weights and derivatives for each node
    for (int i = 0; i < NumOfNodes; i++) {
        // Get activating and inhibiting connections for node i
        Ract = Mact[i];
        Rinh = Minh[i];
        
        // Calculate sum of activation weights
        double sumAlfa = 0.0;
        for (size_t j = 0; j < Ract.size(); j++) {
            sumAlfa += Ract[j];
        }
        double sumAlphaX = 0.0;
        
        // Calculate sum of inhibition weights
        double sumBeta = 0.0;
        for (size_t j = 0; j < Rinh.size(); j++) {
            sumBeta += Rinh[j];
        }
        double sumBetaX = 0.0;
        
        // Calculate weight based on connection types
        if (sumBeta == 0 && sumAlfa > 0) {
            // Node i has no inhibitors and at least 1 activator
            for (size_t k = 0; k < Ract.size(); k++) {
                sumAlphaX += Ract[k] * X_i[k];
            }
            w[i] = ((1.0 + sumAlfa) / sumAlfa) * (sumAlphaX / (1.0 + sumAlphaX));
            
        } else if (sumBeta > 0 && sumAlfa == 0) {
            // Node i has no activators and at least 1 inhibitor
            for (size_t k = 0; k < Rinh.size(); k++) {
                sumBetaX += Rinh[k] * X_i[k];
            }
            w[i] = 1.0 - ((1.0 + sumBeta) / sumBeta) * (sumBetaX / (1.0 + sumBetaX));
            
        } else if (sumBeta > 0 && sumAlfa > 0) {
            // Node i has at least one activator and one inhibitor
            for (size_t k = 0; k < Ract.size(); k++) {
                sumAlphaX += Ract[k] * X_i[k];
            }
            for (size_t k = 0; k < Rinh.size(); k++) {
                sumBetaX += Rinh[k] * X_i[k];
            }
            w[i] = (((1.0 + sumAlfa) / sumAlfa) * (sumAlphaX / (1.0 + sumAlphaX))) * 
                   (1.0 - ((1.0 + sumBeta) / sumBeta) * (sumBetaX / (1.0 + sumBetaX)));
                   
        } else {
            // Node i has no activators nor inhibitors
            w[i] = 0.0;
        }
        
        // Compute the derivative for node i
        f[i] = DELTA * ((-std::exp(0.5 * h) + std::exp(-h * (w[i] - 0.5))) / 
               ((1.0 - std::exp(0.5 * h)) * (1.0 + std::exp(-h * (w[i] - 0.5))))) - 
           gamma[i] * X_i[i];
    }
    
    return f; // Return vector with the computed derivatives
}


    NetworkModule::NetworkModule() {}
    
    // Constructor with cell type and parameters
    NetworkModule::NetworkModule(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        int h)
        : cell_type_(cell_type),
          node_names_(node_names),
          state_(initial_state),
          h_(h),
          num_nodes_(node_names.size())
    {
        // Build node name to index mapping
        for (size_t i = 0; i < node_names.size(); i++) {
            node_name_to_index_[node_names[i]] = i;
        }
    }
    
    // Constructor that takes pre-built cell type matrices
    NetworkModule::NetworkModule(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& cell_Mact,
        const std::vector<std::vector<double>>& cell_Minh,
        int h)
        : cell_type_(cell_type),
          node_names_(node_names),
          state_(initial_state),
          h_(h),
          num_nodes_(node_names.size())
    {
        // Build node name to index mapping
        for (size_t i = 0; i < node_names.size(); i++) {
            node_name_to_index_[node_names[i]] = i;
        }
        
        // Set matrices for this cell type
        cell_type_Mact_[cell_type] = cell_Mact;
        cell_type_Minh_[cell_type] = cell_Minh;
    }
    
    // method to keep the value of a node constant throughout the simulation
    void NetworkModule::SetConstantNode(const std::string& node_name, double value) {
        auto it = node_name_to_index_.find(node_name);
        if (it != node_name_to_index_.end()) {
            state_[it->second] = value;  // Set initial value
            constant_nodes_[node_name] = value;  // Remember it's constant
        } else {
            throw std::runtime_error("Cannot set constant node - name not found: " + node_name);
        }
    }
    
    // Register a new cell type's regulatory matrices
    void NetworkModule::RegisterCellType(
        const std::string& cell_type,
        const std::vector<std::vector<double>>& type_Mact,
        const std::vector<std::vector<double>>& type_Minh) {
        
        cell_type_Mact_[cell_type] = type_Mact;
        cell_type_Minh_[cell_type] = type_Minh;
    }
    
    // Get current cell type
    const std::string& NetworkModule::GetCellType() const {
        return cell_type_;
    }
    
    // Run a single integration step using Runge-Kutta 4th order method
    void NetworkModule::Step() {
        // Ensure we have matrices for this cell type
        if (cell_type_Mact_.find(cell_type_) == cell_type_Mact_.end() ||
            cell_type_Minh_.find(cell_type_) == cell_type_Minh_.end()) {
            throw std::runtime_error("No network matrices defined for cell type: " + cell_type_);
        }
        
        // Get the right matrices for this cell type
        const auto& Mact = cell_type_Mact_[cell_type_];
        const auto& Minh = cell_type_Minh_[cell_type_];
        
        // k1 calculation
        auto k1 = CalculateDerivatives(state_, Mact, Minh, h_, num_nodes_);
        
        // k2 calculation
        std::vector<double> k2factor(num_nodes_);
        for (int i = 0; i < num_nodes_; i++) {
            k2factor[i] = state_[i] + k1[i] * 0.5;
        }
        auto k2 = CalculateDerivatives(k2factor, Mact, Minh, h_, num_nodes_);
        
        // k3 calculation
        std::vector<double> k3factor(num_nodes_);
        for (int i = 0; i < num_nodes_; i++) {
            k3factor[i] = state_[i] + k2[i] * 0.5;
        }
        auto k3 = CalculateDerivatives(k3factor, Mact, Minh, h_, num_nodes_);
        
        // k4 calculation
        std::vector<double> k4factor(num_nodes_);
        for (int i = 0; i < num_nodes_; i++) {
            k4factor[i] = state_[i] + k3[i];
        }
        auto k4 = CalculateDerivatives(k4factor, Mact, Minh, h_, num_nodes_);
        
        // Update state
        for (int i = 0; i < num_nodes_; i++) {
            state_[i] = state_[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
        
        //keep the value of constant nodes unchanged
        for (const auto& pair : constant_nodes_) {
            state_[node_name_to_index_[pair.first]] = pair.second;
        }
    }
    
    // Run full integration for a specified number of steps
    void NetworkModule::RunFullIntegration(int num_steps) {
        for (int time = 0; time < num_steps; time++) {
            Step();
        }
    }
    
    // Get current state vector
    const std::vector<double>& NetworkModule::GetStateVector() const {
        return state_;
    }
    
    // Get state of a specific node by name
    double NetworkModule::GetNodeState(const std::string& node_name) const {
        auto it = node_name_to_index_.find(node_name);
        if (it != node_name_to_index_.end()) {
            return state_[it->second];
        }
        throw std::runtime_error("Node name not found: " + node_name);
    }
    
    // Set state of a specific node by name
    void NetworkModule::SetNodeState(const std::string& node_name, double value) {
        auto it = node_name_to_index_.find(node_name);
        if (it != node_name_to_index_.end()) {
            state_[it->second] = value;
        } else {
            throw std::runtime_error("Node name not found: " + node_name);
        }
    }
    
    // Get node names
    const std::vector<std::string>& NetworkModule::GetNodeNames() const {
        return node_names_;
    }
    
    // Get parameter h
    int NetworkModule::GetH() const {
        return h_;
    }
    
    NetworkCell::NetworkCell() : Cell() {}
    
    // Constructor with network parameters
    NetworkCell::NetworkCell(
        const std::string& cell_type,
        const std::vector<std::string>& node_names,
        const std::vector<double>& initial_state,
        const std::vector<std::vector<double>>& cell_Mact,
        const std::vector<std::vector<double>>& cell_Minh,
        int h) 
        : Cell(),
          network(cell_type, node_names, initial_state, cell_Mact, cell_Minh, h) {}
    
    // Required BioDynaMo methods for Cell
    Cell* NetworkCell::New() const { return new NetworkCell(); }
    Cell* NetworkCell::NewCopy() const  { return new NetworkCell(*this); }
    
    
    void NetworkCell::UpdateNetwork() {
        // Update network dynamics
        network.Step();  
    }
    
   
    NetworkUpdateBehavior::NetworkUpdateBehavior() {}
    
    // Required by BioDynaMo
    Behavior* NetworkUpdateBehavior::New() const  { return new NetworkUpdateBehavior(); }
    Behavior* NetworkUpdateBehavior::NewCopy() const  { return new NetworkUpdateBehavior(*this); }
    
    void NetworkUpdateBehavior::Run(Agent* agent)  {
        // Cast to NetworkCell and update network if successful
        auto* networkEnabled = dynamic_cast<INetworkEnabled*>(agent);               
        if (networkEnabled) {
            networkEnabled->GetNetwork().Step();                                                                          
                                                                                                       
        }
    }
 
 // Behavior to update the ECM network//////////////////////////////////

    ECMUpdateBehavior::ECMUpdateBehavior() {}
    
    Behavior* ECMUpdateBehavior::New() const  { return new ECMUpdateBehavior(); }
    Behavior* ECMUpdateBehavior::NewCopy() const  { return new ECMUpdateBehavior(*this); }
    
    void ECMUpdateBehavior::Run(Agent* agent) {
        // Update the ECM network
        ECMEnvironment::Update();
    }


INetworkEnabled* endothelial_cell_ptr = nullptr;
INetworkEnabled* smc_cell_ptr = nullptr;    

// Initialize static member////////////////////////////////////
NetworkModule* ECMEnvironment::ecm_network_ = nullptr;// Global ECM network module to track extracellular processes        

// Helper function to create initial state from a name->value map
std::vector<double> CreateInitialState(
    const std::vector<std::string>& node_names,
    const std::unordered_map<std::string, double>& initial_values) {
    
    // Create zero-initialized vector
    std::vector<double> state(node_names.size(), 0.0);
    
    // Build node name to index mapping
    std::unordered_map<std::string, int> name_to_idx;
    for (size_t i = 0; i < node_names.size(); i++) {
        name_to_idx[node_names[i]] = i;
    }
    
    // Set initial values by name
    for (const auto& pair : initial_values) {
        auto it = name_to_idx.find(pair.first);
        if (it != name_to_idx.end()) {
            state[it->second] = pair.second;
        }
    }
    
    return state;
}

// Create connection matrix from a list of interactions
std::vector<std::vector<double>> CreateMatrixFromConnections(
    const std::vector<std::string>& node_names,
    const std::vector<std::tuple<std::string, std::string, double>>& connections) {
    
    int num_nodes = node_names.size();
    std::vector<std::vector<double>> matrix(num_nodes, std::vector<double>(num_nodes, 0.0));
    
    // Create node name to index mapping
    std::unordered_map<std::string, int> name_to_idx;
    for (size_t i = 0; i < node_names.size(); i++) {
        name_to_idx[node_names[i]] = i;
    }
    
    // Set connections 
    for (const auto& connection : connections) {
        std::string source = std::get<0>(connection);
        std::string target = std::get<1>(connection);
        double weight = std::get<2>(connection);
        
        auto source_it = name_to_idx.find(source);
        auto target_it = name_to_idx.find(target);
        
        if (source_it != name_to_idx.end() && target_it != name_to_idx.end()) {
            int source_idx = source_it->second;
            int target_idx = target_it->second;
            matrix[target_idx][source_idx] = weight;
        }
    }
    
    return matrix;
}
} // namespace bdm

