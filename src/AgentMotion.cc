#include "AgentMotion.h"

namespace bdm {



void BrownianMotion::Run(Agent* agent) {
    auto* random = Simulation::GetActive()->GetRandom();
    // Cast to Cell if needed
    Cell* cell = dynamic_cast<Cell*>(agent);
    if (!cell) return;
   
    // Get current position
    Real3 current_position = cell->GetPosition();
    Real3 displacement = {random->Uniform(-0.5, 0.5),
                          random->Uniform(-0.5, 0.5),
                          random->Uniform(-0.5, 0.5)
  };
    // Calculate new position
    Real3 new_position = current_position + displacement;
    // Update position - BioDynaMo will handle the grid updates at the appropriate time
    cell->SetPosition(new_position);
}
}  // namespace bdm
