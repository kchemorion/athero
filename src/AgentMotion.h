#ifndef BDM_AGENTMOTION_H_
#define BDM_AGENTMOTION_H_

#include "biodynamo.h"
#include "core/behavior/behavior.h"
#include "Agents.h"
#include "core/agent/cell.h"
#include "core/util/math.h"
#include <unordered_set>
#include <memory>

namespace bdm {

struct BrownianMotion : public Behavior {
    BDM_BEHAVIOR_HEADER(BrownianMotion, Behavior, 1);

    BrownianMotion() = default;
    virtual ~BrownianMotion() = default;

    void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // BDM_AGENTMOTION_H_
