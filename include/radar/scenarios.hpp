#pragma once

#include "radar/scene.hpp"

namespace radar::scenarios {

/// Stationary sensor watching an intersection: crossing, oncoming, receding
/// and turning vehicles, two walking pedestrians, poles and a barrier.
Scene intersection(const SceneConfig& cfg = {});

/// Same environment observed from a moving, turning ego platform.
Scene ego_drive(const SceneConfig& cfg = {});

}  // namespace radar::scenarios
