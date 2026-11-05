#pragma once

#include "radar/scene.hpp"

namespace radar::scenarios {

/// Stationary sensor watching an intersection: crossing, oncoming, receding
/// and turning vehicles, two walking pedestrians, poles and a barrier.
Scene intersection(const SceneConfig& cfg = {});

/// Ego platform accelerating through a gentle left curve (8 -> ~17 m/s over
/// 30 s) with roadside poles along the route, a lead vehicle, oncoming and
/// crossing traffic and a crossing pedestrian. The radar is mounted 3.6 m
/// ahead of the rotation centre unless `cfg` sets a mount.
Scene ego_drive(const SceneConfig& cfg = {}, double duration_s = 30.0);

}  // namespace radar::scenarios
