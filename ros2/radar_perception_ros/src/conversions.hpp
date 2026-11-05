#pragma once

#include <cmath>
#include <cstring>
#include <geometry_msgs/msg/quaternion.hpp>
#include <radar_perception_msgs/msg/radar_cube.hpp>
#include <rclcpp/time.hpp>
#include <string>

#include "radar/geometry.hpp"
#include "radar/types.hpp"

namespace radar_ros {

inline geometry_msgs::msg::Quaternion yaw_to_quaternion(double yaw) {
  geometry_msgs::msg::Quaternion q;
  q.z = std::sin(0.5 * yaw);
  q.w = std::cos(0.5 * yaw);
  return q;
}

inline double quaternion_to_yaw(const geometry_msgs::msg::Quaternion& q) {
  return std::atan2(2.0 * (q.w * q.z + q.x * q.y), 1.0 - 2.0 * (q.y * q.y + q.z * q.z));
}

inline double to_seconds(const builtin_interfaces::msg::Time& t) {
  return static_cast<double>(t.sec) + 1e-9 * static_cast<double>(t.nanosec);
}

inline void cube_to_msg(const radar::DataCube& cube, radar_perception_msgs::msg::RadarCube& msg) {
  msg.num_rx = static_cast<uint32_t>(cube.num_rx());
  msg.num_chirps = static_cast<uint32_t>(cube.num_chirps());
  msg.num_samples = static_cast<uint32_t>(cube.num_samples());
  msg.iq.resize(2 * cube.size());
  // std::complex<float> is layout-compatible with float[2].
  std::memcpy(msg.iq.data(), cube.data(), msg.iq.size() * sizeof(float));
}

inline bool msg_to_cube(const radar_perception_msgs::msg::RadarCube& msg, radar::DataCube& cube) {
  const std::size_t n = static_cast<std::size_t>(msg.num_rx) * msg.num_chirps * msg.num_samples;
  if (msg.iq.size() != 2 * n) return false;
  if (cube.num_rx() != msg.num_rx || cube.num_chirps() != msg.num_chirps ||
      cube.num_samples() != msg.num_samples) {
    cube.resize(msg.num_rx, msg.num_chirps, msg.num_samples);
  }
  std::memcpy(cube.data(), msg.iq.data(), msg.iq.size() * sizeof(float));
  return true;
}

}  // namespace radar_ros
