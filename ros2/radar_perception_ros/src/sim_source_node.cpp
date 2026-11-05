// Simulated sensor source: publishes radar cubes (20 Hz), IMU samples (100 Hz),
// ground truth and /clock. Simulation time advances as fast as the machine
// allows, so downstream nodes should run with use_sim_time:=true.

#include <chrono>
#include <memory>
#include <nav_msgs/msg/odometry.hpp>
#include <radar_perception_msgs/msg/radar_cube.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rosgraph_msgs/msg/clock.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include "conversions.hpp"
#include "radar/imu.hpp"
#include "radar/scenarios.hpp"
#include "radar/simulator.hpp"

namespace radar_ros {

class SimSourceNode : public rclcpp::Node {
 public:
  SimSourceNode() : Node("radar_sim_source") {
    const auto scenario = declare_parameter<std::string>("scenario", "ego_drive");
    duration_s_ = declare_parameter<double>("duration", 30.0);
    const auto seed = static_cast<std::uint32_t>(declare_parameter<int>("seed", 3));
    base_frame_ = declare_parameter<std::string>("base_frame", "base_link");
    radar_frame_ = declare_parameter<std::string>("radar_frame", "radar");
    world_frame_ = declare_parameter<std::string>("world_frame", "world");

    scene_ = std::make_unique<radar::Scene>(scenario == "intersection"
                                                ? radar::scenarios::intersection()
                                                : radar::scenarios::ego_drive({}, duration_s_));
    radar::SimulatorConfig sc;
    sc.seed = seed;
    sim_ = std::make_unique<radar::Simulator>(radar::RadarConfig{}, sc);
    radar::ImuConfig ic;
    ic.seed = seed + 10;
    imu_ = std::make_unique<radar::ImuSimulator>(ic);

    // Cubes are ~2 MB: keep them reliable with a short queue.
    pub_cube_ = create_publisher<radar_perception_msgs::msg::RadarCube>("radar/cube",
                                                                        rclcpp::QoS(2).reliable());
    pub_imu_ = create_publisher<sensor_msgs::msg::Imu>("imu", rclcpp::SensorDataQoS());
    pub_clock_ = create_publisher<rosgraph_msgs::msg::Clock>("/clock", rclcpp::ClockQoS());
    pub_truth_ = create_publisher<nav_msgs::msg::Odometry>("ground_truth/odometry", 10);
    pub_markers_ =
        create_publisher<visualization_msgs::msg::MarkerArray>("ground_truth/objects", 10);
    timer_ = create_wall_timer(std::chrono::milliseconds(1), [this] { step(); });
    RCLCPP_INFO(get_logger(), "Simulating '%s' for %.0f s", scenario.c_str(), duration_s_);
  }

 private:
  void step() {
    if (time_s_ > duration_s_) {
      RCLCPP_INFO(get_logger(), "Simulation finished");
      timer_->cancel();
      return;
    }
    const rclcpp::Time stamp(static_cast<int64_t>(time_s_ * 1e9), RCL_ROS_TIME);
    rosgraph_msgs::msg::Clock clock;
    clock.clock = stamp;
    pub_clock_->publish(clock);

    const radar::ImuSample s = imu_->sample(scene_->ego(), time_s_);
    sensor_msgs::msg::Imu imu;
    imu.header.stamp = stamp;
    imu.header.frame_id = base_frame_;
    imu.linear_acceleration.x = s.accel_body.x();
    imu.linear_acceleration.y = s.accel_body.y();
    imu.angular_velocity.z = s.gyro_z;
    pub_imu_->publish(imu);

    if (tick_ % 5 == 0) publish_radar(stamp);
    scene_->step(kImuDt);
    time_s_ += kImuDt;
    ++tick_;
  }

  void publish_radar(const rclcpp::Time& stamp) {
    const radar::SceneFrame f = scene_->snapshot();
    sim_->simulate_frame(f.targets, cube_);
    auto msg = std::make_unique<radar_perception_msgs::msg::RadarCube>();
    msg->header.stamp = stamp;
    msg->header.frame_id = radar_frame_;
    cube_to_msg(cube_, *msg);
    pub_cube_->publish(std::move(msg));

    nav_msgs::msg::Odometry odo;
    odo.header.stamp = stamp;
    odo.header.frame_id = world_frame_;
    odo.child_frame_id = base_frame_;
    odo.pose.pose.position.x = f.ego.pose.x;
    odo.pose.pose.position.y = f.ego.pose.y;
    odo.pose.pose.orientation = yaw_to_quaternion(f.ego.pose.yaw);
    odo.twist.twist.linear.x = f.ego.speed;
    odo.twist.twist.angular.z = f.ego.yaw_rate;
    pub_truth_->publish(odo);

    visualization_msgs::msg::MarkerArray markers;
    for (const auto& g : f.objects) {
      visualization_msgs::msg::Marker m;
      m.header = odo.header;
      m.ns = "ground_truth";
      m.id = g.id;
      m.type = visualization_msgs::msg::Marker::CUBE;
      m.pose.position.x = g.position_world.x();
      m.pose.position.y = g.position_world.y();
      m.pose.orientation = yaw_to_quaternion(g.heading);
      m.scale.x = std::max(0.3, g.length);
      m.scale.y = std::max(0.3, g.width);
      m.scale.z = 1.5;
      m.color.r = m.color.g = m.color.b = 0.6f;
      m.color.a = 0.4f;
      m.lifetime = rclcpp::Duration::from_seconds(0.2);
      markers.markers.push_back(m);
    }
    pub_markers_->publish(markers);
  }

  static constexpr double kImuDt = 0.01;
  std::unique_ptr<radar::Scene> scene_;
  std::unique_ptr<radar::Simulator> sim_;
  std::unique_ptr<radar::ImuSimulator> imu_;
  radar::DataCube cube_;
  double duration_s_ = 30.0;
  double time_s_ = 0.0;
  std::uint64_t tick_ = 0;
  std::string base_frame_, radar_frame_, world_frame_;

  rclcpp::Publisher<radar_perception_msgs::msg::RadarCube>::SharedPtr pub_cube_;
  rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr pub_imu_;
  rclcpp::Publisher<rosgraph_msgs::msg::Clock>::SharedPtr pub_clock_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr pub_truth_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_markers_;
  rclcpp::TimerBase::SharedPtr timer_;
};

}  // namespace radar_ros

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<radar_ros::SimSourceNode>());
  rclcpp::shutdown();
  return 0;
}
