// Radar perception node: RadarCube -> point cloud, Doppler ego-velocity and
// object tracks.
//
// Subscribes:  radar/cube (RadarCube), ego/odometry (Odometry, optional)
// Publishes:   radar/detections (PointCloud2, radar frame)
//              radar/ego_velocity (TwistWithCovarianceStamped, radar frame)
//              radar/tracks (TrackArray, odom frame), radar/track_markers (MarkerArray)
//              radar/processing_time_ms (Float64)
//
// Without odometry the node tracks in the radar frame (static sensor).

#include <geometry_msgs/msg/twist_with_covariance_stamped.hpp>
#include <memory>
#include <nav_msgs/msg/odometry.hpp>
#include <radar_perception_msgs/msg/radar_cube.hpp>
#include <radar_perception_msgs/msg/track_array.hpp>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <sensor_msgs/point_cloud2_iterator.hpp>
#include <std_msgs/msg/float64.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include "conversions.hpp"
#include "radar/camera.hpp"
#include "radar/clustering.hpp"
#include "radar/ego_velocity.hpp"
#include "radar/pipeline.hpp"
#include "radar/tracker.hpp"

namespace radar_ros {

class RadarPerceptionNode : public rclcpp::Node {
 public:
  RadarPerceptionNode() : Node("radar_perception") {
    radar::PipelineConfig pc;
    pc.rcs_offset_db = declare_parameter<double>("rcs_offset_db", 101.6);
    pc.cfar.pfa = declare_parameter<double>("cfar_pfa", 1e-6);
    pipeline_ = std::make_unique<radar::RadarPipeline>(radar::RadarConfig{}, pc);
    mount_ = {declare_parameter<double>("mount_x", 3.6), declare_parameter<double>("mount_y", 0.0),
              declare_parameter<double>("mount_yaw", 0.0)};
    odom_frame_ = declare_parameter<std::string>("odom_frame", "odom");

    sub_cube_ = create_subscription<radar_perception_msgs::msg::RadarCube>(
        "radar/cube", rclcpp::QoS(2).reliable(),
        [this](radar_perception_msgs::msg::RadarCube::ConstSharedPtr msg) { on_cube(*msg); });
    sub_odom_ = create_subscription<nav_msgs::msg::Odometry>(
        "ego/odometry", 10, [this](nav_msgs::msg::Odometry::ConstSharedPtr msg) { odom_ = msg; });

    pub_points_ = create_publisher<sensor_msgs::msg::PointCloud2>("radar/detections", 10);
    pub_ego_ =
        create_publisher<geometry_msgs::msg::TwistWithCovarianceStamped>("radar/ego_velocity", 10);
    pub_tracks_ = create_publisher<radar_perception_msgs::msg::TrackArray>("radar/tracks", 10);
    pub_markers_ =
        create_publisher<visualization_msgs::msg::MarkerArray>("radar/track_markers", 10);
    pub_latency_ = create_publisher<std_msgs::msg::Float64>("radar/processing_time_ms", 10);
  }

 private:
  /// Radar pose and velocity in the odom frame from the latest odometry.
  radar::SensorState sensor_state() const {
    radar::SensorState s;
    if (!odom_) return s;  // static sensor: track in the radar frame
    const auto& p = odom_->pose.pose;
    const radar::Pose2D body{p.position.x, p.position.y, quaternion_to_yaw(p.orientation)};
    s.pose = body.compose(mount_);
    const radar::Vec2 v_body(odom_->twist.twist.linear.x, odom_->twist.twist.linear.y);
    const double w = odom_->twist.twist.angular.z;
    s.velocity = body.rotation_matrix() * v_body +
                 radar::cross_z(w, body.rotation_matrix() * mount_.translation());
    return s;
  }

  void on_cube(const radar_perception_msgs::msg::RadarCube& msg) {
    const auto t0 = std::chrono::steady_clock::now();
    if (!msg_to_cube(msg, cube_)) {
      RCLCPP_WARN(get_logger(), "Malformed radar cube, dropping");
      return;
    }
    pipeline_->process(cube_, result_);
    const auto& dets = result_.detections;
    publish_points(msg.header, dets);

    const radar::EgoVelocityEstimate ev = ego_velocity_.estimate(dets);
    if (ev.valid) {
      geometry_msgs::msg::TwistWithCovarianceStamped tw;
      tw.header = msg.header;
      tw.twist.twist.linear.x = ev.velocity.x();
      tw.twist.twist.linear.y = ev.velocity.y();
      tw.twist.covariance[0] = ev.covariance(0, 0);
      tw.twist.covariance[1] = ev.covariance(0, 1);
      tw.twist.covariance[6] = ev.covariance(1, 0);
      tw.twist.covariance[7] = ev.covariance(1, 1);
      pub_ego_->publish(tw);
    }

    dbscan_.cluster(dets, clusters_);
    tracker_.update(to_seconds(msg.header.stamp), clusters_, sensor_state());
    publish_tracks(msg.header);

    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std_msgs::msg::Float64 lat;
    lat.data = ms;
    pub_latency_->publish(lat);
    if (++frames_ % 100 == 0) {
      RCLCPP_INFO(get_logger(), "%zu detections, %zu tracks, %.2f ms (signal %.2f ms)", dets.size(),
                  tracker_.tracks().size(), ms, result_.timings.total_ms);
    }
  }

  void publish_points(const std_msgs::msg::Header& header,
                      const std::vector<radar::Detection>& dets) {
    sensor_msgs::msg::PointCloud2 pc;
    pc.header = header;
    pc.height = 1;
    sensor_msgs::PointCloud2Modifier mod(pc);
    mod.setPointCloud2Fields(
        6, "x", 1, sensor_msgs::msg::PointField::FLOAT32, "y", 1,
        sensor_msgs::msg::PointField::FLOAT32, "z", 1, sensor_msgs::msg::PointField::FLOAT32,
        "velocity", 1, sensor_msgs::msg::PointField::FLOAT32, "snr", 1,
        sensor_msgs::msg::PointField::FLOAT32, "rcs", 1, sensor_msgs::msg::PointField::FLOAT32);
    mod.resize(dets.size());
    sensor_msgs::PointCloud2Iterator<float> x(pc, "x"), y(pc, "y"), z(pc, "z"), v(pc, "velocity"),
        snr(pc, "snr"), rcs(pc, "rcs");
    for (const auto& d : dets) {
      *x = d.x_m;
      *y = d.y_m;
      *z = 0.0f;
      *v = d.velocity_mps;
      *snr = d.snr_db;
      *rcs = d.rcs_dbsm;
      ++x, ++y, ++z, ++v, ++snr, ++rcs;
    }
    pub_points_->publish(pc);
  }

  void publish_tracks(const std_msgs::msg::Header& radar_header) {
    radar_perception_msgs::msg::TrackArray arr;
    arr.header.stamp = radar_header.stamp;
    arr.header.frame_id = odom_ ? odom_frame_ : radar_header.frame_id;
    visualization_msgs::msg::MarkerArray markers;
    for (const auto& t : tracker_.tracks()) {
      radar_perception_msgs::msg::Track m;
      m.id = t.id;
      m.confirmed = t.confirmed();
      m.classification = radar::to_string(t.cls);
      m.class_confidence = static_cast<float>(t.class_confidence);
      m.position.x = t.x(0);
      m.position.y = t.x(1);
      m.velocity.x = t.x(2);
      m.velocity.y = t.x(3);
      m.extent = static_cast<float>(t.extent_m);
      arr.tracks.push_back(m);
      if (!t.confirmed()) continue;

      const radar::ObjectSize sz = radar::nominal_size(t.cls);
      visualization_msgs::msg::Marker box;
      box.header = arr.header;
      box.ns = "tracks";
      box.id = t.id;
      box.type = visualization_msgs::msg::Marker::CUBE;
      box.pose.position = m.position;
      box.pose.orientation = yaw_to_quaternion(std::atan2(t.x(3), t.x(2)));
      box.scale.x = sz.length;
      box.scale.y = sz.width;
      box.scale.z = sz.height;
      box.color.a = 0.8f;
      box.color.r = t.cls == radar::ObjectClass::kPedestrian ? 0.92f : 0.16f;
      box.color.g = t.cls == radar::ObjectClass::kPedestrian ? 0.41f : 0.47f;
      box.color.b = t.cls == radar::ObjectClass::kPedestrian ? 0.20f : 0.84f;
      box.lifetime = rclcpp::Duration::from_seconds(0.2);
      markers.markers.push_back(box);
    }
    pub_tracks_->publish(arr);
    pub_markers_->publish(markers);
  }

  std::unique_ptr<radar::RadarPipeline> pipeline_;
  radar::EgoVelocityEstimator ego_velocity_;
  radar::Dbscan dbscan_;
  radar::MultiTargetTracker tracker_;
  radar::DataCube cube_;
  radar::FrameResult result_;
  std::vector<radar::Cluster> clusters_;
  radar::Pose2D mount_;
  std::string odom_frame_;
  nav_msgs::msg::Odometry::ConstSharedPtr odom_;
  std::uint64_t frames_ = 0;

  rclcpp::Subscription<radar_perception_msgs::msg::RadarCube>::SharedPtr sub_cube_;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr sub_odom_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_points_;
  rclcpp::Publisher<geometry_msgs::msg::TwistWithCovarianceStamped>::SharedPtr pub_ego_;
  rclcpp::Publisher<radar_perception_msgs::msg::TrackArray>::SharedPtr pub_tracks_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_markers_;
  rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr pub_latency_;
};

}  // namespace radar_ros

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<radar_ros::RadarPerceptionNode>());
  rclcpp::shutdown();
  return 0;
}
