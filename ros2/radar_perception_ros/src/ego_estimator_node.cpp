// Ego-motion node: IMU + radar Doppler ego-velocity -> odometry.
//
// Subscribes:  imu (Imu), radar/ego_velocity (TwistWithCovarianceStamped)
// Publishes:   ego/odometry (Odometry), TF odom -> base_link, static TF base_link -> radar

#include <tf2_ros/static_transform_broadcaster.h>
#include <tf2_ros/transform_broadcaster.h>

#include <geometry_msgs/msg/transform_stamped.hpp>
#include <geometry_msgs/msg/twist_with_covariance_stamped.hpp>
#include <memory>
#include <nav_msgs/msg/odometry.hpp>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>

#include "conversions.hpp"
#include "radar/ego_estimator.hpp"

namespace radar_ros {

class EgoEstimatorNode : public rclcpp::Node {
 public:
  EgoEstimatorNode() : Node("ego_estimator") {
    radar::EgoEstimatorConfig cfg;
    cfg.sensor_mount = {declare_parameter<double>("mount_x", 3.6),
                        declare_parameter<double>("mount_y", 0.0),
                        declare_parameter<double>("mount_yaw", 0.0)};
    odom_frame_ = declare_parameter<std::string>("odom_frame", "odom");
    base_frame_ = declare_parameter<std::string>("base_frame", "base_link");
    radar_frame_ = declare_parameter<std::string>("radar_frame", "radar");
    mount_ = cfg.sensor_mount;
    estimator_ = std::make_unique<radar::EgoEstimator>(cfg);

    tf_ = std::make_unique<tf2_ros::TransformBroadcaster>(*this);
    static_tf_ = std::make_unique<tf2_ros::StaticTransformBroadcaster>(*this);
    geometry_msgs::msg::TransformStamped mount_tf;
    mount_tf.header.stamp = now();
    mount_tf.header.frame_id = base_frame_;
    mount_tf.child_frame_id = radar_frame_;
    mount_tf.transform.translation.x = mount_.x;
    mount_tf.transform.translation.y = mount_.y;
    mount_tf.transform.rotation = yaw_to_quaternion(mount_.yaw);
    static_tf_->sendTransform(mount_tf);

    pub_odom_ = create_publisher<nav_msgs::msg::Odometry>("ego/odometry", 10);
    sub_imu_ = create_subscription<sensor_msgs::msg::Imu>(
        "imu", rclcpp::SensorDataQoS(),
        [this](sensor_msgs::msg::Imu::ConstSharedPtr msg) { on_imu(*msg); });
    sub_vel_ = create_subscription<geometry_msgs::msg::TwistWithCovarianceStamped>(
        "radar/ego_velocity", 10,
        [this](geometry_msgs::msg::TwistWithCovarianceStamped::ConstSharedPtr msg) {
          on_velocity(*msg);
        });
  }

 private:
  void on_imu(const sensor_msgs::msg::Imu& msg) {
    if (!estimator_->initialized()) return;  // wait for the first radar velocity
    radar::ImuSample s;
    s.time_s = to_seconds(msg.header.stamp);
    s.accel_body = {msg.linear_acceleration.x, msg.linear_acceleration.y};
    s.gyro_z = msg.angular_velocity.z;
    estimator_->propagate(s);
    publish(msg.header.stamp);
  }

  void on_velocity(const geometry_msgs::msg::TwistWithCovarianceStamped& msg) {
    const radar::Vec2 v(msg.twist.twist.linear.x, msg.twist.twist.linear.y);
    if (!estimator_->initialized()) {
      // Start at the odom origin, moving with the measured velocity (lever arm
      // ignored at start-up; the filter corrects it within a few updates).
      const radar::Vec2 v_body = mount_.rotation_matrix() * v;
      estimator_->initialize({}, v_body, to_seconds(msg.header.stamp));
      RCLCPP_INFO(get_logger(), "Initialised at %.2f m/s", v_body.norm());
      return;
    }
    radar::Mat2 cov;
    cov << msg.twist.covariance[0], msg.twist.covariance[1], msg.twist.covariance[6],
        msg.twist.covariance[7];
    if (!estimator_->update_radar_velocity(v, cov)) {
      RCLCPP_DEBUG(get_logger(), "Radar velocity rejected by the innovation gate");
    }
  }

  void publish(const builtin_interfaces::msg::Time& stamp) {
    const radar::Pose2D p = estimator_->pose();
    const radar::Vec2 vb = estimator_->velocity_body();
    const auto& P = estimator_->covariance();

    nav_msgs::msg::Odometry odo;
    odo.header.stamp = stamp;
    odo.header.frame_id = odom_frame_;
    odo.child_frame_id = base_frame_;
    odo.pose.pose.position.x = p.x;
    odo.pose.pose.position.y = p.y;
    odo.pose.pose.orientation = yaw_to_quaternion(p.yaw);
    odo.pose.covariance[0] = P(0, 0);
    odo.pose.covariance[1] = P(0, 1);
    odo.pose.covariance[6] = P(1, 0);
    odo.pose.covariance[7] = P(1, 1);
    odo.pose.covariance[35] = P(4, 4);
    odo.twist.twist.linear.x = vb.x();
    odo.twist.twist.linear.y = vb.y();
    odo.twist.twist.angular.z = estimator_->yaw_rate();
    odo.twist.covariance[0] = P(2, 2);
    odo.twist.covariance[7] = P(3, 3);
    pub_odom_->publish(odo);

    geometry_msgs::msg::TransformStamped tf;
    tf.header = odo.header;
    tf.child_frame_id = base_frame_;
    tf.transform.translation.x = p.x;
    tf.transform.translation.y = p.y;
    tf.transform.rotation = odo.pose.pose.orientation;
    tf_->sendTransform(tf);
  }

  std::unique_ptr<radar::EgoEstimator> estimator_;
  radar::Pose2D mount_;
  std::string odom_frame_, base_frame_, radar_frame_;
  std::unique_ptr<tf2_ros::TransformBroadcaster> tf_;
  std::unique_ptr<tf2_ros::StaticTransformBroadcaster> static_tf_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr pub_odom_;
  rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr sub_imu_;
  rclcpp::Subscription<geometry_msgs::msg::TwistWithCovarianceStamped>::SharedPtr sub_vel_;
};

}  // namespace radar_ros

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<radar_ros::EgoEstimatorNode>());
  rclcpp::shutdown();
  return 0;
}
