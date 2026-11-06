"""Simulated drive through the full radar stack.

    ros2 launch radar_perception_ros demo.launch.py [scenario:=intersection] [duration:=30.0]
        [array:=tdm_mimo_3tx4rx]

The simulator publishes /clock; the other nodes run on simulation time.
"""

from pathlib import Path

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description():
    params = str(Path(get_package_share_directory("radar_perception_ros")) / "config" / "params.yaml")
    scenario = LaunchConfiguration("scenario")
    duration = LaunchConfiguration("duration")
    array = LaunchConfiguration("array")  # antenna layout; source and perception must agree
    return LaunchDescription([
        DeclareLaunchArgument("scenario", default_value="ego_drive"),
        DeclareLaunchArgument("duration", default_value="30.0"),
        DeclareLaunchArgument("array", default_value="ula_8rx"),
        Node(package="radar_perception_ros", executable="sim_source_node", name="radar_sim_source",
             parameters=[params, {"scenario": scenario, "duration": duration,
                                  "array": array}], output="screen"),
        Node(package="radar_perception_ros", executable="radar_perception_node",
             name="radar_perception",
             parameters=[params, {"use_sim_time": True, "array": array}], output="screen"),
        Node(package="radar_perception_ros", executable="ego_estimator_node", name="ego_estimator",
             parameters=[params, {"use_sim_time": True}], output="screen"),
    ])
