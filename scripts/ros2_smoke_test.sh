#!/usr/bin/env bash
# Launches the ROS 2 demo and checks that every output topic produces messages.
# Run from the colcon workspace root after `colcon build`.
set -eo pipefail

source /opt/ros/"${ROS_DISTRO:-jazzy}"/setup.bash
source install/setup.bash
set -u

ros2 launch radar_perception_ros demo.launch.py duration:=20.0 > launch.log 2>&1 &
LAUNCH_PID=$!
cleanup() { kill -INT "$LAUNCH_PID" 2>/dev/null || true; wait "$LAUNCH_PID" 2>/dev/null || true; }
trap cleanup EXIT

status=0
for topic in /radar/detections /radar/ego_velocity /radar/tracks /ego/odometry /radar/processing_time_ms; do
  if timeout 90 ros2 topic echo --once "$topic" > /dev/null 2>&1; then
    echo "ok      $topic"
  else
    echo "MISSING $topic"
    status=1
  fi
done

if timeout 30 ros2 topic echo --once /radar/processing_time_ms --field data 2>/dev/null | head -1 | grep -q .; then
  echo "processing time sample: $(timeout 30 ros2 topic echo --once /radar/processing_time_ms --field data | head -1) ms"
fi

if [ "$status" -ne 0 ]; then
  echo "---- launch log ----"
  cat launch.log
fi
exit "$status"
