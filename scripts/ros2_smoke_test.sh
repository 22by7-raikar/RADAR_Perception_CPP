#!/usr/bin/env bash
# Launches the ROS 2 demo and checks that every output topic produces messages.
# Run from the colcon workspace root after `colcon build`.
set -eo pipefail

source /opt/ros/"${ROS_DISTRO:-jazzy}"/setup.bash
source install/setup.bash
set -u

# Own process group so the whole launch tree can be stopped. (Background jobs
# in a non-interactive shell ignore SIGINT, so SIGINT alone would not stop it.)
setsid ros2 launch radar_perception_ros demo.launch.py duration:=20.0 > launch.log 2>&1 &
LAUNCH_PID=$!
cleanup() {
  kill -TERM -- "-$LAUNCH_PID" 2>/dev/null || true
  for _ in $(seq 1 10); do
    kill -0 "$LAUNCH_PID" 2>/dev/null || return 0
    sleep 1
  done
  kill -KILL -- "-$LAUNCH_PID" 2>/dev/null || true
}
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

sample=$(timeout 30 ros2 topic echo --once --field data /radar/processing_time_ms 2>/dev/null | head -n 1 || true)
echo "processing time per radar frame: ${sample:-n/a} ms"

if [ "$status" -ne 0 ]; then
  echo "---- launch log ----"
  cat launch.log
fi
exit "$status"
