#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "radar/calibration.hpp"
#include "radar/camera.hpp"
#include "radar/fusion.hpp"

namespace radar {

namespace {

constexpr double kDeg = kPi / 180.0;

Track make_track(int id, Vec2 p, Vec2 v, ObjectClass cls, double conf = 0.8) {
  Track t;
  t.id = id;
  t.status = TrackStatus::kConfirmed;
  t.x << p.x(), p.y(), v.x(), v.y();
  t.cls = cls;
  t.class_confidence = conf;
  return t;
}

/// Camera detection exactly where the object is (radar frame at camera time).
CameraDetection detection_for(const CameraModel& cam, Vec2 centre, double heading, ObjectClass cls,
                              double radar_height) {
  const ObjectSize s = nominal_size(cls);
  CameraDetection d;
  EXPECT_TRUE(cam.project_box(centre, heading, s.length, s.width, s.height, -radar_height, d.box));
  d.cls = cls;
  d.score = 0.9;
  return d;
}

}  // namespace

TEST(Camera, AxesAndProjection) {
  CameraModel cam;
  cam.T_cam_from_radar = make_camera_extrinsic(Vec3::Zero(), 0.0, 0.0, 0.0);
  const Vec3 p = cam.T_cam_from_radar * Vec3(10.0, 2.0, 1.0);
  EXPECT_TRUE(p.isApprox(Vec3(-2.0, -1.0, 10.0)));  // right, down, forward
  Eigen::Vector2d uv;
  ASSERT_TRUE(cam.project(Vec3(20.0, 0.0, 0.0), uv));
  EXPECT_NEAR(uv.x(), cam.K.cx, 1e-9);
  EXPECT_NEAR(uv.y(), cam.K.cy, 1e-9);
  EXPECT_FALSE(cam.project(Vec3(-5.0, 0.0, 0.0), uv));  // behind
  ASSERT_TRUE(cam.project(Vec3(20.0, 2.0, 0.0), uv));   // left of boresight
  EXPECT_LT(uv.x(), cam.K.cx);
}

TEST(Camera, PitchMisalignmentShiftsImageVertically) {
  CameraModel level, pitched;
  level.T_cam_from_radar = make_camera_extrinsic(Vec3::Zero(), 0.0, 0.0, 0.0);
  pitched.T_cam_from_radar = make_camera_extrinsic(Vec3::Zero(), 0.0, 2.0 * kDeg, 0.0);
  Eigen::Vector2d a, b;
  level.project(Vec3(30.0, 0.0, 0.0), a);
  pitched.project(Vec3(30.0, 0.0, 0.0), b);
  // Pitching down (positive about +y with z up) moves the scene up in the image.
  EXPECT_NEAR(std::abs(b.y() - a.y()), level.K.fy * std::tan(2.0 * kDeg), 0.5);
}

TEST(Camera, BoxIou) {
  const BoundingBox a{0, 0, 10, 10}, b{5, 0, 15, 10}, c{20, 20, 30, 30};
  EXPECT_NEAR(a.iou(b), 50.0 / 150.0, 1e-12);
  EXPECT_DOUBLE_EQ(a.iou(c), 0.0);
  EXPECT_DOUBLE_EQ(a.iou(a), 1.0);
}

TEST(Calibration, RecoversExtrinsicFromRoughGuess) {
  CameraModel truth;
  truth.T_cam_from_radar =
      make_camera_extrinsic(Vec3(-1.8, 0.1, 0.8), 0.4 * kDeg, 2.0 * kDeg, -0.7 * kDeg);
  const Pose3 rough =
      make_camera_extrinsic(Vec3(-1.6, 0.0, 0.95), -1.5 * kDeg, -1.0 * kDeg, 2.5 * kDeg);
  std::vector<Correspondence> data;
  std::mt19937 rng(1);
  std::normal_distribution<double> n(0.0, 1.0);
  for (double r : {8.0, 15.0, 25.0}) {
    for (double az : {-20.0, 0.0, 20.0}) {
      for (double z : {-0.3, 0.6}) {
        Correspondence c;
        c.p_radar = Vec3(r * std::cos(az * kDeg), r * std::sin(az * kDeg), z);
        ASSERT_TRUE(truth.project(c.p_radar, c.uv));
        c.uv += 0.5 * Eigen::Vector2d(n(rng), n(rng));
        data.push_back(c);
      }
    }
  }
  const auto res = calibrate_extrinsic(truth.K, data, rough);
  EXPECT_TRUE(res.converged);
  EXPECT_GT(res.initial_rms_px, 20.0);
  EXPECT_LT(res.final_rms_px, 1.0);
  const PoseError e = pose_difference(truth.T_cam_from_radar, res.T_cam_from_radar);
  EXPECT_LT(e.rotation_rad, 0.1 * kDeg);
  EXPECT_LT(e.translation_m, 0.05);
}

TEST(Fusion, AssociatesProjectedTrackWithCameraBox) {
  CameraModel cam;
  cam.T_cam_from_radar = make_camera_extrinsic(Vec3(-1.8, 0.0, 0.8), 0.0, 0.0, 0.0);
  RadarCameraFusion fusion(cam);
  SensorState s;
  // Radar sees the rear surface of a vehicle whose centre is 2.3 m further.
  const auto trk = make_track(1, {30.0, 1.0}, {0.0, 0.0}, ObjectClass::kVehicle);
  std::vector<CameraDetection> dets{
      detection_for(cam, {32.3, 1.0}, 0.0, ObjectClass::kVehicle, 0.6),
      detection_for(cam, {20.0, -6.0}, 0.0, ObjectClass::kPedestrian, 0.6)};
  const auto fused = fusion.fuse({trk}, 0.0, s, 0.0, dets);
  ASSERT_EQ(fused.size(), 1u);
  EXPECT_EQ(fused[0].camera_index, 0);
  EXPECT_GT(fused[0].iou, 0.5);
}

TEST(Fusion, TimeCompensationKeepsFastCrossingObjectAssociated) {
  CameraModel cam;
  cam.T_cam_from_radar = make_camera_extrinsic(Vec3(-1.8, 0.0, 0.8), 0.0, 0.0, 0.0);
  FusionConfig sync, nosync;
  nosync.compensate_time = false;
  RadarCameraFusion with(cam, sync), without(cam, nosync);
  SensorState s;
  // Track at radar time 0; camera frame 0.25 s later, object moved 6 m left
  // (more than a car length, so the stale box no longer overlaps).
  const auto trk = make_track(1, {20.0, -1.5}, {0.0, 24.0}, ObjectClass::kVehicle);
  const std::vector<CameraDetection> dets{
      detection_for(cam, {21.0, 4.5}, kPi / 2.0, ObjectClass::kVehicle, 0.6)};
  EXPECT_EQ(with.fuse({trk}, 0.0, s, 0.25, dets)[0].camera_index, 0);
  EXPECT_EQ(without.fuse({trk}, 0.0, s, 0.25, dets)[0].camera_index, -1);
}

TEST(Fusion, CameraEvidenceCorrectsRadarStaticPedestrianButNotOneBadFrame) {
  CameraModel cam;
  cam.T_cam_from_radar = make_camera_extrinsic(Vec3(-1.8, 0.0, 0.8), 0.0, 0.0, 0.0);
  RadarCameraFusion fusion(cam);
  SensorState s;
  const auto ped = make_track(1, {15.0, 2.0}, {0.0, 0.0}, ObjectClass::kStatic, 0.9);
  const auto car = make_track(2, {25.0, -3.0}, {8.0, 0.0}, ObjectClass::kVehicle, 0.9);
  const auto ped_det = detection_for(cam, {15.0, 2.0}, 0.0, ObjectClass::kPedestrian, 0.6);
  auto car_det = detection_for(cam, {27.3, -3.0}, 0.0, ObjectClass::kVehicle, 0.6);

  std::vector<FusedObject> out;
  for (int k = 0; k < 6; ++k) {
    // One frame where the camera mislabels the car.
    car_det.cls = k == 3 ? ObjectClass::kPedestrian : ObjectClass::kVehicle;
    out = fusion.fuse({ped, car}, 0.0, s, 0.0, {ped_det, car_det});
    for (const auto& f : out) {
      if (f.track_id == 2) {
        EXPECT_EQ(f.cls, ObjectClass::kVehicle) << "frame " << k;
      }
    }
  }
  for (const auto& f : out) {
    if (f.track_id == 1) {
      EXPECT_EQ(f.cls, ObjectClass::kPedestrian);
    }
  }
}

TEST(CameraSimulator, DetectsVehiclesAndPedestriansOnly) {
  CameraModel cam;
  cam.T_cam_from_radar = make_camera_extrinsic(Vec3(-1.8, 0.0, 0.8), 0.0, 0.0, 0.0);
  CameraDetectorConfig cfg;
  cfg.miss_probability = 0.0;
  cfg.false_positives_per_frame = 0.0;
  cfg.class_confusion = 0.0;
  CameraSimulator sim(cam, 0.6, cfg);
  Scene scene;
  scene.add(SceneObject::vehicle(1, {30.0, 2.0}, 0.0, 5.0));
  scene.add(SceneObject::pedestrian(2, {15.0, -3.0}, kPi / 2.0, 1.2));
  scene.add(SceneObject::pole(3, {20.0, 6.0}));
  const auto dets = sim.detect(scene.snapshot());
  ASSERT_EQ(dets.size(), 2u);
  for (const auto& d : dets) {
    EXPECT_NE(d.truth_id, 3);
    EXPECT_GT(d.box.area(), 100.0);
  }
}

}  // namespace radar
