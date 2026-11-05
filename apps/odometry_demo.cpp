// Ego-motion demo: radar Doppler ego-velocity (RANSAC) fused with a simulated
// IMU in an error-state Kalman filter, compared with IMU-only dead reckoning
// and radar-only odometry. Object tracking runs on the estimated ego state.
//
//   odometry_demo [--duration S] [--out DIR] [--seed N]

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "radar/clustering.hpp"
#include "radar/ego_estimator.hpp"
#include "radar/ego_velocity.hpp"
#include "radar/imu.hpp"
#include "radar/pipeline.hpp"
#include "radar/scenarios.hpp"
#include "radar/simulator.hpp"
#include "radar/tracker.hpp"
#include "radar/tracking_metrics.hpp"

namespace {

struct Args {
  double duration_s = 30.0;
  std::string out_dir = "out_odometry";
  std::uint32_t seed = 3;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + k);
      return argv[++i];
    };
    if (k == "--duration")
      a.duration_s = std::stod(next());
    else if (k == "--out")
      a.out_dir = next();
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::stoul(next()));
    else
      throw std::invalid_argument("unknown argument " + k);
  }
  return a;
}

struct ErrorStats {
  double sum_sq = 0.0;
  double final_err = 0.0;
  double max_err = 0.0;
  std::size_t n = 0;
  void add(double e) {
    sum_sq += e * e;
    final_err = e;
    max_err = std::max(max_err, e);
    ++n;
  }
  double rmse() const { return n ? std::sqrt(sum_sq / static_cast<double>(n)) : 0.0; }
};

double calibrate(const radar::RadarConfig& radar, radar::RadarPipeline& pipeline) {
  radar::SimulatorConfig sc;
  sc.seed = 999;
  radar::Simulator sim(radar, sc);
  radar::FrameResult res;
  double sum = 0.0;
  for (int i = 0; i < 10; ++i) {
    pipeline.process(sim.simulate_frame({{50.0, 0.0, 0.0, 10.0}}), res);
    double best = -1e9;
    for (const auto& d : res.detections) best = std::max(best, static_cast<double>(d.snr_db));
    sum += best;
  }
  return radar::calibrate_rcs_offset(sum / 10.0, 50.0, 10.0);
}

}  // namespace

int main(int argc, char** argv) {
  using namespace radar;
  try {
    const Args args = parse(argc, argv);
    Scene scene = scenarios::ego_drive({}, args.duration_s);
    const Pose2D mount = scene.config().sensor_mount;

    RadarConfig radar;
    SimulatorConfig sc;
    sc.seed = args.seed;
    Simulator sim(radar, sc);
    RadarPipeline pipeline(radar);
    pipeline.set_rcs_offset_db(calibrate(radar, pipeline));

    ImuConfig imu_cfg;
    imu_cfg.seed = args.seed + 10;
    ImuSimulator imu(imu_cfg);
    EgoEstimatorConfig ecfg;
    ecfg.sensor_mount = mount;
    EgoEstimator fused(ecfg), imu_only(ecfg);
    RadarOdometry radar_only(mount);
    EgoVelocityEstimator ego_vel;

    const EgoState& e0 = scene.ego();
    fused.initialize(e0.pose, e0.velocity(), 0.0);
    imu_only.initialize(e0.pose, e0.velocity(), 0.0);
    radar_only.initialize(e0.pose, 0.0);

    Dbscan dbscan;
    MultiTargetTracker trk_est, trk_truth;
    TrackingEvaluator eval_est, eval_truth;

    std::filesystem::create_directories(args.out_dir);
    std::ofstream f_traj(args.out_dir + "/trajectory.csv");
    f_traj << "time,true_x,true_y,true_yaw,fused_x,fused_y,fused_yaw,imu_x,imu_y,imu_yaw,"
              "radar_x,radar_y,radar_yaw,true_vsx,true_vsy,meas_vsx,meas_vsy,inliers,"
              "detections,accepted,gyro_bias_est,gyro_bias_true\n";

    const double imu_dt = 1.0 / imu_cfg.rate_hz;
    const int radar_every = 5;  // 20 Hz radar
    const auto steps = static_cast<int>(args.duration_s / imu_dt);
    ErrorStats err_fused, err_imu, err_radar, err_yaw_fused, err_yaw_imu, err_yaw_radar, err_vel;
    std::size_t frames = 0, valid = 0, accepted = 0;
    double inlier_frac_sum = 0.0;
    DataCube adc;
    FrameResult res;
    std::vector<Cluster> clusters;

    for (int k = 0; k <= steps; ++k) {
      const double t = k * imu_dt;
      const ImuSample s = imu.sample(scene.ego(), t);
      fused.propagate(s);
      imu_only.propagate(s);

      if (k % radar_every == 0) {
        const SceneFrame sf = scene.snapshot();
        sim.simulate_frame(sf.targets, adc);
        pipeline.process(adc, res);
        ++frames;

        const EgoVelocityEstimate ev = ego_vel.estimate(res.detections);
        const Vec2 true_vs = sf.sensor_pose.rotation_matrix().transpose() * sf.sensor_velocity;
        bool ok = false;
        if (ev.valid) {
          ++valid;
          inlier_frac_sum +=
              static_cast<double>(ev.inliers) / static_cast<double>(res.detections.size());
          err_vel.add((ev.velocity - true_vs).norm());
          ok = fused.update_radar_velocity(ev.velocity, ev.covariance);
          accepted += ok;
          radar_only.update(ev.velocity, t);
        }

        dbscan.cluster(res.detections, clusters);
        trk_est.update(t, clusters, fused.sensor_state());
        trk_truth.update(t, clusters, {sf.sensor_pose, sf.sensor_velocity});
        // Score the estimated-ego tracker relative to the sensor, so odometry
        // drift does not count as a tracking error.
        const Pose2D est_to_true = sf.sensor_pose.compose(fused.sensor_state().pose.inverse());
        eval_est.add_frame(sf.objects, trk_est.tracks(), est_to_true);
        eval_truth.add_frame(sf.objects, trk_truth.tracks());

        const Pose2D tp = sf.ego.pose, fp = fused.pose(), ip = imu_only.pose(),
                     rp = radar_only.pose();
        err_fused.add((fp.translation() - tp.translation()).norm());
        err_imu.add((ip.translation() - tp.translation()).norm());
        err_radar.add((rp.translation() - tp.translation()).norm());
        err_yaw_fused.add(std::abs(wrap_angle(fp.yaw - tp.yaw)));
        err_yaw_imu.add(std::abs(wrap_angle(ip.yaw - tp.yaw)));
        err_yaw_radar.add(std::abs(wrap_angle(rp.yaw - tp.yaw)));

        f_traj << t << ',' << tp.x << ',' << tp.y << ',' << tp.yaw << ',' << fp.x << ',' << fp.y
               << ',' << fp.yaw << ',' << ip.x << ',' << ip.y << ',' << ip.yaw << ',' << rp.x << ','
               << rp.y << ',' << rp.yaw << ',' << true_vs.x() << ',' << true_vs.y() << ','
               << (ev.valid ? ev.velocity.x() : NAN) << ',' << (ev.valid ? ev.velocity.y() : NAN)
               << ',' << ev.inliers << ',' << res.detections.size() << ',' << ok << ','
               << fused.gyro_bias() << ',' << imu.gyro_bias() << '\n';
      }
      scene.step(imu_dt);
    }

    const double deg = 180.0 / kPi;
    std::printf("Drive: %.0f s, %zu radar frames, %.0f m travelled\n", args.duration_s, frames,
                8.0 * args.duration_s + 0.15 * args.duration_s * args.duration_s);
    std::printf("\nRadar ego-velocity (RANSAC)\n");
    std::printf("  valid %zu / %zu frames, mean inlier fraction %.2f, fused updates accepted %zu\n",
                valid, frames, valid ? inlier_frac_sum / static_cast<double>(valid) : 0.0,
                accepted);
    std::printf("  velocity error RMSE %.3f m/s\n", err_vel.rmse());
    std::printf("\nPosition error            RMSE     final      max   | yaw final\n");
    std::printf("  IMU only            %8.2f m %8.2f m %8.2f m | %6.2f deg\n", err_imu.rmse(),
                err_imu.final_err, err_imu.max_err, err_yaw_imu.final_err * deg);
    std::printf("  radar only          %8.2f m %8.2f m %8.2f m | %6.2f deg\n", err_radar.rmse(),
                err_radar.final_err, err_radar.max_err, err_yaw_radar.final_err * deg);
    std::printf("  radar + IMU (ESKF)  %8.2f m %8.2f m %8.2f m | %6.2f deg\n", err_fused.rmse(),
                err_fused.final_err, err_fused.max_err, err_yaw_fused.final_err * deg);
    std::printf("\nGyro bias: true %.4f rad/s, estimated %.4f rad/s\n", imu.gyro_bias(),
                fused.gyro_bias());
    std::printf("Accel bias: true (%.3f, %.3f), estimated (%.3f, %.3f) m/s^2\n",
                imu.accel_bias().x(), imu.accel_bias().y(), fused.accel_bias().x(),
                fused.accel_bias().y());
    const auto& me = eval_est.metrics();
    const auto& mt = eval_truth.metrics();
    std::printf("\nTracking from the moving platform (scored relative to the sensor)\n");
    std::printf(
        "                                     MOTA   recall  false tracks  vel err  class acc\n");
    std::printf("  ego from radar + IMU estimate     %6.3f  %6.3f  %8zu   %5.2f m/s  %5.1f %%\n",
                me.mota(), me.recall(), me.false_positives, me.mean_velocity_error_mps(),
                100.0 * me.class_accuracy());
    std::printf("  ego from ground truth             %6.3f  %6.3f  %8zu   %5.2f m/s  %5.1f %%\n",
                mt.mota(), mt.recall(), mt.false_positives, mt.mean_velocity_error_mps(),
                100.0 * mt.class_accuracy());
    std::printf("Outputs written to %s/\n", args.out_dir.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
