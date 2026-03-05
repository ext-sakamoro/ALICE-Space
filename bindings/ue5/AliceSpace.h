// AliceSpace.h — UE5 C FFI header for alice-space
// Author: Moroya Sakamoto

#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Orbit
double alice_space_orbital_period(double semi_major_axis_km, double mu);
double alice_space_orbital_velocity(double r_km, double a_km, double mu);
double alice_space_light_delay_s(double distance_km);
void   alice_space_delta_v_hohmann(double r1_km, double r2_km, double mu, double* out_dv);

// LinkBudget
double alice_space_friis_path_loss_db(double distance_km, double frequency_ghz);
double alice_space_link_budget_margin(double distance_km, double frequency_ghz, double data_rate_bps);

// Propagator
void   alice_space_rk4_single(const double* state_in, double mu, double dt_s, double* state_out);
void   alice_space_two_body_accel(const double* r, double mu, double* out_accel);

// Comm
double alice_space_comm_latency_s(double distance_km);
double alice_space_comm_bits_per_window(double bandwidth_bps, double window_s);

// Constellation
double alice_space_walker_plane_spacing_rad(uint32_t num_planes);
double alice_space_walker_ground_track_period_s(double semi_major_axis_km, double mu);

#ifdef __cplusplus
}
#endif
