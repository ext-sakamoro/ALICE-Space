// AliceSpace.cs — Unity P/Invoke bindings for alice-space
// Author: Moroya Sakamoto

using System;
using System.Runtime.InteropServices;

namespace Alice.Space
{
    public static class AliceSpaceNative
    {
#if UNITY_IOS && !UNITY_EDITOR
        private const string Lib = "__Internal";
#else
        private const string Lib = "alice_space";
#endif

        // Orbit
        [DllImport(Lib)] public static extern double alice_space_orbital_period(double semiMajorAxisKm, double mu);
        [DllImport(Lib)] public static extern double alice_space_orbital_velocity(double rKm, double aKm, double mu);
        [DllImport(Lib)] public static extern double alice_space_light_delay_s(double distanceKm);
        [DllImport(Lib)] public static extern void alice_space_delta_v_hohmann(double r1Km, double r2Km, double mu, double[] outDv);

        // LinkBudget
        [DllImport(Lib)] public static extern double alice_space_friis_path_loss_db(double distanceKm, double frequencyGhz);
        [DllImport(Lib)] public static extern double alice_space_link_budget_margin(double distanceKm, double frequencyGhz, double dataRateBps);

        // Propagator
        [DllImport(Lib)] public static extern void alice_space_rk4_single(double[] stateIn, double mu, double dtS, double[] stateOut);
        [DllImport(Lib)] public static extern void alice_space_two_body_accel(double[] r, double mu, double[] outAccel);

        // Comm
        [DllImport(Lib)] public static extern double alice_space_comm_latency_s(double distanceKm);
        [DllImport(Lib)] public static extern double alice_space_comm_bits_per_window(double bandwidthBps, double windowS);

        // Constellation
        [DllImport(Lib)] public static extern double alice_space_walker_plane_spacing_rad(uint numPlanes);
        [DllImport(Lib)] public static extern double alice_space_walker_ground_track_period_s(double semiMajorAxisKm, double mu);
    }
}
