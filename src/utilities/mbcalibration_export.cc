/*--------------------------------------------------------------------
 *    The MB-system:  mbcalibration_export.cc
 *
 *    Copyright (c) 2026 by
 *    Sebastian Rodriguez (srodriguez@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *
 *    See README.md file for copying and redistribution conditions.
 *--------------------------------------------------------------------
 *
 * mbcalibration_export reads MB-System processed .mb89 files together
 * with their ancillary files (.fnv, .baa) and exports per-ping calibration
 * data to a NetCDF4 file.
 *
 * Output schema (NetCDF4):
 *   Dimensions:  n_ping (unlimited), max_beam (fixed)
 *   Per-ping scalars (float64): ping_time, nav_lon, nav_lat, sonar_depth,
 *       altitude, heading, speed, roll, pitch, heave
 *   Per-ping int (int32): n_beams
 *   Per-beam (float64, [n_ping, max_beam]):
 *       bath_depth, bath_lon, bath_lat
 *   Per-beam amplitude (float32): amplitude
 *   Per-beam flags (uint8): beam_flag, beam_ok
 *
 * roll / pitch come from .fnv — corrected INS attitude AFTER all platform
 * transforms (lever arm + boresight).  This is the attitude MB-System used
 * to place bath_lon/bath_lat in the world frame, so it is the only attitude
 * that correctly inverts the transformation back to sensor frame.
 *
 * Author:   Sebastian Rodriguez
 * Date:     2026-05-08
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>
#include <limits>
#include <unistd.h>
#include <string>
#include <vector>

#include <netcdf.h>

#include "mb_define.h"
#include "mb_format.h"
#include "mb_io.h"
#include "mb_status.h"

/* ======================================================================
 * Ancillary file readers (formerly mb_ancillary_io.h/.cc)
 * Merged here to minimise library footprint — used only by this tool.
 * ====================================================================== */

struct MbFnvRecord {
    double time_d;       /* Unix epoch seconds */
    double navlon;       /* degrees East */
    double navlat;       /* degrees North */
    double heading;      /* degrees true */
    double speed;        /* km/hr */
    double sensordepth;  /* meters, positive down */
    double roll;         /* degrees */
    double pitch;        /* degrees */
    double heave;        /* meters */
};

/* .baa / .bsa: big-endian double time_d + float roll + float pitch (16 bytes) */
struct MbBaaRecord {
    double time_d;
    float  roll;
    float  pitch;
};

/* .bah: big-endian double time_d + float heading (12 bytes) */
struct MbBahRecord {
    double time_d;
    float  heading;
};

/* .bas: big-endian double time_d + float sensordepth (12 bytes) */
struct MbBasRecord {
    double time_d;
    float  sensordepth;
};

static bool anc_read_dff(FILE* fp, double* d, float* f1, float* f2) {
    unsigned char buf[16];
    if (fread(buf, 1, 16, fp) != 16) return false;
    mb_get_binary_double(true, buf,      d);
    mb_get_binary_float (true, buf + 8,  f1);
    mb_get_binary_float (true, buf + 12, f2);
    return true;
}

static bool anc_read_df(FILE* fp, double* d, float* f) {
    unsigned char buf[12];
    if (fread(buf, 1, 12, fp) != 12) return false;
    mb_get_binary_double(true, buf,     d);
    mb_get_binary_float (true, buf + 8, f);
    return true;
}

/* .fnv format: 19 columns per line (see mbpreprocess.cc:3428 for format string).
 * Columns 1-6 = calendar time (unused), 7=time_d, 8-9=lon/lat, 10=heading,
 * 11=speed, 12=sensordepth, 13=roll, 14=pitch, 15=heave, 16-19=swath limits. */
static std::vector<MbFnvRecord> mb_read_fnv(const std::string& path) {
    std::vector<MbFnvRecord> recs;
    FILE* fp = fopen(path.c_str(), "r");
    if (!fp) return recs;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        int yr, mo, dy, hr, mn; double sec;
        MbFnvRecord r; double p1, p2, p3, p4;
        int n = sscanf(line,
            "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &yr, &mo, &dy, &hr, &mn, &sec,
            &r.time_d, &r.navlon, &r.navlat, &r.heading, &r.speed,
            &r.sensordepth, &r.roll, &r.pitch, &r.heave, &p1, &p2, &p3, &p4);
        if (n >= 15) recs.push_back(r);
    }
    fclose(fp);
    return recs;
}

static std::vector<MbBaaRecord> mb_read_baa(const std::string& path) {
    std::vector<MbBaaRecord> recs;
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return recs;
    fseek(fp, 0, SEEK_END); long sz = ftell(fp); rewind(fp);
    if (sz > 0 && sz % 16 == 0) recs.reserve((size_t)(sz / 16));
    MbBaaRecord r;
    while (anc_read_dff(fp, &r.time_d, &r.roll, &r.pitch)) recs.push_back(r);
    fclose(fp);
    return recs;
}

static std::vector<MbBahRecord> mb_read_bah(const std::string& path) {
    std::vector<MbBahRecord> recs;
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return recs;
    fseek(fp, 0, SEEK_END); long sz = ftell(fp); rewind(fp);
    if (sz > 0 && sz % 12 == 0) recs.reserve((size_t)(sz / 12));
    MbBahRecord r;
    while (anc_read_df(fp, &r.time_d, &r.heading)) recs.push_back(r);
    fclose(fp);
    return recs;
}

static std::vector<MbBasRecord> mb_read_bas(const std::string& path) {
    std::vector<MbBasRecord> recs;
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return recs;
    fseek(fp, 0, SEEK_END); long sz = ftell(fp); rewind(fp);
    if (sz > 0 && sz % 12 == 0) recs.reserve((size_t)(sz / 12));
    MbBasRecord r;
    while (anc_read_df(fp, &r.time_d, &r.sensordepth)) recs.push_back(r);
    fclose(fp);
    return recs;
}

template <typename T>
static size_t anc_find_lower(const std::vector<T>& recs, double t) {
    size_t lo = 0, hi = recs.size() - 1;
    while (lo < hi) {
        size_t mid = lo + (hi - lo + 1) / 2;
        if (recs[mid].time_d <= t) lo = mid; else hi = mid - 1;
    }
    return lo;
}

static bool mb_interp_attitude(const std::vector<MbBaaRecord>& recs,
                               double t, double* roll, double* pitch) {
    if (recs.empty()) {
        *roll = *pitch = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    if (t <= recs.front().time_d) { *roll = recs.front().roll; *pitch = recs.front().pitch; return true; }
    if (t >= recs.back().time_d)  { *roll = recs.back().roll;  *pitch = recs.back().pitch;  return true; }
    size_t i = anc_find_lower(recs, t);
    double t0 = recs[i].time_d, t1 = recs[i+1].time_d;
    double f = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
    *roll  = recs[i].roll  + f * (recs[i+1].roll  - recs[i].roll);
    *pitch = recs[i].pitch + f * (recs[i+1].pitch - recs[i].pitch);
    return true;
}

static bool mb_interp_heading(const std::vector<MbBahRecord>& recs,
                              double t, double* heading) {
    if (recs.empty()) { *heading = std::numeric_limits<double>::quiet_NaN(); return false; }
    if (t <= recs.front().time_d) { *heading = recs.front().heading; return true; }
    if (t >= recs.back().time_d)  { *heading = recs.back().heading;  return true; }
    size_t i = anc_find_lower(recs, t);
    double t0 = recs[i].time_d, t1 = recs[i+1].time_d;
    double f = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
    double h0 = recs[i].heading, h1 = recs[i+1].heading, diff = h1 - h0;
    if (diff >  180.0) diff -= 360.0;
    if (diff < -180.0) diff += 360.0;
    double h = h0 + f * diff;
    if (h < 0.0) h += 360.0; else if (h >= 360.0) h -= 360.0;
    *heading = h;
    return true;
}

static bool mb_interp_sensordepth(const std::vector<MbBasRecord>& recs,
                                  double t, double* sd) {
    if (recs.empty()) { *sd = std::numeric_limits<double>::quiet_NaN(); return false; }
    if (t <= recs.front().time_d) { *sd = recs.front().sensordepth; return true; }
    if (t >= recs.back().time_d)  { *sd = recs.back().sensordepth;  return true; }
    size_t i = anc_find_lower(recs, t);
    double t0 = recs[i].time_d, t1 = recs[i+1].time_d;
    double f = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
    *sd = recs[i].sensordepth + f * (recs[i+1].sensordepth - recs[i].sensordepth);
    return true;
}

/* Look up corrected roll, pitch, and heave from .fnv records by nearest time_d.
 *
 * .fnv has exactly one record per survey ping, written with the same time_d
 * that is stored in the .mb89 ping record.  This is a 1:1 match, not an
 * interpolation problem — interpolating between adjacent ping records would
 * blend attitudes that belong to distinct pings and introduce geometric error.
 *
 * anc_find_lower() returns the largest index i where recs[i].time_d <= t.
 * We then pick whichever neighbor (i or i+1) is closer in time. */
static bool mb_interp_fnv_attitude(const std::vector<MbFnvRecord>& recs,
                                   double t,
                                   double* roll, double* pitch, double* heave) {
    if (recs.empty()) {
        *roll = *pitch = *heave = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    if (recs.size() == 1 || t <= recs.front().time_d) {
        *roll = recs.front().roll; *pitch = recs.front().pitch; *heave = recs.front().heave;
        return true;
    }
    if (t >= recs.back().time_d) {
        *roll = recs.back().roll; *pitch = recs.back().pitch; *heave = recs.back().heave;
        return true;
    }
    size_t i = anc_find_lower(recs, t);
    /* Pick the closer of i and i+1 */
    if (i + 1 < recs.size()) {
        double dt_lo = t - recs[i].time_d;
        double dt_hi = recs[i + 1].time_d - t;
        if (dt_hi < dt_lo) i += 1;
    }
    *roll  = recs[i].roll;
    *pitch = recs[i].pitch;
    *heave = recs[i].heave;
    return true;
}

/* End of ancillary file readers */


#define NC_ERR(call) do { \
    int _r = (call); \
    if (_r != NC_NOERR) { \
        fprintf(stderr, "NetCDF error in %s:%d: %s\n", \
                __FILE__, __LINE__, nc_strerror(_r)); \
        return _r; \
    } \
} while(0)

/* Extrinsics and preprocessing flags to embed as optimizer starting point.
 * Extrinsics = full 6-DOF sensor placement: [S, F, U meters, H, R, P degrees]
 * matching the --sensor-offsets=S/F/U/H/R/P convention of mbmakeplatform. */
struct ExtrinsicsInfo {
    char   sensor2_model[MB_LONGNAME_LENGTH] = {};
    /* TX extrinsics (offset 0): [S,F,U m, H,R,P deg] */
    double tx_extrinsics[6] = {0, 0, 0, 0, 0, 0};
    /* RX extrinsics (offset 1): [S,F,U m, H,R,P deg] */
    double rx_extrinsics[6] = {0, 0, 0, 0, 0, 0};
    bool   platform_loaded  = false;

    /* Kluge flags that have geometric impact on the data */
    bool   kluge_flipsign_pitch   = false;
    bool   kluge_flipsign_roll    = false;
    double kluge_soundspeed_tweak = 1.0;
    bool   attitude_zero_heave    = false;
};

/* Read sensor 2 extrinsics from a platform file.
 * target_sensor is the .plf sensor index for the multibeam (typically 2). */
static void load_platform_extrinsics(const char *plf_path, int target_sensor,
                                     ExtrinsicsInfo& ex) {
    int verbose = 0, error;
    void *platform_ptr = nullptr;
    if (mb_platform_read(verbose, (char *)plf_path, &platform_ptr, &error) != MB_SUCCESS
        || platform_ptr == nullptr)
        return;

    struct mb_platform_struct *platform = (struct mb_platform_struct *)platform_ptr;
    if (target_sensor < 0 || target_sensor >= platform->num_sensors) {
        mb_platform_deall(verbose, &platform_ptr, &error);
        return;
    }

    struct mb_sensor_struct *sensor = &platform->sensors[target_sensor];
    snprintf(ex.sensor2_model, sizeof(ex.sensor2_model), "%s", sensor->model);

    if (sensor->num_offsets >= 1) {
        struct mb_sensor_offset_struct *tx = &sensor->offsets[0];
        ex.tx_extrinsics[0] = tx->position_offset_x;      /* S */
        ex.tx_extrinsics[1] = tx->position_offset_y;      /* F */
        ex.tx_extrinsics[2] = tx->position_offset_z;      /* U */
        ex.tx_extrinsics[3] = tx->attitude_offset_heading; /* H */
        ex.tx_extrinsics[4] = tx->attitude_offset_roll;    /* R */
        ex.tx_extrinsics[5] = tx->attitude_offset_pitch;   /* P */
    }
    if (sensor->num_offsets >= 2) {
        struct mb_sensor_offset_struct *rx = &sensor->offsets[1];
        ex.rx_extrinsics[0] = rx->position_offset_x;      /* S */
        ex.rx_extrinsics[1] = rx->position_offset_y;      /* F */
        ex.rx_extrinsics[2] = rx->position_offset_z;      /* U */
        ex.rx_extrinsics[3] = rx->attitude_offset_heading; /* H */
        ex.rx_extrinsics[4] = rx->attitude_offset_roll;    /* R */
        ex.rx_extrinsics[5] = rx->attitude_offset_pitch;   /* P */
    }
    ex.platform_loaded = true;
    mb_platform_deall(verbose, &platform_ptr, &error);
}

struct CalWriter {
    int  ncid      = -1;
    int  ping_count = 0;
    int  max_beam  = 0;

    /* 1D variable IDs */
    int vid_ping_time, vid_nav_lon, vid_nav_lat, vid_sonar_depth;
    int vid_altitude, vid_heading, vid_speed;
    int vid_roll, vid_pitch, vid_heave, vid_n_beams;

    /* 2D variable IDs [n_ping, max_beam] */
    int vid_bath_depth, vid_bath_lon, vid_bath_lat;
    int vid_amplitude, vid_beam_flag, vid_beam_ok;

    /* scratch buffers for type conversion */
    std::vector<float>         amp_buf;
    std::vector<unsigned char> ok_buf;
};

/* Add a units attribute to a variable */
static int add_units(int ncid, int varid, const char* units) {
    return nc_put_att_text(ncid, varid, "units", strlen(units), units);
}

/* Add deflate + chunking to a 2D variable */
static int set_compression(int ncid, int varid, int max_beam) {
    size_t chunks[2] = {1, (size_t)max_beam};  /* one ping at a time */
    nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunks);
    return nc_def_var_deflate(ncid, varid, 1 /*shuffle*/, 1, 4 /*level*/);
}

static int calwriter_open(CalWriter& w, const std::string& path,
                           int max_beam, const std::string& platform_file,
                           const std::string& command_line,
                           const ExtrinsicsInfo& ex) {
    w.max_beam = max_beam;
    w.amp_buf.resize(max_beam);
    w.ok_buf.resize(max_beam);

    NC_ERR(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &w.ncid));

    /* --- Dimensions --- */
    int n_ping_dim, beam_dim;
    NC_ERR(nc_def_dim(w.ncid, "n_ping",   NC_UNLIMITED, &n_ping_dim));
    NC_ERR(nc_def_dim(w.ncid, "max_beam", (size_t)max_beam, &beam_dim));

    int dims1[1] = {n_ping_dim};
    int dims2[2] = {n_ping_dim, beam_dim};

    /* --- 1D per-ping variables --- */
    NC_ERR(nc_def_var(w.ncid, "ping_time",   NC_DOUBLE, 1, dims1, &w.vid_ping_time));
    NC_ERR(nc_def_var(w.ncid, "nav_lon",     NC_DOUBLE, 1, dims1, &w.vid_nav_lon));
    NC_ERR(nc_def_var(w.ncid, "nav_lat",     NC_DOUBLE, 1, dims1, &w.vid_nav_lat));
    NC_ERR(nc_def_var(w.ncid, "sonar_depth", NC_DOUBLE, 1, dims1, &w.vid_sonar_depth));
    NC_ERR(nc_def_var(w.ncid, "altitude",    NC_DOUBLE, 1, dims1, &w.vid_altitude));
    NC_ERR(nc_def_var(w.ncid, "heading",     NC_DOUBLE, 1, dims1, &w.vid_heading));
    NC_ERR(nc_def_var(w.ncid, "speed",       NC_DOUBLE, 1, dims1, &w.vid_speed));
    NC_ERR(nc_def_var(w.ncid, "roll",    NC_DOUBLE, 1, dims1, &w.vid_roll));
    NC_ERR(nc_def_var(w.ncid, "pitch",   NC_DOUBLE, 1, dims1, &w.vid_pitch));
    NC_ERR(nc_def_var(w.ncid, "heave",       NC_DOUBLE, 1, dims1, &w.vid_heave));
    NC_ERR(nc_def_var(w.ncid, "n_beams",     NC_INT,    1, dims1, &w.vid_n_beams));

    /* Units for 1D variables */
    add_units(w.ncid, w.vid_ping_time,   "seconds since 1970-01-01T00:00:00Z");
    add_units(w.ncid, w.vid_nav_lon,     "degrees_east");
    add_units(w.ncid, w.vid_nav_lat,     "degrees_north");
    add_units(w.ncid, w.vid_sonar_depth, "m");
    add_units(w.ncid, w.vid_altitude,    "m");
    add_units(w.ncid, w.vid_heading,     "degrees_true");
    add_units(w.ncid, w.vid_speed,       "km/hr");
    add_units(w.ncid, w.vid_roll,    "degrees");
    add_units(w.ncid, w.vid_pitch,   "degrees");
    add_units(w.ncid, w.vid_heave,       "m");

    /* Long-name hints for the optimizer */
    nc_put_att_text(w.ncid, w.vid_roll,  "long_name",
                    strlen("corrected roll from .fnv (after lever arm and boresight)"),
                             "corrected roll from .fnv (after lever arm and boresight)");
    nc_put_att_text(w.ncid, w.vid_pitch, "long_name",
                    strlen("corrected pitch from .fnv (after lever arm and boresight)"),
                             "corrected pitch from .fnv (after lever arm and boresight)");

    /* --- 2D per-beam variables ---
     * bath_acrosstrack / bath_alongtrack omitted: mb_read fills bathlon/bathlat
     * (world-frame) but not the relative acrosstrack/alongtrack distances.
     * World-frame positions are sufficient for the calibration optimizer. */
    NC_ERR(nc_def_var(w.ncid, "bath_depth", NC_DOUBLE, 2, dims2, &w.vid_bath_depth));
    NC_ERR(nc_def_var(w.ncid, "bath_lon",   NC_DOUBLE, 2, dims2, &w.vid_bath_lon));
    NC_ERR(nc_def_var(w.ncid, "bath_lat",   NC_DOUBLE, 2, dims2, &w.vid_bath_lat));
    NC_ERR(nc_def_var(w.ncid, "amplitude",  NC_FLOAT,  2, dims2, &w.vid_amplitude));
    NC_ERR(nc_def_var(w.ncid, "beam_flag",  NC_UBYTE,  2, dims2, &w.vid_beam_flag));
    NC_ERR(nc_def_var(w.ncid, "beam_ok",    NC_UBYTE,  2, dims2, &w.vid_beam_ok));

    add_units(w.ncid, w.vid_bath_depth, "m");
    add_units(w.ncid, w.vid_bath_lon,   "degrees_east");
    add_units(w.ncid, w.vid_bath_lat,   "degrees_north");
    add_units(w.ncid, w.vid_amplitude,  "dB");

    /* Compression on 2D variables */
    set_compression(w.ncid, w.vid_bath_depth, max_beam);
    set_compression(w.ncid, w.vid_bath_lon,   max_beam);
    set_compression(w.ncid, w.vid_bath_lat,   max_beam);
    set_compression(w.ncid, w.vid_amplitude,  max_beam);
    set_compression(w.ncid, w.vid_beam_flag,  max_beam);
    set_compression(w.ncid, w.vid_beam_ok,    max_beam);

    /* --- Global attributes --- */
    time_t now = time(nullptr);
    char tstr[64];
    strftime(tstr, sizeof(tstr), "%Y-%m-%dT%H:%M:%SZ", gmtime(&now));

    nc_put_att_text(w.ncid, NC_GLOBAL, "mbsystem_version",
                    strlen(MB_VERSION), MB_VERSION);
    nc_put_att_text(w.ncid, NC_GLOBAL, "creation_time",
                    strlen(tstr), tstr);
    nc_put_att_text(w.ncid, NC_GLOBAL, "command_line",
                    command_line.size(), command_line.c_str());
    nc_put_att_text(w.ncid, NC_GLOBAL, "attitude_source",
                    strlen(".fnv ancillary file"), ".fnv ancillary file");
    nc_put_att_text(w.ncid, NC_GLOBAL, "frame_convention",
                    strlen("starboard-forward-up"), "starboard-forward-up");
    nc_put_att_text(w.ncid, NC_GLOBAL, "pitch_sign_convention",
                    strlen("bow-up-positive"), "bow-up-positive");
    nc_put_att_text(w.ncid, NC_GLOBAL, "roll_pitch_source",
                    strlen("corrected attitude after platform transforms (from .fnv ancillary file)"),
                             "corrected attitude after platform transforms (from .fnv ancillary file)");
    int target_sensor = 2;
    nc_put_att_int(w.ncid, NC_GLOBAL, "platform_target_sensor", NC_INT, 1, &target_sensor);

    /* --- Platform file and sensor 2 extrinsics (optimizer starting point) ---
     * Extrinsics = [S, F, U meters, H, R, P degrees] per array element,
     * matching the mbmakeplatform --sensor-offsets=S/F/U/H/R/P convention. */
    nc_put_att_text(w.ncid, NC_GLOBAL, "platform_file",
                    platform_file.size(), platform_file.c_str());
    if (ex.platform_loaded) {
        nc_put_att_text(w.ncid, NC_GLOBAL, "sensor2_model",
                        strlen(ex.sensor2_model), ex.sensor2_model);
        /* TX extrinsics (offset 0): [S, F, U m, H, R, P deg] */
        nc_put_att_double(w.ncid, NC_GLOBAL, "sensor2_tx_extrinsics",
                          NC_DOUBLE, 6, ex.tx_extrinsics);
        /* RX extrinsics (offset 1): [S, F, U m, H, R, P deg] */
        nc_put_att_double(w.ncid, NC_GLOBAL, "sensor2_rx_extrinsics",
                          NC_DOUBLE, 6, ex.rx_extrinsics);
    }

    /* --- Kluge flags that affect beam geometry / depths --- */
    int v;
    v = ex.kluge_flipsign_pitch ? 1 : 0;
    nc_put_att_int(w.ncid, NC_GLOBAL, "kluge_flipsign_pitch", NC_INT, 1, &v);
    v = ex.kluge_flipsign_roll ? 1 : 0;
    nc_put_att_int(w.ncid, NC_GLOBAL, "kluge_flipsign_roll",  NC_INT, 1, &v);
    nc_put_att_double(w.ncid, NC_GLOBAL, "kluge_soundspeed_tweak",
                      NC_DOUBLE, 1, &ex.kluge_soundspeed_tweak);
    v = ex.attitude_zero_heave ? 1 : 0;
    nc_put_att_int(w.ncid, NC_GLOBAL, "attitude_zero_heave",  NC_INT, 1, &v);

    NC_ERR(nc_enddef(w.ncid));
    return NC_NOERR;
}

static int calwriter_add_ping(CalWriter& w,
    double time_d, double nav_lon, double nav_lat,
    double sonar_depth, double altitude, double heading, double speed,
    double roll, double pitch, double heave,
    int n_beams,
    const double *bath_depth, const double *bath_lon,
    const double *bath_lat,   const double *amp_d,
    const char   *beam_flag) {

    size_t idx = (size_t)w.ping_count;
    size_t start1[1] = {idx};
    size_t count1[1] = {1};
    size_t start2[2] = {idx, 0};
    size_t count2[2] = {1, (size_t)w.max_beam};

    /* 1D scalars */
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_ping_time,   start1, count1, &time_d));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_nav_lon,     start1, count1, &nav_lon));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_nav_lat,     start1, count1, &nav_lat));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_sonar_depth, start1, count1, &sonar_depth));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_altitude,    start1, count1, &altitude));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_heading,     start1, count1, &heading));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_speed,       start1, count1, &speed));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_roll,    start1, count1, &roll));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_pitch,   start1, count1, &pitch));
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_heave,       start1, count1, &heave));
    NC_ERR(nc_put_vara_int   (w.ncid, w.vid_n_beams,     start1, count1, &n_beams));

    /* 2D beam arrays: [0, n_beams-1] actual data, [n_beams, max_beam-1] fill. */
    const int mb = w.max_beam;
    int n = std::min(n_beams, mb);

    std::vector<double> dbuf(mb, NC_FILL_DOUBLE);
    std::fill(w.amp_buf.begin(), w.amp_buf.end(), NC_FILL_FLOAT);
    std::fill(w.ok_buf.begin(),  w.ok_buf.end(),  0);

    /* bath_depth */
    std::copy(bath_depth, bath_depth + n, dbuf.begin());
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_bath_depth, start2, count2, dbuf.data()));

    /* bath_lon */
    std::copy(bath_lon, bath_lon + n, dbuf.begin());
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_bath_lon, start2, count2, dbuf.data()));

    /* bath_lat */
    std::copy(bath_lat, bath_lat + n, dbuf.begin());
    NC_ERR(nc_put_vara_double(w.ncid, w.vid_bath_lat, start2, count2, dbuf.data()));

    /* amplitude (float32) */
    for (int ib = 0; ib < n; ib++)
        w.amp_buf[ib] = (float)amp_d[ib];
    NC_ERR(nc_put_vara_float(w.ncid, w.vid_amplitude, start2, count2, w.amp_buf.data()));

    /* beam_flag and beam_ok (uint8) */
    std::vector<unsigned char> flag_buf(mb, 0);
    for (int ib = 0; ib < n; ib++) {
        flag_buf[ib]  = (unsigned char)beam_flag[ib];
        w.ok_buf[ib]  = mb_beam_ok(beam_flag[ib]) ? 1 : 0;
    }
    NC_ERR(nc_put_vara_uchar(w.ncid, w.vid_beam_flag, start2, count2, flag_buf.data()));
    NC_ERR(nc_put_vara_uchar(w.ncid, w.vid_beam_ok,   start2, count2, w.ok_buf.data()));

    w.ping_count++;
    return NC_NOERR;
}

static void calwriter_close(CalWriter& w) {
    if (w.ncid >= 0) {
        nc_close(w.ncid);
        w.ncid = -1;
    }
}

/* ======================================================================
 * CLI
 * ====================================================================== */

constexpr char program_name[] = "mbcalibration_export";

constexpr char help_message[] =
    "mbcalibration_export exports per-ping calibration data from processed\n"
    "MB-System .mb89 files to NetCDF4 for automated extrinsics\n"
    "calibration. Ancillary files (.fnv, .baa) must exist alongside the .mb89\n"
    "files (produced by mbpreprocess --output-sensor-fnv).";

constexpr char usage_message[] =
    "mbcalibration_export -I datalist -O output.nc\n"
    "                     [--platform-file=path]\n"
    "                     [--kluge-flipsign-pitch]\n"
    "                     [--kluge-flipsign-roll]\n"
    "                     [--kluge-soundspeed-tweak=factor]\n"
    "                     [--attitude-zero-heave]\n"
    "                     [--include-raw] [-V] [-H]\n"
    "\n"
    "Kluge flags must match those used in the mbpreprocess command that\n"
    "produced the .mb89 files — they are recorded as metadata for the optimizer.";

static struct option long_options[] = {
    {"platform-file",          required_argument, nullptr, 0},
    {"include-raw",            no_argument,       nullptr, 0},
    {"kluge-flipsign-pitch",   no_argument,       nullptr, 0},
    {"kluge-flipsign-roll",    no_argument,       nullptr, 0},
    {"kluge-soundspeed-tweak", required_argument, nullptr, 0},
    {"attitude-zero-heave",    no_argument,       nullptr, 0},
    {"help",                   no_argument,       nullptr, 'H'},
    {"verbose",                no_argument,       nullptr, 'V'},
    {nullptr,                  0,                 nullptr, 0}
};

/* Build a command-line string from argv for metadata */
static std::string reconstruct_cmdline(int argc, char **argv) {
    std::string s;
    for (int i = 0; i < argc; i++) {
        if (i > 0) s += ' ';
        s += argv[i];
    }
    return s;
}

/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {
    int verbose = 0;
    int error = MB_ERROR_NO_ERROR;

    /* MB-System defaults */
    int format;
    int pings;
    int lonflip;
    double bounds[4];
    int btime_i[7];
    int etime_i[7];
    double speedmin;
    double timegap;
    mb_defaults(verbose, &format, &pings, &lonflip, bounds,
                btime_i, etime_i, &speedmin, &timegap);

    /* Command-line parameters */
    char read_file[MB_PATH_MAXLINE]     = "";
    char output_file[MB_PATH_MAXLINE]   = "";
    char platform_file[MB_PATH_MAXLINE] = "";
    bool include_raw = false;
    ExtrinsicsInfo ex;

    /* Parse arguments */
    {
        bool help = false;
        int c;
        int option_index = 0;
        while ((c = getopt_long(argc, argv, "HhI:i:O:o:Vv",
                                long_options, &option_index)) != -1) {
            switch (c) {
            case 0: {
                const char *opt = long_options[option_index].name;
                if (strcmp(opt, "platform-file") == 0)
                    snprintf(platform_file, sizeof(platform_file), "%s", optarg);
                else if (strcmp(opt, "include-raw") == 0)
                    include_raw = true;
                else if (strcmp(opt, "kluge-flipsign-pitch") == 0)
                    ex.kluge_flipsign_pitch = true;
                else if (strcmp(opt, "kluge-flipsign-roll") == 0)
                    ex.kluge_flipsign_roll = true;
                else if (strcmp(opt, "kluge-soundspeed-tweak") == 0)
                    sscanf(optarg, "%lf", &ex.kluge_soundspeed_tweak);
                else if (strcmp(opt, "attitude-zero-heave") == 0)
                    ex.attitude_zero_heave = true;
                break;
            }
            case 'H': case 'h':
                help = true;
                break;
            case 'I': case 'i':
                snprintf(read_file, sizeof(read_file), "%s", optarg);
                break;
            case 'O': case 'o':
                snprintf(output_file, sizeof(output_file), "%s", optarg);
                break;
            case 'V': case 'v':
                verbose++;
                break;
            default:
                fprintf(stderr, "\nOption -%c not recognized\n", c);
                fprintf(stderr, "\nUsage: %s\n", usage_message);
                exit(MB_ERROR_BAD_PARAMETER);
            }
        }

        if (help) {
            fprintf(stderr, "\nProgram %s\nMB-system Version %s\n\n",
                    program_name, MB_VERSION);
            fprintf(stderr, "%s\n\nUsage: %s\n\n", help_message, usage_message);
            exit(MB_ERROR_NO_ERROR);
        }
    }

    /* Validate required arguments */
    if (strlen(read_file) == 0) {
        fprintf(stderr, "\nProgram %s: input datalist required (-I)\n", program_name);
        fprintf(stderr, "Usage: %s\n\n", usage_message);
        exit(MB_ERROR_BAD_PARAMETER);
    }
    if (strlen(output_file) == 0) {
        fprintf(stderr, "\nProgram %s: output file required (-O)\n", program_name);
        fprintf(stderr, "Usage: %s\n\n", usage_message);
        exit(MB_ERROR_BAD_PARAMETER);
    }

    /* Load platform file extrinsics if provided */
    if (strlen(platform_file) > 0) {
        load_platform_extrinsics(platform_file, 2, ex);
        if (!ex.platform_loaded)
            fprintf(stderr, "Warning: could not read platform file %s\n",
                    platform_file);
    }

    if (verbose > 0) {
        fprintf(stderr, "\nProgram %s\nMB-system Version %s\n",
                program_name, MB_VERSION);
        fprintf(stderr, "\nInput:           %s\n", read_file);
        fprintf(stderr, "Output:          %s\n", output_file);
        fprintf(stderr, "Platform file:   %s%s\n",
                strlen(platform_file) > 0 ? platform_file : "(none)",
                ex.platform_loaded ? " (loaded)" : "");
        if (ex.platform_loaded) {
            fprintf(stderr, "  Sensor 2 (%s) TX extrinsics S/F/U m: %.4f %.4f %.4f  H/R/P deg: %.4f %.4f %.4f\n",
                    ex.sensor2_model,
                    ex.tx_extrinsics[0], ex.tx_extrinsics[1], ex.tx_extrinsics[2],
                    ex.tx_extrinsics[3], ex.tx_extrinsics[4], ex.tx_extrinsics[5]);
            fprintf(stderr, "  Sensor 2 RX extrinsics S/F/U m: %.4f %.4f %.4f  H/R/P deg: %.4f %.4f %.4f\n",
                    ex.rx_extrinsics[0], ex.rx_extrinsics[1], ex.rx_extrinsics[2],
                    ex.rx_extrinsics[3], ex.rx_extrinsics[4], ex.rx_extrinsics[5]);
        }
        fprintf(stderr, "Kluge pitch flip: %s  roll flip: %s  ss_tweak: %.6f  zero_heave: %s\n\n",
                ex.kluge_flipsign_pitch ? "yes" : "no",
                ex.kluge_flipsign_roll  ? "yes" : "no",
                ex.kluge_soundspeed_tweak,
                ex.attitude_zero_heave  ? "yes" : "no");
    }

    /* Determine input type */
    if (format == 0)
        mb_get_format(verbose, read_file, nullptr, &format, &error);
    const bool read_datalist = format < 0;

    char file[MB_PATH_MAXLINE];
    char dfile[MB_PATH_MAXLINE];
    void *datalist = nullptr;
    double file_weight;
    bool read_data;

    if (read_datalist) {
        const int look_processed = MB_DATALIST_LOOK_UNSET;
        if (mb_datalist_open(verbose, &datalist, read_file,
                             look_processed, &error) != MB_SUCCESS) {
            fprintf(stderr, "\nUnable to open data list file: %s\n", read_file);
            exit(MB_ERROR_OPEN_FAIL);
        }
        read_data = mb_datalist_read(verbose, datalist, file, dfile,
                                     &format, &file_weight, &error) == MB_SUCCESS;
    } else {
        snprintf(file, sizeof(file), "%s", read_file);
        read_data = true;
    }

    /* Beam arrays — allocated by MB-System memory manager */
    char   *beamflag       = nullptr;
    double *bath           = nullptr;
    double *bathacrosstrack = nullptr;
    double *bathalongtrack  = nullptr;
    double *bathlon        = nullptr;
    double *bathlat        = nullptr;
    double *amp            = nullptr;
    double *ss             = nullptr;
    double *ssacrosstrack  = nullptr;
    double *ssalongtrack   = nullptr;

    /* Initialise the NetCDF writer — opened on first file so we know max_beam */
    CalWriter writer;
    bool writer_open = false;

    int n_files = 0;
    int n_pings_total  = 0;
    int n_beams_total  = 0;

    /* ----------------------------------------------------------------
     * Main file loop
     * ---------------------------------------------------------------- */
    while (read_data) {
        n_files++;

        std::string fnv_path = std::string(file) + ".fnv";
        std::string baa_path = std::string(file) + ".baa";

        auto fnv_recs = mb_read_fnv(fnv_path);
        auto baa_recs = mb_read_baa(baa_path);

        if (fnv_recs.empty()) {
            fprintf(stderr, "Warning: no .fnv records for %s — skipping\n", file);
            goto next_file;
        }
        if (baa_recs.empty() && verbose > 0)
            fprintf(stderr, "  Note: no .baa records for %s (not needed for corrected attitude)\n", file);

        if (verbose > 0) {
            fprintf(stderr, "File %d: %s\n", n_files, file);
            fprintf(stderr, "  .fnv pings: %zu\n", fnv_recs.size());
        }

        {
            double btime_d, etime_d;
            int beams_bath, beams_amp, pixels_ss;
            void *mbio_ptr = nullptr;

            if (mb_read_init(verbose, file, format, pings, lonflip, bounds,
                             btime_i, etime_i, speedmin, timegap, &mbio_ptr,
                             &btime_d, &etime_d,
                             &beams_bath, &beams_amp, &pixels_ss,
                             &error) != MB_SUCCESS) {
                char *message;
                mb_error(verbose, error, &message);
                fprintf(stderr, "Warning: cannot open %s: %s\n", file, message);
                goto next_file;
            }

            /* Open NetCDF output file using max_beam from the first file's init */
            if (!writer_open) {
                std::string cmdline = reconstruct_cmdline(argc, argv);
                if (calwriter_open(writer, output_file, beams_bath,
                                   platform_file, cmdline, ex) != NC_NOERR) {
                    fprintf(stderr, "Fatal: cannot create output file %s\n",
                            output_file);
                    mb_close(verbose, &mbio_ptr, &error);
                    exit(1);
                }
                writer_open = true;
                if (verbose > 0)
                    fprintf(stderr, "Output file: %s  (max_beam=%d)\n",
                            output_file, beams_bath);
            }

            /* Register beam arrays */
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(char),   (void **)&beamflag, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bath, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathacrosstrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathalongtrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathlon, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathlat, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE,
                              sizeof(double), (void **)&amp, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ss, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ssacrosstrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ssalongtrack, &error);

            if (error != MB_ERROR_NO_ERROR) {
                fprintf(stderr, "Warning: array allocation failed for %s — skipping\n",
                        file);
                mb_close(verbose, &mbio_ptr, &error);
                goto next_file;
            }

            /* ---- Ping loop (mb_read returns only survey pings + bathlon/bathlat) ---- */
            int n_pings_file = 0;
            int kind;
            int rpings;
            int time_i[7];
            double time_d, navlon, navlat, speed, heading;
            double distance, altitude, sensordepth;
            char comment[MB_COMMENT_MAXLINE];

            while (error <= MB_ERROR_NO_ERROR) {
                int status = mb_read(verbose, mbio_ptr, &kind, &rpings,
                                     time_i, &time_d, &navlon, &navlat,
                                     &speed, &heading,
                                     &distance, &altitude, &sensordepth,
                                     &beams_bath, &beams_amp, &pixels_ss,
                                     beamflag, bath, amp,
                                     bathlon, bathlat,
                                     ss, ssacrosstrack, ssalongtrack,
                                     comment, &error);

                if (error == MB_ERROR_TIME_GAP) {
                    error = MB_ERROR_NO_ERROR;
                    status = MB_SUCCESS;
                }
                if (status != MB_SUCCESS)
                    continue;

                /* Get corrected attitude and heave from .fnv.
                 * These are the values MB-System used to place bath_lon/bath_lat
                 * in the world frame — required to correctly invert Tsw. */
                double fnv_roll, fnv_pitch, heave;
                mb_interp_fnv_attitude(fnv_recs, time_d, &fnv_roll, &fnv_pitch, &heave);

                /* Write to NetCDF */
                calwriter_add_ping(writer,
                    time_d, navlon, navlat, sensordepth, altitude,
                    heading, speed,
                    fnv_roll, fnv_pitch, heave,
                    beams_bath,
                    bath, bathlon, bathlat,
                    amp, beamflag);

                n_pings_file++;

                for (int ib = 0; ib < beams_bath; ib++)
                    if (mb_beam_ok(beamflag[ib]))
                        n_beams_total++;
            }

            n_pings_total += n_pings_file;

            if (verbose > 0)
                fprintf(stderr, "  Survey pings written: %d\n", n_pings_file);

            mb_close(verbose, &mbio_ptr, &error);
            error = MB_ERROR_NO_ERROR;
        }

next_file:
        if (read_datalist)
            read_data = mb_datalist_read(verbose, datalist, file, dfile,
                                         &format, &file_weight, &error) == MB_SUCCESS;
        else
            read_data = false;
    }

    if (read_datalist)
        mb_datalist_close(verbose, &datalist, &error);

    calwriter_close(writer);

    fprintf(stderr, "\nProgram %s completed:\n", program_name);
    fprintf(stderr, "  Files processed:    %d\n", n_files);
    fprintf(stderr, "  Survey pings:       %d\n", n_pings_total);
    fprintf(stderr, "  Good beam returns:  %d\n", n_beams_total);
    fprintf(stderr, "  Output:             %s\n", output_file);

    return MB_ERROR_NO_ERROR;
}
