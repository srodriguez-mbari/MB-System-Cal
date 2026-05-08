/*--------------------------------------------------------------------
 *    The MB-system:  mb_ancillary_io.cc
 *
 *    Copyright (c) 2026 by
 *    Sebastian Rodriguez (seroma09@gmail.com)
 *      Monterey Bay Aquarium Research Institute
 *
 *    See README.md file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include "mb_ancillary_io.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>

/* MB-System binary byte-swap helpers */
extern "C" {
#include "mb_define.h"
}

/* --------------------------------------------------------------------------
 * Internal helpers
 * -------------------------------------------------------------------------- */

/* Read one big-endian double + float + float record (16 bytes) from fp.
 * Returns false on short read. */
static bool read_double_float_float(FILE* fp, double* d, float* f1, float* f2) {
    unsigned char buf[16];
    if (fread(buf, 1, 16, fp) != 16)
        return false;
    mb_get_binary_double(true, buf,     d);
    mb_get_binary_float (true, buf + 8, f1);
    mb_get_binary_float (true, buf + 12, f2);
    return true;
}

/* Read one big-endian double + float record (12 bytes) from fp. */
static bool read_double_float(FILE* fp, double* d, float* f) {
    unsigned char buf[12];
    if (fread(buf, 1, 12, fp) != 12)
        return false;
    mb_get_binary_double(true, buf,    d);
    mb_get_binary_float (true, buf + 8, f);
    return true;
}

/* --------------------------------------------------------------------------
 * .fnv reader
 * --------------------------------------------------------------------------
 * Format (one data line per survey ping, confirmed from mbpreprocess.cc:3428):
 *   col 1-6  : yyyy mm dd hh mm ss.ssssss  (calendar time, unused)
 *   col 7    : time_d        (Unix epoch seconds, tab-separated)
 *   col 8-9  : navlon navlat (degrees)
 *   col 10   : heading       (degrees)
 *   col 11   : speed         (km/hr)
 *   col 12   : sensordepth   (meters, positive down; header labels it "draft")
 *   col 13   : roll          (degrees)
 *   col 14   : pitch         (degrees)
 *   col 15   : heave         (meters)
 *   col 16-19: portlon portlat stbdlon stbdlat (skipped)
 * Header/comment lines start with '#'. */
std::vector<MbFnvRecord> mb_read_fnv(const std::string& path) {
    std::vector<MbFnvRecord> records;

    FILE* fp = fopen(path.c_str(), "r");
    if (!fp)
        return records;

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;

        int year, month, day, hour, minute;
        double seconds;
        MbFnvRecord r;
        double portlon, portlat, stbdlon, stbdlat;

        int n = sscanf(line,
            "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &year, &month, &day, &hour, &minute, &seconds,
            &r.time_d, &r.navlon, &r.navlat,
            &r.heading, &r.speed, &r.sensordepth,
            &r.roll, &r.pitch, &r.heave,
            &portlon, &portlat, &stbdlon, &stbdlat);

        if (n < 15) {
            /* tolerate missing portlon/portlat/stbdlon/stbdlat */
            fprintf(stderr, "mb_read_fnv: short parse (%d fields) in %s\n",
                    n, path.c_str());
            continue;
        }
        records.push_back(r);
    }

    fclose(fp);
    return records;
}

/* --------------------------------------------------------------------------
 * .baa reader (async attitude, big-endian: double + float + float = 16 bytes)
 * Same layout as .bsa (sync attitude per ping). */
std::vector<MbBaaRecord> mb_read_baa(const std::string& path) {
    std::vector<MbBaaRecord> records;

    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp)
        return records;

    /* Pre-size the vector based on file size */
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);
    if (fsize > 0 && (fsize % 16) == 0)
        records.reserve(static_cast<size_t>(fsize / 16));

    MbBaaRecord r;
    while (read_double_float_float(fp, &r.time_d, &r.roll, &r.pitch))
        records.push_back(r);

    fclose(fp);
    return records;
}

/* --------------------------------------------------------------------------
 * .bah reader (async heading, big-endian: double + float = 12 bytes) */
std::vector<MbBahRecord> mb_read_bah(const std::string& path) {
    std::vector<MbBahRecord> records;

    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp)
        return records;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);
    if (fsize > 0 && (fsize % 12) == 0)
        records.reserve(static_cast<size_t>(fsize / 12));

    MbBahRecord r;
    while (read_double_float(fp, &r.time_d, &r.heading))
        records.push_back(r);

    fclose(fp);
    return records;
}

/* --------------------------------------------------------------------------
 * .bas reader (async sensor depth, big-endian: double + float = 12 bytes) */
std::vector<MbBasRecord> mb_read_bas(const std::string& path) {
    std::vector<MbBasRecord> records;

    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp)
        return records;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);
    if (fsize > 0 && (fsize % 12) == 0)
        records.reserve(static_cast<size_t>(fsize / 12));

    MbBasRecord r;
    while (read_double_float(fp, &r.time_d, &r.sensordepth))
        records.push_back(r);

    fclose(fp);
    return records;
}

/* --------------------------------------------------------------------------
 * Interpolation helpers
 * -------------------------------------------------------------------------- */

/* Find the index of the last record with time_d <= t using binary search.
 * Returns 0 if all records are after t (clamp to first).
 * Returns records.size()-1 if all records are before t (clamp to last). */
template <typename T>
static size_t find_lower(const std::vector<T>& records, double t) {
    /* Binary search on time_d field */
    size_t lo = 0, hi = records.size() - 1;
    while (lo < hi) {
        size_t mid = lo + (hi - lo + 1) / 2;
        if (records[mid].time_d <= t)
            lo = mid;
        else
            hi = mid - 1;
    }
    return lo;
}

bool mb_interp_attitude(const std::vector<MbBaaRecord>& records,
                        double time_d, double* roll, double* pitch) {
    if (records.empty()) {
        *roll  = std::numeric_limits<double>::quiet_NaN();
        *pitch = std::numeric_limits<double>::quiet_NaN();
        return false;
    }

    /* Clamp to edges */
    if (time_d <= records.front().time_d) {
        *roll  = records.front().roll;
        *pitch = records.front().pitch;
        return true;
    }
    if (time_d >= records.back().time_d) {
        *roll  = records.back().roll;
        *pitch = records.back().pitch;
        return true;
    }

    size_t i = find_lower(records, time_d);
    /* i is the last record with time_d <= target; i+1 is guaranteed to exist */
    double t0 = records[i].time_d;
    double t1 = records[i + 1].time_d;
    double frac = (t1 > t0) ? (time_d - t0) / (t1 - t0) : 0.0;

    *roll  = records[i].roll  + frac * (records[i + 1].roll  - records[i].roll);
    *pitch = records[i].pitch + frac * (records[i + 1].pitch - records[i].pitch);
    return true;
}

bool mb_interp_heading(const std::vector<MbBahRecord>& records,
                       double time_d, double* heading) {
    if (records.empty()) {
        *heading = std::numeric_limits<double>::quiet_NaN();
        return false;
    }

    if (time_d <= records.front().time_d) { *heading = records.front().heading; return true; }
    if (time_d >= records.back().time_d)  { *heading = records.back().heading;  return true; }

    size_t i = find_lower(records, time_d);
    double t0 = records[i].time_d;
    double t1 = records[i + 1].time_d;
    double frac = (t1 > t0) ? (time_d - t0) / (t1 - t0) : 0.0;

    /* Interpolate heading with 0/360 wraparound */
    double h0 = records[i].heading;
    double h1 = records[i + 1].heading;
    double diff = h1 - h0;
    if (diff >  180.0) diff -= 360.0;
    if (diff < -180.0) diff += 360.0;
    double h = h0 + frac * diff;
    if (h < 0.0)   h += 360.0;
    if (h >= 360.0) h -= 360.0;
    *heading = h;
    return true;
}

bool mb_interp_sensordepth(const std::vector<MbBasRecord>& records,
                           double time_d, double* sensordepth) {
    if (records.empty()) {
        *sensordepth = std::numeric_limits<double>::quiet_NaN();
        return false;
    }

    if (time_d <= records.front().time_d) { *sensordepth = records.front().sensordepth; return true; }
    if (time_d >= records.back().time_d)  { *sensordepth = records.back().sensordepth;  return true; }

    size_t i = find_lower(records, time_d);
    double t0 = records[i].time_d;
    double t1 = records[i + 1].time_d;
    double frac = (t1 > t0) ? (time_d - t0) / (t1 - t0) : 0.0;

    *sensordepth = records[i].sensordepth
                 + frac * (records[i + 1].sensordepth - records[i].sensordepth);
    return true;
}
