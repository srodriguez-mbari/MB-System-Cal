/*--------------------------------------------------------------------
 *    The MB-system:  mb_ancillary_io.h
 *
 *    Copyright (c) 2026 by
 *    Sebastian Rodriguez (seroma09@gmail.com)
 *      Monterey Bay Aquarium Research Institute
 *
 *    See README.md file for copying and redistribution conditions.
 *--------------------------------------------------------------------
 *
 * Readers for the binary and ASCII ancillary files written by mbpreprocess:
 *
 *   .fnv  ASCII  per-ping integrated nav (19 columns)
 *   .baa  binary async attitude at source sensor rate (double+float+float)
 *   .bah  binary async heading  (double+float)
 *   .bas  binary async sensor depth (double+float)
 *   .bsa  binary sync attitude per ping (double+float+float) — same layout as .baa
 *
 * All binary files are big-endian.
 * All time_d values are Unix epoch seconds.
 *
 * These readers are intentionally kept in src/utilities/ rather than src/mbio/
 * because src/mbio is an all-C library and these use C++ (std::vector).
 */

#pragma once

#include <cmath>
#include <string>
#include <vector>

/* Per-ping record from .fnv ASCII ancillary file.
 * Nav reflects the TARGET SENSOR position (lever-arm corrected),
 * not the raw INS fix. */
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

/* Async attitude record from .baa (or sync from .bsa) — 16 bytes on disk.
 * Format: big-endian double time_d + float roll + float pitch. */
struct MbBaaRecord {
    double time_d;  /* Unix epoch seconds */
    float  roll;    /* degrees */
    float  pitch;   /* degrees */
};

/* Async heading record from .bah — 12 bytes on disk.
 * Format: big-endian double time_d + float heading. */
struct MbBahRecord {
    double time_d;  /* Unix epoch seconds */
    float  heading; /* degrees true */
};

/* Async sensor depth record from .bas — 12 bytes on disk.
 * Format: big-endian double time_d + float sensordepth. */
struct MbBasRecord {
    double time_d;      /* Unix epoch seconds */
    float  sensordepth; /* meters, positive down */
};

/* Read .fnv ASCII ancillary file.
 * Returns records in file order (time-sorted by mbpreprocess).
 * Returns empty vector if file does not exist — not treated as an error. */
std::vector<MbFnvRecord> mb_read_fnv(const std::string& path);

/* Read .baa binary ancillary file (async attitude: time_d, roll, pitch).
 * Also works for .bsa (sync attitude — identical binary layout).
 * Returns empty vector if file does not exist. */
std::vector<MbBaaRecord> mb_read_baa(const std::string& path);

/* Read .bah binary ancillary file (async heading: time_d, heading).
 * Returns empty vector if file does not exist. */
std::vector<MbBahRecord> mb_read_bah(const std::string& path);

/* Read .bas binary ancillary file (async sensor depth: time_d, sensordepth).
 * Returns empty vector if file does not exist. */
std::vector<MbBasRecord> mb_read_bas(const std::string& path);

/* Linearly interpolate roll and pitch from baa records at time_d.
 * Clamps to edge values when time_d is outside the record range.
 * Sets *roll and *pitch to NaN and returns false if records is empty. */
bool mb_interp_attitude(const std::vector<MbBaaRecord>& records,
                        double time_d, double* roll, double* pitch);

/* Linearly interpolate heading from bah records at time_d.
 * Handles 0/360 wraparound. Returns false if records is empty. */
bool mb_interp_heading(const std::vector<MbBahRecord>& records,
                       double time_d, double* heading);

/* Linearly interpolate sensor depth from bas records at time_d.
 * Returns false if records is empty. */
bool mb_interp_sensordepth(const std::vector<MbBasRecord>& records,
                           double time_d, double* sensordepth);
