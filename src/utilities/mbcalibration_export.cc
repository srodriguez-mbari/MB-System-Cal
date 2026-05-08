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
 * data to a NetCDF4 file for use by the automated boresight and lever-arm
 * calibration pipeline.
 *
 * Per-ping output:
 *   - Nav: lon, lat, sensor depth, heading, speed (from .fnv)
 *   - Raw attitude: roll, pitch, heave (from .baa, before boresight)
 *   - Per-beam: world-frame bath position, beam flags, amplitude
 *
 * File-level metadata includes the platform file path, active sensor
 * indices, kluge flags, and the MB-System version used.
 *
 * Author:   Sebastian Rodriguez
 * Date:     2026-05-08
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>
#include <unistd.h>
#include <string>

#include "mb_define.h"
#include "mb_format.h"
#include "mb_status.h"
#include "mb_ancillary_io.h"

constexpr char program_name[] = "mbcalibration_export";

constexpr char help_message[] =
    "mbcalibration_export exports per-ping calibration data from processed\n"
    "MB-System .mb89 files to NetCDF4 for automated boresight and lever-arm\n"
    "calibration. Ancillary files (.fnv, .baa) must exist alongside the .mb89\n"
    "files (produced by mbpreprocess --output-sensor-fnv).";

constexpr char usage_message[] =
    "mbcalibration_export -I datalist -O output.nc\n"
    "                     [--platform-file=path]\n"
    "                     [--include-raw]\n"
    "                     [-V] [-H]";

static struct option long_options[] = {
    {"platform-file",  required_argument, nullptr, 0},
    {"include-raw",    no_argument,       nullptr, 0},
    {"help",           no_argument,       nullptr, 'H'},
    {"verbose",        no_argument,       nullptr, 'V'},
    {nullptr,          0,                 nullptr, 0}
};

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
    char read_file[MB_PATH_MAXLINE] = "";
    char output_file[MB_PATH_MAXLINE] = "";
    char platform_file[MB_PATH_MAXLINE] = "";
    bool include_raw = false;

    /* Parse arguments */
    {
        bool help = false;
        int c;
        int option_index = 0;
        while ((c = getopt_long(argc, argv, "HhI:i:O:o:Vv",
                                long_options, &option_index)) != -1) {
            switch (c) {
            case 0:
                if (strcmp(long_options[option_index].name, "platform-file") == 0)
                    snprintf(platform_file, sizeof(platform_file), "%s", optarg);
                else if (strcmp(long_options[option_index].name, "include-raw") == 0)
                    include_raw = true;
                break;
            case 'H':
            case 'h':
                help = true;
                break;
            case 'I':
            case 'i':
                snprintf(read_file, sizeof(read_file), "%s", optarg);
                break;
            case 'O':
            case 'o':
                snprintf(output_file, sizeof(output_file), "%s", optarg);
                break;
            case 'V':
            case 'v':
                verbose++;
                break;
            default:
                fprintf(stderr, "\nOption -%c not recognized\n", c);
                fprintf(stderr, "\nUsage: %s\n", usage_message);
                exit(MB_ERROR_BAD_PARAMETER);
            }
        }

        if (help) {
            fprintf(stderr, "\nProgram %s\nMB-system Version %s\n\n", program_name, MB_VERSION);
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

    if (verbose > 0) {
        fprintf(stderr, "\nProgram %s\nMB-system Version %s\n", program_name, MB_VERSION);
        fprintf(stderr, "\nInput:           %s\n", read_file);
        fprintf(stderr, "Output:          %s\n", output_file);
        fprintf(stderr, "Platform file:   %s\n", strlen(platform_file) > 0 ? platform_file : "(none)");
        fprintf(stderr, "Include raw:     %s\n\n", include_raw ? "yes" : "no");
    }

    /* Determine whether input is a datalist or a single file */
    if (format == 0)
        mb_get_format(verbose, read_file, nullptr, &format, &error);
    const bool read_datalist = format < 0;

    /* Open input */
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
            fprintf(stderr, "\nProgram <%s> Terminated\n", program_name);
            exit(MB_ERROR_OPEN_FAIL);
        }
        read_data = mb_datalist_read(verbose, datalist, file, dfile,
                                     &format, &file_weight, &error) == MB_SUCCESS;
    } else {
        snprintf(file, sizeof(file), "%s", read_file);
        read_data = true;
    }

    /* Per-file beam arrays (allocated by MB-System's memory manager) */
    char   *beamflag = nullptr;
    double *bath = nullptr;
    double *bathacrosstrack = nullptr;
    double *bathalongtrack = nullptr;
    double *amp = nullptr;
    double *ss = nullptr;
    double *ssacrosstrack = nullptr;
    double *ssalongtrack = nullptr;
    /* bathlon/bathlat added in Task 6 (mb_read, not mb_get_all, populates these) */

    /* Totals across all files */
    int n_files = 0;
    int n_pings_total = 0;
    int n_beams_total = 0;

    /* ----------------------------------------------------------------
     * Main file loop
     * ---------------------------------------------------------------- */
    while (read_data) {
        n_files++;

        /* Build ancillary file paths from the .mb89 path */
        std::string fnv_path = std::string(file) + ".fnv";
        std::string baa_path = std::string(file) + ".baa";

        /* Load ancillary files into memory */
        auto fnv_recs = mb_read_fnv(fnv_path);
        auto baa_recs = mb_read_baa(baa_path);

        if (fnv_recs.empty()) {
            fprintf(stderr, "Warning: no .fnv records for %s — skipping\n", file);
            goto next_file;
        }
        if (baa_recs.empty()) {
            fprintf(stderr, "Warning: no .baa records for %s — skipping\n", file);
            goto next_file;
        }

        if (verbose > 0) {
            fprintf(stderr, "File %d: %s\n", n_files, file);
            fprintf(stderr, "  .fnv pings: %zu   .baa records: %zu\n",
                    fnv_recs.size(), baa_recs.size());
        }

        {
            /* Open .mb89 for reading */
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
                fprintf(stderr, "Warning: cannot open %s: %s — skipping\n",
                        file, message);
                goto next_file;
            }

            /* Register beam arrays with MB-System's memory manager.
             * bathlon/bathlat deferred to Task 6 (require mb_read not mb_get_all). */
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(char),   (void **)&beamflag, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bath, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathacrosstrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY,
                              sizeof(double), (void **)&bathalongtrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE,
                              sizeof(double), (void **)&amp, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ss, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ssacrosstrack, &error);
            mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN,
                              sizeof(double), (void **)&ssalongtrack, &error);

            if (error != MB_ERROR_NO_ERROR) {
                fprintf(stderr, "Warning: memory allocation failed for %s — skipping\n", file);
                mb_close(verbose, &mbio_ptr, &error);
                goto next_file;
            }

            /* ---- Ping loop ---- */
            int n_pings_file = 0;
            void *store_ptr;
            char comment[MB_COMMENT_MAXLINE];

            while (error <= MB_ERROR_NO_ERROR) {
                int kind;
                int time_i[7];
                double time_d, navlon, navlat, speed, heading;
                double distance, altitude, sensordepth;

                int status = mb_get_all(verbose, mbio_ptr, &store_ptr,
                                        &kind, time_i, &time_d,
                                        &navlon, &navlat, &speed, &heading,
                                        &distance, &altitude, &sensordepth,
                                        &beams_bath, &beams_amp, &pixels_ss,
                                        beamflag, bath, amp,
                                        bathacrosstrack, bathalongtrack,
                                        ss, ssacrosstrack, ssalongtrack,
                                        comment, &error);

                /* Time gaps are not errors here */
                if (error == MB_ERROR_TIME_GAP) {
                    error = MB_ERROR_NO_ERROR;
                    status = MB_SUCCESS;
                }

                if (status != MB_SUCCESS || kind != MB_DATA_DATA)
                    continue;

                /* Interpolate raw attitude from .baa at this ping's time */
                double raw_roll, raw_pitch;
                mb_interp_attitude(baa_recs, time_d, &raw_roll, &raw_pitch);

                /* TODO Task 5: write ping + beams to NetCDF4 output */
                (void)raw_roll;
                (void)raw_pitch;

                n_pings_file++;

                /* Count good beams */
                for (int ib = 0; ib < beams_bath; ib++) {
                    if (mb_beam_ok(beamflag[ib]))
                        n_beams_total++;
                }
            }

            n_pings_total += n_pings_file;

            if (verbose > 0)
                fprintf(stderr, "  Survey pings read: %d\n", n_pings_file);

            mb_close(verbose, &mbio_ptr, &error);
            error = MB_ERROR_NO_ERROR;
        }

next_file:
        /* Advance datalist */
        if (read_datalist)
            read_data = mb_datalist_read(verbose, datalist, file, dfile,
                                         &format, &file_weight, &error) == MB_SUCCESS;
        else
            read_data = false;
    }

    if (read_datalist)
        mb_datalist_close(verbose, &datalist, &error);

    /* Summary */
    fprintf(stderr, "\nProgram %s completed:\n", program_name);
    fprintf(stderr, "  Files processed:    %d\n", n_files);
    fprintf(stderr, "  Survey pings:       %d\n", n_pings_total);
    fprintf(stderr, "  Good beam returns:  %d\n", n_beams_total);
    fprintf(stderr, "  Output file:        %s (not yet written — Task 5 pending)\n",
            output_file);

    return MB_ERROR_NO_ERROR;
}
