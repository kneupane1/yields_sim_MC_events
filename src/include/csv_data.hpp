#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
        short electron_sector;
        float w;
        float q2;
        float pip_theta;
        float pip_phi;
        float mm2;

        // Static functions can be called without making a new struct
        static std::string header() {
                // Make a string for the header of the csv file
                // return "electron_sector,w,q2,pip_theta,pip_phi,mm2";
                return "pid,mc_pid,mc_npart,mc_event,mc_helicity,mc_weight";

        }

        friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
                os << std::setprecision(10);
                os << data.electron_sector << ",";
                os << data.w << ",";
                os << data.q2 << ",";
                os << data.pip_theta << ",";
                os << data.pip_phi << ",";
                os << data.mm2;

                return os;
        }
};

#endif
