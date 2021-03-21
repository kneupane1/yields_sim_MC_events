#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
        short electron_sector;
        float w;
        float q2;
        float pim_mom_mPim;
        float pim_theta_mPim;
        float pim_phi_mPim;
        float mm2_mPim;
        float weight_mPim;
        // float pim_mom_exclusive;
        // float pim_theta_exclusive;
        // float pim_phi_exclusive;
        // float mm2_exclusive;
        // float mm2_exclusive_at_zero;
        // float weight_exclusive;
        // Static functions can be called without making a new struct
        static std::string header() {
                // Make a string for the header of the csv file
                return "electron_sector,w,q2,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight";
                //return "w,q2,pim_mom_exclusive,pim_theta_exclusive,pim_phi_exclusive,mm2_exclusive,mm2_exclusive_at_zero,weight";

        }

        friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
                os << std::setprecision(10);
                os << data.electron_sector << ",";
                os << data.w << ",";
                os << data.q2 << ",";
                os << data.pim_mom_mPim << ",";
                os << data.pim_theta_mPim << ",";
                os << data.pim_phi_mPim << ",";
                os << data.mm2_mPim<<",";
                os << data.weight_mPim<<",";
                // os << data.w << ",";
                // os << data.q2 << ",";
                // os << data.pim_mom_exclusive << ",";
                // os << data.pim_theta_exclusive << ",";
                // os << data.pim_phi_exclusive << ",";
                // os << data.mm2_exclusive << ",";
                // os << data.mm2_exclusive_at_zero<<",";
                // os << data.weight_exclusive<<",";
                return os;
        }
};

#endif
