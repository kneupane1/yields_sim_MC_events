#include "eff_corr.hpp"
#include <iostream>
EffCorr::~EffCorr() {} ///Since your class does not involve dynamic memory allocation (e.g., using new),
// there's no need to define a custom destructor (e.g., ~EffCorr() {}). The compiler-generated destructor
//  will work fine in this case, and you can omit the empty destructor.

////////////////////////////////// PROTON EFF FACTORS

float EffCorr::PROT_EFF_CORR_FACT(float mom_mes, float theta, float phi)
{
    if (theta >= 0 && theta < 20)
        for (int p = 0; p < FD_SEC; p++)
        {
            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                /////////// Map the measured momentum and missing momentum /////////////
                float mom = mom_mes / ((FDProtMap[0][p][0]) * pow(mom_mes, 2) + (FDProtMap[0][p][1]) * mom_mes + (FDProtMap[0][p][2]));
                // std::cout << " Prot mom_mes is : " << mom_mes << "  mom_miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (FDProtCoef[0][p][1]) << std::endl;
                // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
                //           << ((FDProtCoef[0][p][0]) * pow(mom, 2) + (FDProtCoef[0][p][1]) * mom + (FDProtCoef[0][p][2])) << std::endl;

                return (FDProtCoef[0][p][0]) * pow(mom, 2) + (FDProtCoef[0][p][1]) * mom + (FDProtCoef[0][p][2]);
            }
        }

    else if (theta >= 20 && theta < 40)
    {

        for (int p = 0; p < FD_SEC; p++)
        {

            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {

                float mom = mom_mes / ((FDProtMap[1][p][0]) * pow(mom_mes, 2) + (FDProtMap[1][p][1]) * mom_mes + (FDProtMap[1][p][2]));

                // std::cout << " Prot mom_mes is : " << mom_mes << "  mom_miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (FDProtCoef[1][p][0]) << std::endl;
                // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
                //           << ((FDProtCoef[1][p][0]) * pow(mom, 2) + (FDProtCoef[1][p][1]) * mom + (FDProtCoef[1][p][2])) << std::endl;

                return (FDProtCoef[1][p][0]) * pow(mom, 2) + (FDProtCoef[1][p][1]) * mom + (FDProtCoef[1][p][2]);
            }
        }
    }
    else if (theta >= 40 && theta <= 180)
    {
        for (int p = 0; p < CD_SEC; p++)
        {

            if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((CDProtMap[p][0]) * pow(mom_mes, 2) + (CDProtMap[p][1]) * mom_mes + (CDProtMap[p][2]));

                // std::cout << " Prot mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (CDProtCoef[p][0]) << std::endl;
                // std::cout << " cd sec is : " << p + 1 << "  eff corr val =  "
                //           << ((CDProtCoef[p][0]) * pow(mom, 2) + (CDProtCoef[p][1]) * mom + (CDProtCoef[p][2])) << std::endl;

                return (CDProtCoef[p][0]) * pow(mom, 2) + (CDProtCoef[p][1]) * mom + (CDProtCoef[p][2]);
            }
        }
    }

    return 1.0; // to test if you are getting good result or not
}

////////////////////////////////// PIP EFF FACTORS

float EffCorr::PIP_EFF_CORR_FACT(float mom_mes, float theta, float phi)
{
    if (theta >= 0 && theta < 20)
        for (int p = 0; p < FD_SEC; p++)
        {

            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((FDPipMap[0][p][0]) * pow(mom_mes, 2) + (FDPipMap[0][p][1]) * mom_mes + (FDPipMap[0][p][2]));

                // std::cout << " Pip mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (FDPipCoef[0][p][0]) << std::endl;
                // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
                //           << (FDPipCoef[0][p][0]) * pow(mom, 3) + (FDPipCoef[0][p][1]) * pow(mom, 2) + (FDPipCoef[0][p][2]) * pow(mom, 1) + (FDPipCoef[0][p][3]) << std::endl;

                return (FDPipCoef[0][p][0]) * pow(mom, 3) + (FDPipCoef[0][p][1]) * pow(mom, 2) + (FDPipCoef[0][p][2]) * pow(mom, 1) + (FDPipCoef[0][p][3]);
            }
        }

    else if (theta >= 20 && theta < 40)
    {

        for (int p = 0; p < FD_SEC; p++)
        {
            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((FDPipMap[1][p][0]) * pow(mom_mes, 2) + (FDPipMap[1][p][1]) * mom_mes + (FDPipMap[1][p][2]));

                // std::cout << " Pip mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (FDPipCoef[1][p][0]) << std::endl;
                // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
                //           << (FDPipCoef[1][p][0]) * pow(mom, 3) + (FDPipCoef[1][p][1]) * pow(mom, 2) + (FDPipCoef[1][p][2]) * pow(mom, 1) + (FDPipCoef[1][p][3]) << std::endl;

                return (FDPipCoef[1][p][0]) * pow(mom, 3) + (FDPipCoef[1][p][1]) * pow(mom, 2) + (FDPipCoef[1][p][2]) * pow(mom, 1) + (FDPipCoef[1][p][3]);
            }
        }
    }
    else if (theta >= 40 && theta <= 180)
    {
        for (int p = 0; p < CD_SEC; p++)
        {
            if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((CDPipMap[p][0]) * pow(mom_mes, 2) + (CDPipMap[p][1]) * mom_mes + (CDPipMap[p][2]));

                // std::cout << " Pip mom_mes is : " << mom_mes << "  mom miss is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (CDPipCoef[p][0]) << std::endl;
                // std::cout << " cd is : " << p + 1 << "  eff corr val =  "
                //           << ((CDPipCoef[p][0]) * pow(mom, 2) + (CDPipCoef[p][1]) * mom + (CDPipCoef[p][2])) << std::endl;

                return (CDPipCoef[p][0]) * pow(mom, 2) + (CDPipCoef[p][1]) * mom + (CDPipCoef[p][2]);
            }
        }
    }
    return 1.0; // to test if you are getting good result or not
}

////////////////////////////////// PIM EFF FACTORS

float EffCorr::PIM_EFF_CORR_FACT(float mom_mes, float theta, float phi)
{
    if (theta >= 0 && theta < 20)
        for (int p = 0; p < FD_SEC; p++)
        {
            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((FDPimMap[0][p][0]) * pow(mom_mes, 2) + (FDPimMap[0][p][1]) * mom_mes + (FDPimMap[0][p][2]));

                return (FDPimCoef[0][p][0]) * pow(mom, 2) + (FDPimCoef[0][p][1]) * mom + (FDPimCoef[0][p][2]);
            }
        }

    else if (theta >= 20 && theta < 40)
    {

        for (int p = 0; p < FD_SEC; p++)
        {
            if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((FDPimMap[1][p][0]) * pow(mom_mes, 2) + (FDPimMap[1][p][1]) * mom_mes + (FDPimMap[1][p][2]));

                // std::cout << "mom is : " << mom << " theta : " << theta << " phi : " << phi
                //           << "  1st coeff =  " << (FDPimCoef[1][p][0]) << std::endl;
                // std::cout << " fdsec is : " << p + 1 << "  eff corr val =  "
                //           << ((FDPimCoef[1][p][0]) * pow(mom, 2) + (FDPimCoef[1][p][1]) * mom + (FDPimCoef[1][p][2])) << std::endl;

                return (FDPimCoef[1][p][0]) * pow(mom, 2) + (FDPimCoef[1][p][1]) * mom + (FDPimCoef[1][p][2]);
            }
        }
    }
    else if (theta >= 40 && theta <= 180)
    {
        for (int p = 0; p < CD_SEC; p++)
        {
            if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
            {
                float mom = mom_mes / ((CDPimMap[p][0]) * pow(mom_mes, 2) + (CDPimMap[p][1]) * mom_mes + (CDPimMap[p][2]));

                return (CDPimCoef[p][0]) * pow(mom, 2) + (CDPimCoef[p][1]) * mom + (FDPimCoef[1][p][2]);
            }
        }
    }
    return 1.0; // to test if you are getting good result or not
}
float EffCorr::EFF_CORR_FACT(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip, float mom_pim, float theta_pim, float phi_pim)
{
    return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
            (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)) *
            (EffCorr::PIM_EFF_CORR_FACT(mom_pim, theta_pim, phi_pim)));
};
float EffCorr::EFF_CORR_FACT1(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip)
{
    return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
            (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)));
};

////// another approach ///////////////////
////// another approach ///////////////////
////// another approach ///////////////////
////// another approach ///////////////////

// float EffCorr::PROT_EFF_CORR_FACT(float mom, float theta, float phi)
// {
//     if (theta >= 0 && theta < 20)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDL_P; ++i)
//                 {
//                     if (mom > mom_bin_ranges_prot_fdl[i] && mom <= mom_bin_ranges_prot_fdl[i + 1])
//                     {
//                         return eff_corr_val_prot_0_20[i][p];
//                     }
//                 }
//             }
//         }
//     }

//     else if (theta >= 20 && theta < 40)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDH_P; ++i)
//                 {
//                     if (mom > mom_bin_ranges_prot_fdh[i] && mom <= mom_bin_ranges_prot_fdh[i + 1])
//                     {
//                         return eff_corr_val_prot_20_40[i][p];
//                     }
//                 }
//             }
//         }
//     }
//     else if (theta >= 40 && theta <= 180)
//     {

//         for (int p = 0; p < CD_SEC; p++)
//         {
//             if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_CD_P; ++i)
//                 {
//                     if (p != 1)
//                     {
//                         if (mom > mom_bin_ranges_prot_cd[i] && mom <= mom_bin_ranges_prot_cd[i + 1])
//                         {
//                             return eff_corr_val_prot_40_180[i][p];
//                         }
//                     }
//                     else
//                     {
//                         if (mom > mom_bin_ranges_prot_cd1[i] && mom <= mom_bin_ranges_prot_cd1[i + 1])
//                         {
//                             return eff_corr_val_prot_40_180[i][p];
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return 1.0; // to test if you are getting good result or not
// }

////////////////////////////////// PIP EFF FACTORS

// float EffCorr::PIP_EFF_CORR_FACT(float mom, float theta, float phi)
// {
//     if (theta >= 0 && theta < 20)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDL_PIP; ++i)
//                 {
//                     if (mom > mom_bin_ranges_pip_fdl[i] && mom <= mom_bin_ranges_pip_fdl[i + 1])
//                     {
//                         return eff_corr_val_pip_0_20[i][p];
//                     }
//                 }
//             }
//         }
//     }

//     else if (theta >= 20 && theta < 40)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDH_PIP; ++i)
//                 {
//                     if (mom > mom_bin_ranges_pip_fdh[i] && mom <= mom_bin_ranges_pip_fdh[i + 1])
//                     {
//                         return eff_corr_val_pip_20_40[i][p];
//                     }
//                 }
//             }
//         }
//     }

//     else if (theta >= 40 && theta <= 180)
//     {

//         for (int p = 0; p < CD_SEC; p++)
//         {
//             if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_CD_PIP; ++i)
//                 {
//                     if (p != 1)
//                     {
//                         if (mom > mom_bin_ranges_pip_cd[i] && mom <= mom_bin_ranges_pip_cd[i + 1])
//                         {
//                             return eff_corr_val_pip_40_180[i][p];
//                         }
//                     }
//                     else
//                     {
//                         if (mom > mom_bin_ranges_pip_cd1[i] && mom <= mom_bin_ranges_pip_cd1[i + 1])
//                         {
//                             return eff_corr_val_pip_40_180[i][p];
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return 1.0; // to test if you are getting good result or not
// }

// ////////////////////////////////// PIM EFF FACTORS

// float EffCorr::PIM_EFF_CORR_FACT(float mom, float theta, float phi)
// {
//     if (theta >= 0 && theta < 20)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDL_PIM; ++i)
//                 {
//                     if (mom > mom_bin_ranges_pim_fdl[i] && mom <= mom_bin_ranges_pim_fdl[i + 1])
//                     {
//                         return eff_corr_val_pim_0_20[i][p];
//                     }
//                 }
//             }
//         }
//     }

//     else if (theta >= 20 && theta < 40)
//     {
//         for (int p = 0; p < FD_SEC; p++)
//         {
//             if (phi > fd_phi_bin_ranges[p] && phi <= fd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_FDH_PIM; ++i)
//                 {
//                     if (mom > mom_bin_ranges_pim_fdh[i] && mom <= mom_bin_ranges_pim_fdh[i + 1])
//                     {
//                         return eff_corr_val_pim_20_40[i][p];
//                     }
//                 }
//             }
//         }
//     }

//     else if (theta >= 40 && theta < 180)
//     {

//         for (int p = 0; p < CD_SEC; p++)
//         {
//             if (phi > cd_phi_bin_ranges[p] && phi <= cd_phi_bin_ranges[p + 1])
//             {
//                 for (int i = 0; i < MOM_BIN_SIZE_CD_PIM; ++i)
//                 {
//                     if (mom > mom_bin_ranges_pim_cd[i] && mom <= mom_bin_ranges_pim_cd[i + 1])
//                     {
//                         return eff_corr_val_pim_40_180[i][p];
//                     }
//                 }
//             }
//         }
//     }
//     return 1.0; // to test if you are getting good result or not
// }
// float EffCorr::EFF_CORR_FACT(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip, float mom_pim, float theta_pim, float phi_pim)
// {
//     return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
//             (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)) *
//             (EffCorr::PIM_EFF_CORR_FACT(mom_pim, theta_pim, phi_pim)));
// };
// float EffCorr::EFF_CORR_FACT1(float mom_p, float theta_p, float phi_p, float mom_pip, float theta_pip, float phi_pip)
// {
//     return ((EffCorr::PROT_EFF_CORR_FACT(mom_p, theta_p, phi_p)) *
//             (EffCorr::PIP_EFF_CORR_FACT(mom_pip, theta_pip, phi_pip)));
// };
