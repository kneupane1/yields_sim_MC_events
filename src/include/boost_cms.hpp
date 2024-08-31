#ifndef BOOST_CMS_HPP
#define BOOST_CMS_HPP

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRotation.h>
#include <cmath>
#include <vector>

namespace boost_cms
{

    // Function to perform rotation and boosting to the center-of-mass system
    TLorentzVector boostToCMS(const TLorentzVector &particle, const TLorentzVector &gamma, const TLorentzVector &e_prime, float Q2);

    // Function to calculate invariant mass from a vector of TLorentzVectors
    double calculateInvariantMass(const TLorentzVector &particle1, const TLorentzVector &particle2);

    // Function to calculate the theta angle of a particle
    double calculateTheta(const TLorentzVector &particle);

    // Function to calculate the phi angle of a particle
    double calculatePhi(const TLorentzVector &particle);

    // Function to calculate the alpha angle between two particles
    // double calculateAlpha(const TLorentzVector &particle1, const TLorentzVector &particle2, const TVector3 &V3_anti_z);
    double calculateAlpha(const TVector3 &vect1, const TVector3 &vect2, const TVector3 &boostedVect);

} // namespace boost_cms

#endif
