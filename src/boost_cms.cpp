#include "boost_cms.hpp"
#include "constants.hpp"
namespace boost_cms
{

    TLorentzVector boostToCMS(const TLorentzVector &particle, const TLorentzVector &gamma, const TLorentzVector &e_prime, float Q2)
    {
        // Create a copy of the particle to apply transformations
        TLorentzVector boosted_particle = particle;
        TLorentzVector boosted_gamma = gamma;

        // Determine the unit vectors for rotation
        TVector3 uz = gamma.Vect().Unit();                         // uit vector along virtual photon
        TVector3 ux = (gamma.Vect().Cross(e_prime.Vect())).Unit(); // unit vector along e cross e'
        ux.Rotate(3. * PI / 2, uz);                                // rotating ux by 3pi/2 with uz as axis of roration

        // Apply rotation
        TRotation rot;
        rot.SetZAxis(uz, ux).Invert(); // setting TRotation rot

        boosted_gamma.Transform(rot);

        boosted_particle.Transform(rot);

        // Calculate beta for the boost
        // float beta = (sqrt(gamma.E() * gamma.E() + Q2)) / (gamma.E() + MASS_P);
        float beta = (sqrt(boosted_gamma.E() * boosted_gamma.E() + Q2)) / (boosted_gamma.E() + MASS_P);

        // Apply the boost
        boosted_particle.Boost(0, 0, -beta);

        return boosted_particle;
    }

    double calculateInvariantMass(const TLorentzVector &particle1, const TLorentzVector &particle2)
    {
        TLorentzVector resonance;

        resonance += particle1;
        resonance += particle2;

        return resonance.M();
    }

    double calculateTheta(const TLorentzVector &particle)
    {
        return particle.Theta() * (180.0 / PI);
    }

    double calculatePhi(const TLorentzVector &particle)
    {
        double phi = particle.Phi() * (180.0 / PI);
        if (phi < 0)
        {
            phi += 360.0;
        }
        return phi;
    }

    double calculateAlpha(const TVector3 &vect1, const TVector3 &vect2, const TVector3 &boostedVect)
    {

        TVector3 V3_anti_z(0, 0, -1);

        float a_gamma = sqrt(1.0 / (1 - pow(vect1 * V3_anti_z, 2)));
        float b_gamma = -(vect1 * V3_anti_z) * a_gamma;
        TVector3 Vect3_gamma = a_gamma * V3_anti_z + b_gamma * vect1;

        float a_beta = sqrt(1.0 / (1 - pow(vect1 * vect2, 2)));
        float b_beta = -(vect1 * vect2) * a_beta;
        TVector3 Vect3_beta = a_beta * vect2 + b_beta * vect1;

        double alpha = (180.0 / PI) * acos(Vect3_gamma * Vect3_beta);
        if (Vect3_gamma.Cross(Vect3_beta) * boostedVect < 0)
        {
            alpha = 360.0 - alpha;
        }

        return alpha;
    }
} // namespace boost_cms
