#ifndef VELOCITYRESCALING_H
#define VELOCITYRESCALING_H

#include "modifier.h"

namespace bomd {

class VelocityRescaling : public Modifier
{
public:
    VelocityRescaling(MolecularSystem *molecularSystem);

    void apply();

    double rescalingFactor() const;
    void setRescalingFactor(const double& rescalingFactor);

private:
    double m_rescalingFactor;
};

}
#endif // VELOCITYRESCALING_H
