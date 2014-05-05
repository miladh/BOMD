#include "velocityrescaling.h"
#include "../molecularSystem/molecularsystem.h"

using namespace bomd;

VelocityRescaling::VelocityRescaling(MolecularSystem *molecularSystem) :
    Modifier(molecularSystem)
{
}

void VelocityRescaling::apply() {

    for(Atom* atom : m_molecularSystem->atoms()) {
        rowvec coreVelocity = atom->coreVelocity() * m_rescalingFactor;
        atom->setCoreVelocity(coreVelocity);
    }
}

double VelocityRescaling::rescalingFactor() const
{
    return m_rescalingFactor;
}

void VelocityRescaling::setRescalingFactor(const double &rescalingFactor)
{
    m_rescalingFactor = rescalingFactor;
}


