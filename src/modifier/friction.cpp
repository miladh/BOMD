#include "friction.h"
#include "../molecularSystem/molecularsystem.h"

using namespace bomd;

Friction::Friction(MolecularSystem *molecularSystem) :
    Modifier(molecularSystem)
{
}

void Friction::apply() {

    for(Atom* atom : m_molecularSystem->atoms()) {
        rowvec coreVelocity = atom->coreVelocity() * m_frictionConstant;
        atom->setCoreVelocity(coreVelocity);
    }
}

double Friction::frictionConstant() const
{
    return m_frictionConstant;
}

void Friction::setFrictionConstant(const double &frictionConstant)
{
    m_frictionConstant = frictionConstant;
}


