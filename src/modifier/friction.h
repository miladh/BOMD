#ifndef FRICTION_H
#define FRICTION_H

#include "modifier.h"

namespace bomd {

class Friction : public Modifier
{
public:
    Friction(MolecularSystem *molecularSystem);

    void apply();

    double frictionConstant() const;
    void setFrictionConstant(const double& frictionConstant);

private:
    double m_frictionConstant;
};

}
#endif // FRICTION_H
