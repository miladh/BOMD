#ifndef MODIFIER_H
#define MODIFIER_H


namespace bomd{
class MolecularSystem;

class Modifier
{
public:
    Modifier(MolecularSystem* molecularSystem);
    virtual void apply() = 0;

protected:
    MolecularSystem* m_molecularSystem;
};

}
#endif // MODIFIER_H
