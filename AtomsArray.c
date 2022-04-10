#include <stdio.h>
#include <stdlib.h>
#include "AtomsArray.h"
#include "ChemicalEntities.h"

int GetAtomElementCode(Atom *atom) {
	unsigned char molId;
	char *atomName;
	molId = (allResiduesArray[(atom->fromResidueId)].moleculeId);
	atomName = (char *)MoleculesList[molId].atomsCodes[(atom->atomIdInsideResidue)];
	return (int)(ElementIdLookupTable[CHAR2ID(atomName[0])]);
}
