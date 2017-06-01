typedef struct Atom {
	unsigned int fromResidueId;
	unsigned char atomIdInsideResidue;
	char element;
	float coordX;
	float coordY;
	float coordZ;
} Atom;

typedef struct Residue {
	unsigned char moleculeId;
	unsigned char fromModelId;
	char fromChain;
	// TODO: remove this?
	signed char numExtraAtoms;
	unsigned int startPositionInAtomsArray;
} Residue;

// NOTE: position 0 in both arrays bellow is not used
Atom *allAtomsArray;
Residue *allResiduesArray;
unsigned int totalNumAtoms;
unsigned int totalNumResidues;

int GetAtomElementCode(Atom *atom);
