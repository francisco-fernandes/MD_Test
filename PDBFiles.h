typedef struct PDBInfo {
	char *fileName;
	unsigned int numAtoms;
	unsigned int numResidues;
	unsigned int numChains;
	unsigned int numModels;
	unsigned int startPositionInAtomsArray;
	unsigned int startPositionInResiduesArray;
} PDBInfo;
