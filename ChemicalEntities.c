#include <stdio.h>
#include <stdlib.h>
#include "ChemicalEntities.h"

int GetElementCode(const char *atomName) {
	return (int)(ElementIdLookupTable[CHAR2ID(atomName[0])]);
}

int GetMoleculeCode(char *moleculeName) {
	static const int charPos[3] = { 0, 2, 1 };
	int i = 0;
	int code = 0;
	do {
		code = (-code);
		code = (int)(MoleculeIdLookupTable[code][CHAR2ID(moleculeName[charPos[i]])]);
		i++;
	} while (code < 0);
	return code;
}
