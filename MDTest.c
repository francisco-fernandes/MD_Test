//#ifdef _MSC_VER
//#define _CRT_SECURE_NO_WARNINGS
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "PDBFiles.h"
#include "Statistics.h"
#include "AtomsArray.h"
#include "SVGGraphics.h"
// TODO: put definitions inside "ChemicalEntities.c"
// TODO: create "PDBFormat.h/c"
#include "ChemicalEntities.h"

#define VERSION "0.4"
#define PAUSE_AT_EXIT 1


PDBInfo *allPDBsArray = NULL;
int numPDBs = 0;

// NOTE: position 0 in both arrays bellow is not used
//Atom *allAtomsArray = NULL;
//Residue *allResiduesArray = NULL;
//unsigned int totalNumAtoms = 0;
//unsigned int totalNumResidues = 0;
unsigned int maxAtomsArraySize = 0;
unsigned int maxResiduesArraySize = 0;

int EqualString(char *str1, const char *str2){
	while ((*str1) != '\0' && (*str1) == (*str2)){
		str1++;
		str2++;
	}
	if ((*str1) == (*str2)) return 1;
	return 0;
}

void ExitErrorMessage(char *msg, int newline){
	if (newline) putchar('\n');
	fprintf(stderr,"> ERROR: %s\n", msg);
	#if PAUSE_AT_EXIT == 1
	getchar();
	#endif
	exit(-1);
}



// TODO: ignore missing header if not found
// TODO: create PDBFile struct with array of Model structs with array of Chain structs
// TODO: only load the first model
// TODO: update total number of warnings to show at the end
// TODO: show warning when unsupported HETATM atoms are found
// TODO: allow atoms' names of residues to be given non-sequentially, and perform extra sanity checks (missing atoms, repeated atoms)
//		 by first implementing GetAtomIdInsideResidueFromAtomName()
// TODO: when atoms are missing (e.g hydrogens in X-ray) fill data with position based on angles, etc.
int LoadPDBFile(int pdbFileId, int onlyFirstModel, int showWarnings){
	FILE *pdbFile;
	char c, prevc, pdbID[5], recordType[7];
	int i, numItemsRead;
	int numModels, numChains, *allChainsNumResidues, numChainResidues, numModelAtoms, numPdbAtoms, numPdbResidues;
	char *allChainsCharIds;
	int atomId, residueId, prevAtomId, prevResidueId, chainId;
	char atomName[5], residueName[4], chainCharId, elementSymbol[3], prevChainCharId;
	float coordX, coordY, coordZ;
	int residueCode, numResidueAtoms, atomCountInsideResidue;
	const char (*atomsNamesInsideResidue)[MAX_NUM_ATOMS_PER_AA][5];
	int validAtom;
	int numWarnings;
	//float occupancy, tempFactor;
	printf("> Loading structures from file <%s> ...\n", allPDBsArray[pdbFileId].fileName);
	fflush(stdout);
	if ((pdbFile = fopen(allPDBsArray[pdbFileId].fileName, "r")) == NULL){
		fprintf(stderr,"\n> ERROR: Cannot open PDB file\n");
		return 0;
	}
	for (i = 0; i < 6; i++){
		c = fgetc(pdbFile);
		recordType[i] = c;
	}
	recordType[6] = '\0';
	if (!EqualString(recordType, "HEADER")){
		fprintf(stderr,"\n> ERROR: Invalid PDB file\n");
		return 0;
	}
	for (i = 6; i != 62; i++){
		c = fgetc(pdbFile);
		if (c == EOF || c == '\n'){
			fprintf(stderr,"\n> ERROR: PDB id not found\n");
			return 0;
		}
	}
	for (i = 0; i < 4; i++){
		c = fgetc(pdbFile);
		pdbID[i] = c;
	}
	pdbID[4] = '\0';
	printf("\t[%s: \"", pdbID);
	while (c != EOF && c != '\n') c = fgetc(pdbFile);
	for (i = 0; i < 6; i++){
		c = fgetc(pdbFile);
		recordType[i] = c;
	}
	recordType[6] = '\0';
	if (!EqualString(recordType, "TITLE ")){
		fprintf(stderr,"\n> ERROR: PDB title not found\n");
		return 0;
	}
	while (c == ' ') c = fgetc(pdbFile);
	prevc = ' ';
	while (c != EOF && c != '\n'){
		if (prevc != ' ' && c >= 65 && c <= 90) c = (char)(c + 32);
		prevc = c;
		c = fgetc(pdbFile);
		if (c == ' ' && prevc == ' ') break;
		putchar(prevc);
	}
	printf("\" (");
	do {
		while (c != EOF && c != '\n') c = fgetc(pdbFile);
		numItemsRead = fscanf(pdbFile, " %6[^ ] ", recordType);
		c = fgetc(pdbFile);
	} while (numItemsRead == 1 && !EqualString(recordType, "EXPDTA"));
	if (numItemsRead != 1){
		fprintf(stderr,"\n> ERROR: PDB experimental technique not found\n");
		return 0;
	}
	while (c == ' ') c = fgetc(pdbFile);
	while (c != EOF && c != '\n'){
		prevc = c;
		c = fgetc(pdbFile);
		if (c == ' ' && prevc == ' ') break;
		putchar(prevc);
	}
	printf(")] ... \n");
	fflush(stdout);
	while (c != EOF){ // go to first ATOM record
		while (c != EOF && c != '\n') c = fgetc(pdbFile);
		c = fgetc(pdbFile);
		if (c != 'A') continue;
		prevc = fgetc(pdbFile);
		if (prevc != 'T') continue;
		ungetc(prevc, pdbFile);
		ungetc(c, pdbFile);
		break;
	}
	if (c == EOF){
		fprintf(stderr,"> ERROR: No ATOM records found\n");
		return 0;
	}
	allPDBsArray[pdbFileId].startPositionInAtomsArray = (totalNumAtoms + 1);
	allPDBsArray[pdbFileId].startPositionInResiduesArray = (totalNumResidues + 1);
	numModels = 0;
	numChains = 0;
	chainId = (-1);
	allChainsCharIds = NULL;
	allChainsNumResidues = NULL;
	numModelAtoms = 0;
	numPdbAtoms = 0;
	numPdbResidues = 0;
	numChainResidues = 0;
	prevAtomId = INT_MAX;
	prevResidueId = 0;
	prevChainCharId = '\0';
	residueCode = 0;
	numResidueAtoms = 0;
	atomCountInsideResidue = 0;
	atomsNamesInsideResidue = NULL;
	validAtom = 0;
	numWarnings = 0;
	while (1){ // process each line
		for (i = 0; i < 6; i++){
			c = fgetc(pdbFile);
			recordType[i] = c;
		}
		recordType[6] = '\0';
		if (c == EOF){ // end of file, save last chain and model stats
			// TODO: this code is duplicated bellow inside the main loop; find way to re-use it
			if ((prevResidueId != 0) && (atomCountInsideResidue != numResidueAtoms)){ // check last residue in file
				if (atomCountInsideResidue < numResidueAtoms){ // less atoms
					if (numModels == 1) fprintf(stderr, "> WARNING: Residue #%d (\"%s\") has less atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
					for (; atomCountInsideResidue != numResidueAtoms; atomCountInsideResidue++){ // fill remaining invalid atom entries for this residue
						totalNumAtoms++;
						allAtomsArray[totalNumAtoms].fromResidueId = 0;
						allAtomsArray[totalNumAtoms].atomIdInsideResidue = 0;
						allAtomsArray[totalNumAtoms].element = '\0';
						allAtomsArray[totalNumAtoms].coordX = 0.0;
						allAtomsArray[totalNumAtoms].coordY = 0.0;
						allAtomsArray[totalNumAtoms].coordZ = 0.0;
					}
				}
				else { // more atoms
					if (numModels == 1) fprintf(stderr, "> WARNING: Residue #%d (\"%s\") has more atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
				}
			}
			if (chainId != (-1) && allChainsNumResidues[chainId] == 0) allChainsNumResidues[chainId] = numChainResidues;
			numPdbAtoms += numModelAtoms;
			break;
		} // end of file check
		c = recordType[0]; // check type of each line
		if ((c == 'A' || c == 'H') && (recordType[1]!='N')){ // "ATOM", "HETATM" (but not "ANISOU")
			// TODO: load residueName as fixed length string to work with one-letter DNA codes
			numItemsRead = fscanf(pdbFile, " %d %4s %3s %c %d %f %f %f %*f %*f %2s ", &atomId, atomName, residueName, &chainCharId, &residueId, &coordX, &coordY, &coordZ, elementSymbol);
			if (numItemsRead != 9){
				fprintf(stderr,"> ERROR: Invalid \"%s\" record\n", recordType);
				return 0;
			}
			if (atomId < prevAtomId){ // new model
				numPdbAtoms += numModelAtoms;
				numModels++;
				numModelAtoms = 0;
				prevChainCharId = '\0';
			}
			if (chainCharId != prevChainCharId){ // different chain
				if (chainId != (-1) && allChainsNumResidues[chainId] == 0) allChainsNumResidues[chainId] = numChainResidues; // save previous chain stats
				for (chainId = 0; chainId < numChains; chainId++){ // get current chain id
					if (allChainsCharIds[chainId] == chainCharId) break;
				}
				if (chainId == numChains){ // new chain
					numChains++;
					allChainsCharIds = (char *)realloc(allChainsCharIds, numChains*sizeof(char));
					allChainsNumResidues = (int *)realloc(allChainsNumResidues, numChains*sizeof(int));
					allChainsCharIds[chainId] = chainCharId;
					allChainsNumResidues[chainId] = 0;
				}
				numChainResidues = 0;
				//prevResidueId = 0; // do not enable this because it is needed next to check the number of atoms of the previous chain
			}
			if (residueId != prevResidueId){ // new residue
				if (prevResidueId != 0){ // check previous residue (if this is not the 1st atom in the file or model)
					if (atomCountInsideResidue != numResidueAtoms){ // different number of atoms
						if (atomCountInsideResidue < numResidueAtoms){ // less atoms
							if (atomCountInsideResidue == (numResidueAtoms - 1) && (residueCode == 5 || residueCode == 9 || residueCode == 12)) { // special cases for "HG" from Cysteine, "HE2" from Histidine and for "HZ3" from Lysine (last atom in all)
								if (numModels == 1) fprintf(stderr, "> INFO: [Model:%d,Chain:%c] Hydrogen atom #%d (\"%s\") is absent from residue #%d (\"%s\")\n", numModels, chainCharId, prevAtomId, (*atomsNamesInsideResidue)[atomCountInsideResidue], prevResidueId, (MoleculesList[residueCode].abbrev));
							}
							else if (numModels == 1) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has less atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
							// TODO: when atom array filling is not sequencial, check missing atoms starting at the position of the 1st atom of the residue
							for (; atomCountInsideResidue != numResidueAtoms; atomCountInsideResidue++){ // fill remaining invalid atom entries for this residue
								totalNumAtoms++;
								allAtomsArray[totalNumAtoms].fromResidueId = prevResidueId;
								// TODO: when the atom is in fact invalid, this should be an invalid index (-1?) instead of pointing to the ordered id
								allAtomsArray[totalNumAtoms].atomIdInsideResidue = atomCountInsideResidue;
								allAtomsArray[totalNumAtoms].element = '\0';
								allAtomsArray[totalNumAtoms].coordX = 0.0;
								allAtomsArray[totalNumAtoms].coordY = 0.0;
								allAtomsArray[totalNumAtoms].coordZ = 0.0;
								numModelAtoms++;
							}
						}
						else { // more atoms
							// TODO: remove from here (and 'if' clause above) since this is already checked bellow, and no more atoms are added if the limit is reached
							if (numModels == 1) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has more atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
						}
					}
				}
				residueCode = GetMoleculeCode(residueName);
				if (residueCode == 0){
					// TODO: do not fail, simply discard molecule
					fprintf(stderr,"> ERROR: Unknown residue #%d (\"%s\")\n", residueId, residueName);
					return 0;
				}
				atomsNamesInsideResidue = &(MoleculesList[residueCode].atomsCodes);
				numResidueAtoms = (MoleculesList[residueCode].numberOfAtoms);
				atomCountInsideResidue = 0;
				numChainResidues++;
				numPdbResidues++;
				validAtom = 1;
				totalNumResidues++;
				if (totalNumResidues == maxResiduesArraySize){
					maxResiduesArraySize += ARRAY_GROW_SIZE;
					allResiduesArray = (Residue *)realloc(allResiduesArray, maxResiduesArraySize*sizeof(Residue));
				}
				allResiduesArray[totalNumResidues].fromChain = chainCharId;
				allResiduesArray[totalNumResidues].fromModelId = numModels;
				allResiduesArray[totalNumResidues].moleculeId = residueCode;
				allResiduesArray[totalNumResidues].numExtraAtoms = 0;
				allResiduesArray[totalNumResidues].startPositionInAtomsArray = (totalNumAtoms + 1); // will start in the next position
				if ((totalNumAtoms+numResidueAtoms) >= maxAtomsArraySize){ // alloc space for all atoms of this residue too if needed
					maxAtomsArraySize += ARRAY_GROW_SIZE;
					allAtomsArray = (Atom *)realloc(allAtomsArray, maxAtomsArraySize*sizeof(Atom));
				}
			} else { // another atom from same residue
				if (atomCountInsideResidue < numResidueAtoms) validAtom = 1; // check if there are more atoms than expected
				else {
					if (numModels == 1) fprintf(stderr, "> WARNING: [Model:%d,Chain:%c] Residue #%d (\"%s\") has more atoms than expected (%d)\n", numModels, chainCharId, prevResidueId, (MoleculesList[residueCode].abbrev), numResidueAtoms);
					validAtom = 0;
				}
			}
			// TODO: restructure this valid atom check, e.g.: if(validAtom) {...} else {...}
			if (validAtom && !EqualString(atomName, (*atomsNamesInsideResidue)[atomCountInsideResidue])){ // invalid atom name
				if (residueCode == 9 && atomCountInsideResidue == 14 && EqualString(atomName, "HD2")) { // special case for "HD1" from Histidine (14th atom)
					if (numModels == 1) fprintf(stderr, "> INFO: [Model:%d,Chain:%c] Hydrogen atom #%d (\"%s\") is absent from residue #%d (\"%s\")\n", numModels, chainCharId, atomId, (*atomsNamesInsideResidue)[atomCountInsideResidue], residueId, residueName);
					totalNumAtoms++;
					allAtomsArray[totalNumAtoms].fromResidueId = totalNumResidues; // fill an empty "fake" atom for "HD1" so "HD2" will be filled up correctly next
					allAtomsArray[totalNumAtoms].atomIdInsideResidue = atomCountInsideResidue;
					allAtomsArray[totalNumAtoms].element = '\0';
					allAtomsArray[totalNumAtoms].coordX = 0.0;
					allAtomsArray[totalNumAtoms].coordY = 0.0;
					allAtomsArray[totalNumAtoms].coordZ = 0.0;
					atomCountInsideResidue++;
					numModelAtoms++;
				} else if (numChainResidues == 1){ // check if first residue contains extra initial hydrogens
					if (EqualString(atomName, "H1")) { // atom "H" from 1st residue might be named "H1" instead of "HA" (only the name is changed)
						if (numModels == 1) fprintf(stderr, "> INFO: [Model:%d,Chain:%c] Atom #%d from residue #%d (\"%s\") is named \"%s\"\n", numModels, chainCharId, atomId, residueId, residueName, atomName);
					}
					else if (EqualString(atomName, "H2") || EqualString(atomName, "H3")) { // initial hydrogens "H2" and "H3" might be present
						if (numModels == 1) fprintf(stderr, "> INFO: [Model:%d,Chain:%c] Residue #%d (\"%s\") has a terminal atom #%d \"%s\"\n", numModels, chainCharId, residueId, residueName, atomId, atomName);
						// TODO: load this extra atom too
						validAtom = 0;
					}
					
				}
				else if (EqualString(atomName, "OXT") || EqualString(atomName, "HXT")){ // terminal oxygen or hydrogen might be present
					if (numModels == 1) fprintf(stderr,"> INFO: [Model:%d,Chain:%c] Residue #%d (\"%s\") has a terminal atom #%d \"%s\"\n", numModels, chainCharId, residueId, residueName, atomId, atomName);
					// TODO: if an extra 1st or last atom exists, it will mess the hardcoded atoms bonds ids
					validAtom = 0;
					/*
					atomCountInsideResidue--;
					numModelAtoms--;
					*/
					/*
					allResiduesArray[totalNumResidues].numExtraAtoms++;
					numResidueAtoms++;
					totalNumAtoms++;
					if (totalNumAtoms >= maxAtomsArraySize){
						maxAtomsArraySize += ARRAY_GROW_SIZE;
						allAtomsArray = (Atom *)realloc(allAtomsArray, maxAtomsArraySize*sizeof(Atom));
					}
					*/
				}
				else {
					// TODO: instead of failing, search entire list of atoms' names for correct atom
					if (numModels == 1) fprintf(stderr,"> WARNING: [Model:%d,Chain:%c] Atom #%d (\"%s\") does not match atom #%d (\"%s\") from residue #%d (\"%s\")\n", numModels, chainCharId, atomId, atomName, atomCountInsideResidue, (*atomsNamesInsideResidue)[atomCountInsideResidue], residueId, residueName);
					validAtom = 0;
				}
			}
			// TODO: if a non-valid atom is found (and is not added to the list of atoms), the atoms' ids will not match the real ids from the PDB file
			if (validAtom){
				totalNumAtoms++;
				allAtomsArray[totalNumAtoms].fromResidueId = totalNumResidues;
				allAtomsArray[totalNumAtoms].atomIdInsideResidue = atomCountInsideResidue;
				allAtomsArray[totalNumAtoms].element = atomName[0];
				allAtomsArray[totalNumAtoms].coordX = coordX;
				allAtomsArray[totalNumAtoms].coordY = coordY;
				allAtomsArray[totalNumAtoms].coordZ = coordZ;
				atomCountInsideResidue++;
				numModelAtoms++;
			}
			/*
			else if(atomCountInsideResidue < numResidueAtoms){ // if not valid but still not exceeded the number of atoms, keep its space
				totalNumAtoms++;
				allAtomsArray[totalNumAtoms].fromResidueId = 0;
				allAtomsArray[totalNumAtoms].atomIdInsideResidue = 0;
				allAtomsArray[totalNumAtoms].element = '\0';
				allAtomsArray[totalNumAtoms].coordX = 0.0;
				allAtomsArray[totalNumAtoms].coordY = 0.0;
				allAtomsArray[totalNumAtoms].coordZ = 0.0;
			}
			*/
			prevAtomId = atomId;
			prevResidueId = residueId;
			prevChainCharId = chainCharId;
			continue;
		} // end of "ATOM" record
		else if (c == 'A' || c == 'T' || c == 'M' || c == 'E' || c == 'C'){ // "ANISOU", "TER", "MODEL", "MASTER", "ENDMDL", "END", "CONECT"
			while (c != EOF && c != '\n') c = fgetc(pdbFile);
			continue;
		}
		else {
			fprintf(stderr,"> ERROR: Unknown record type \"%s\"\n", recordType);
			return 0;
		}
	} // end of processing file
	if (numChains == 0 || numModels == 0){
		fprintf(stderr,"> ERROR: No chains/models found\n");
		return 0;
	}
	if (numPdbAtoms == 0){
		fprintf(stderr,"> ERROR: No atoms found\n");
		return 0;
	}
	allPDBsArray[pdbFileId].numAtoms = numPdbAtoms;
	allPDBsArray[pdbFileId].numResidues = numPdbResidues;
	allPDBsArray[pdbFileId].numChains = numChains;
	allPDBsArray[pdbFileId].numModels = numModels;
	allResiduesArray = (Residue *)realloc(allResiduesArray, (totalNumResidues + 1)*sizeof(Residue)); // +1 to account for the unused 0-th pos
	allAtomsArray = (Atom *)realloc(allAtomsArray, (totalNumAtoms + 1)*sizeof(Atom));
	printf("\t(%d chains, %d models, %d residues/model, %d atoms/model)", numChains, numModels, (numPdbResidues/numModels), (numPdbAtoms/numModels));
	//if (numWarnings > maxNumWarnings) printf("\t(%d warnings)", numWarnings);
	printf(" OK\n");
	fflush(stdout);
	if (allChainsCharIds != NULL) free(allChainsCharIds);
	if (allChainsNumResidues != NULL) free(allChainsNumResidues);
	return numPdbAtoms;
}


void SaveAtomsToJson() {
	unsigned int atomId, realNumAtoms;
	int elemId, resId, molId;
	float minX, minY, minZ, maxX, maxY, maxZ;
	FILE *jsonFile;
	if ((jsonFile = fopen("PDBAtoms.json", "w")) == NULL) {
		fprintf(stderr, "> ERROR: Cannot create JSON file\n");
		return (-1);
	}
	minX = FLT_MAX;
	minY = FLT_MAX;
	minZ = FLT_MAX;
	maxX = -(FLT_MAX);
	maxY = -(FLT_MAX);
	maxZ = -(FLT_MAX);
	realNumAtoms = 0;
	fprintf(jsonFile, "var pdbAtoms = [\n");
	for (atomId = 1; atomId <= totalNumAtoms; atomId++) {
		if (allAtomsArray[atomId].element == '\0') continue;
		resId = (int)allAtomsArray[atomId].fromResidueId;
		if ((atomId != totalNumAtoms) && (allResiduesArray[resId].fromModelId) != (allResiduesArray[(allAtomsArray[(atomId + 1)].fromResidueId)].fromModelId)) break; // only get first model
		if (realNumAtoms != 0) fprintf(jsonFile, ",\n");
		realNumAtoms++;
		molId = (int)allResiduesArray[resId].moleculeId;
		if (allAtomsArray[atomId].coordX < minX) minX = allAtomsArray[atomId].coordX;
		if (allAtomsArray[atomId].coordX > maxX) maxX = allAtomsArray[atomId].coordX;
		if (allAtomsArray[atomId].coordY < minY) minY = allAtomsArray[atomId].coordY;
		if (allAtomsArray[atomId].coordY > maxY) maxY = allAtomsArray[atomId].coordY;
		if (allAtomsArray[atomId].coordZ < minZ) minZ = allAtomsArray[atomId].coordZ;
		if (allAtomsArray[atomId].coordZ > maxZ) maxZ = allAtomsArray[atomId].coordZ;
		fprintf(jsonFile, "{ \"atomId\": %u, \"atomName\": \"%s\", ", realNumAtoms, MoleculesList[molId].atomsCodes[allAtomsArray[atomId].atomIdInsideResidue]);
		fprintf(jsonFile, "\"resName\": \"%s\", \"resId\": %d, ", MoleculesList[molId].abbrev, resId);
		fprintf(jsonFile, "\"x\": %.3f, \"y\": %.3f, \"z\": %.3f, ", allAtomsArray[atomId].coordX, allAtomsArray[atomId].coordY, allAtomsArray[atomId].coordZ);
		fprintf(jsonFile, "\"atomElem\": \"%c\"}", allAtomsArray[atomId].element);
	}
	fprintf(jsonFile, "\n];\n");
	fprintf(jsonFile, "var numPdbAtoms = %u;\n", realNumAtoms);
	fprintf(jsonFile, "var minCoords = [%.3f, %.3f, %.3f];\n", minX, minY, minZ);
	fprintf(jsonFile, "var maxCoords = [%.3f, %.3f, %.3f];\n", maxX, maxY, maxZ);
	if (fclose(jsonFile) == EOF) {
		fprintf(stderr, "> ERROR: Cannot write JSON file\n");
		return (-1);
	}
	return 0;
}

// TODO: support loading from multiple PDB files; separated atoms array per model; 
int main(int argc, char *argv[]){
	int i;
	printf("[ MDTest v%s ]\n", VERSION);
	if (argc<2){
		printf("\nUsage:\n");
		printf("\t%s <options> <PDB_file(s)>\n", argv[0]);
		printf("\nOptions:\n");
		printf("\t-s\tGet bond length statistics\n");
		printf("\t-f\tGet interatomic force statistics\n");
		printf("\t-i\tCreate SVG image\n");
		printf("\t-j\tConvert PDB to JSON format\n");
		printf("\n");
		return (-1);
	}
	numPDBs = 0;
	allPDBsArray = NULL;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			allPDBsArray = (PDBInfo *)realloc(allPDBsArray, (numPDBs + 1) * sizeof(PDBInfo));
			allPDBsArray[numPDBs].fileName = argv[i];
			allPDBsArray[numPDBs].numAtoms = 0;
			allPDBsArray[numPDBs].numResidues = 0;
			allPDBsArray[numPDBs].numChains = 0;
			allPDBsArray[numPDBs].numModels = 0;
			allPDBsArray[numPDBs].startPositionInAtomsArray = 0;
			allPDBsArray[numPDBs].startPositionInResiduesArray = 0;
			numPDBs++;
		}
	}
	if (numPDBs == 0) ExitErrorMessage("Please provide at least one PDB file name", 0);
	allAtomsArray = (Atom *)calloc(1, sizeof(Atom));
	allResiduesArray = (Residue *)calloc(1, sizeof(Residue));
	maxAtomsArraySize = 1;
	maxResiduesArraySize = 1;
	totalNumAtoms = 0;
	totalNumResidues = 0;
	for (i = 0; i < numPDBs; i++) LoadPDBFile(i,1,1);
	// TODO: set better check than the number of atoms
	if(totalNumAtoms==0) ExitErrorMessage("No valid atoms were found", 0);
	if (argv[1][0] == '-') {
		if ((argv[1][1] == 'S') || (argv[1][1] == 's')) GetBondLengthStats();
		if ((argv[1][1] == 'F') || (argv[1][1] == 'f')) GetInteratomicForceStats();
		if ((argv[1][1] == 'I') || (argv[1][1] == 'i')) DrawAtoms();
		if ((argv[1][1] == 'J') || (argv[1][1] == 'j')) SaveAtomsToJson();
	}
	printf("> Done!\n");
	if (allAtomsArray != NULL) free(allAtomsArray);
	if (allResiduesArray != NULL) free(allResiduesArray);
	if (allPDBsArray != NULL) free(allPDBsArray);
	#if PAUSE_AT_EXIT == 1
	getchar();
	#endif
	return 0;
}
