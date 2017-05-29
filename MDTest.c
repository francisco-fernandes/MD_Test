#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "ChemicalEntities.h"

#define VERSION "0.1"
#define PAUSE_AT_EXIT 1
//#define CHECK_MOLECULES 1
#define ARRAY_GROW_SIZE 64

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

typedef struct PDBInfo {
	char *fileName;
	unsigned int numAtoms;
	unsigned int numResidues;
	unsigned int numChains;
	unsigned int numModels;
	unsigned int startPositionInAtomsArray;
	unsigned int startPositionInResiduesArray;
} PDBInfo;

PDBInfo *allPDBsArray = NULL;
int numPDBs = 0;

// NOTE: position 0 in both arrays bellow is not used
Atom *allAtomsArray = NULL;
Residue *allResiduesArray = NULL;
unsigned int totalNumAtoms = 0;
unsigned int totalNumResidues = 0;
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

int GetElementCode(const char *atomName){
	return (int)( ElementIdLookupTable[ CHAR2ID(atomName[0]) ] );
}

int GetMoleculeCode(char *moleculeName){
	static const int charPos[3] = { 0, 2, 1 };
	int i = 0;
	int code = 0;
	do {
		code = (-code);
		code = (int)( MoleculeIdLookupTable[ code ][ CHAR2ID(moleculeName[charPos[i]]) ] );
		i++;
	} while (code < 0);
	return code;
}

// TODO: allow atoms' names of residues to be given non-sequentially, and perform extra sanity checks (missing atoms, repeated atoms)
//		 by first implementing GetAtomIdInsideResidueFromAtomName()
// TODO: when atoms are missing (e.g hydrogens in X-ray) fill data with position based on angles, etc.
int LoadPDBFile(int pdbFileId, int showWarnings){
	FILE *pdbFile;
	char c, prevc, pdbID[5], recordType[7];
	int i, numItemsRead;
	int numModels, numChains, *allChainsNumResidues, numChainResidues, numModelAtoms, numPdbAtoms, numPdbResidues;
	char *allChainsCharIds;
	int atomId, residueId, prevAtomId, prevResidueId, chainId;
	char atomName[5], residueName[4], chainCharId, elementSymbol[3], prevChainCharId;
	float coordX, coordY, coordZ;
	int residueCode, numResidueAtoms, atomCountInsideResidue;
	const char(*atomsNamesInsideResidue)[MAX_NUM_ATOMS_PER_AA][5];
	int validAtom;
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
	prevResidueId = INT_MAX;
	prevChainCharId = '\0';
	residueCode = 0;
	numResidueAtoms = 0;
	atomCountInsideResidue = 0;
	atomsNamesInsideResidue = NULL;
	validAtom = 0;
	while (1){ // process each line
		for (i = 0; i < 6; i++){
			c = fgetc(pdbFile);
			recordType[i] = c;
		}
		recordType[6] = '\0';
		if (c == EOF){ // end of file, save last chain and model stats
			if ((prevResidueId != 0) && (atomCountInsideResidue != numResidueAtoms)){ // check last residue in file
				if (atomCountInsideResidue < numResidueAtoms){ // less atoms
					if(showWarnings) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has less atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
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
					if (showWarnings) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has more atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
				}
			}
			if (chainId != (-1) && allChainsNumResidues[chainId] == 0) allChainsNumResidues[chainId] = numChainResidues;
			numPdbAtoms += numModelAtoms;
			break;
		} // end of file check
		c = recordType[0]; // check type of each line
		if ((c == 'A' || c == 'H') && (recordType[1]!='N')){ // "ATOM", "HETATM" (but not "ANISOU")
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
				if (prevResidueId != 0){ // check previous residue (if this is not the 1st atom in the file)
					if (atomCountInsideResidue != numResidueAtoms){ // different number of atoms
						if (atomCountInsideResidue < numResidueAtoms){ // less atoms
							if (showWarnings) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has less atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
							// TODO: when atom array filling is not sequencial, check missing atoms starting at the position of the 1st atom of the residue
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
							if (showWarnings) fprintf(stderr,"> WARNING: Residue #%d (\"%s\") has more atoms than expected (%d out of %d)\n", prevResidueId, (MoleculesList[residueCode].abbrev), atomCountInsideResidue, numResidueAtoms);
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
				if (atomCountInsideResidue < numResidueAtoms) validAtom = 1;
				else validAtom = 0; // more atoms than expected
			}
			if (validAtom && !EqualString(atomName, (*atomsNamesInsideResidue)[atomCountInsideResidue])){ // invalid atom name
				//if (atomName[0] != ((*atomsNamesInsideResidue)[atomCountInsideResidue][0])){ // compare 1st char of atoms' names
				if (residueId == 1 && EqualString(atomName, "H1")){ // atom "H" from 1st residue might be named "H1"
					if (showWarnings) fprintf(stderr,"> INFO: Atom #%d from residue #%d (\"%s\") is named \"%s\"\n", atomId, residueId, residueName, atomName);
				}
				else if (EqualString(atomName, "OXT") || EqualString(atomName, "HXT")){ // terminal oxygen or hydrogen might be present
					if (showWarnings) fprintf(stderr,"> INFO: Residue #%d (\"%s\") has a terminal atom #%d \"%s\"\n", residueId, residueName, atomId, atomName);
					// TODO: if an extra 1st or last atom exists, it will mess the hardcoded atoms bonds ids
					validAtom = 0;
					atomCountInsideResidue--;
					numModelAtoms--;
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
					fprintf(stderr,"> WARNING: Atom #%d (\"%s\") does not match atom #%d (\"%s\") from residue #%d (\"%s\")\n", atomId, atomName, atomCountInsideResidue, (*atomsNamesInsideResidue)[atomCountInsideResidue], residueId, residueName);
					validAtom = 0;
				}
			}
			if (validAtom){
				totalNumAtoms++;
				allAtomsArray[totalNumAtoms].fromResidueId = totalNumResidues;
				allAtomsArray[totalNumAtoms].atomIdInsideResidue = atomCountInsideResidue;
				allAtomsArray[totalNumAtoms].element = atomName[0];
				allAtomsArray[totalNumAtoms].coordX = coordX;
				allAtomsArray[totalNumAtoms].coordY = coordY;
				allAtomsArray[totalNumAtoms].coordZ = coordZ;
			}
			else if(atomCountInsideResidue < numResidueAtoms){ // if not valid but still not exceeded the number of atoms, keep its space
				totalNumAtoms++;
				allAtomsArray[totalNumAtoms].fromResidueId = 0;
				allAtomsArray[totalNumAtoms].atomIdInsideResidue = 0;
				allAtomsArray[totalNumAtoms].element = '\0';
				allAtomsArray[totalNumAtoms].coordX = 0.0;
				allAtomsArray[totalNumAtoms].coordY = 0.0;
				allAtomsArray[totalNumAtoms].coordZ = 0.0;
			}
			atomCountInsideResidue++;
			numModelAtoms++;
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
	printf("\t(%d chains, %d models, %d atoms) OK\n", numChains, numModels, numPdbAtoms);
	fflush(stdout);
	if (allChainsCharIds != NULL) free(allChainsCharIds);
	if (allChainsNumResidues != NULL) free(allChainsNumResidues);
	return numPdbAtoms;
}

typedef struct BondStats {
	signed char moleculeId;
	signed char srcAtomId;
	signed char destAtomId;
	float minLength;
	float maxLength;
	double avgLength;
	unsigned int numOccurrences;
} BondStats;

typedef struct ElementBondStats {
	signed char srcElementId;
	signed char destElementId;
	signed char scrElementNumBonds;
	signed char destElementNumBonds;
	float minLength;
	float maxLength;
	double avgLength;
	unsigned int numOccurrences;
	signed char minAtMoleculeId;
	signed char minAtSrcAtomId;
	signed char minAtDestAtomId;
	signed char maxAtMoleculeId;
	signed char maxAtSrcAtomId;
	signed char maxAtDestAtomId;
} ElementBondStats;

// TODO: for element bond stats, besides the number of atoms the current atom is connected to,
//		 also consider which elements are those connected atoms
void GetBondLengthStats() {
	unsigned int numAtoms;
	unsigned int numMoleculeAtoms;
	unsigned int nm, na, nb, nrb, nr, nba;
	const signed char(*bondsArray)[MAX_NUM_ATOMS_PER_AA][MAX_NUM_BONDS_PER_ATOM];
	signed char bondedAtomId;
#ifdef CHECK_MOLECULES
	unsigned char elementId;
	char atomElement, bondedAtomElement;
	char atomPlacement, bondedAtomPlacement;
	const unsigned char numBondsPerAtomType[(NUM_ELEMENTS + 1)] = { 0,1,4,2,3,2,4 }; // _,H,C,O,N,S,P
	unsigned char inBondsCount[MAX_NUM_ATOMS_PER_AA], outBondsCount[MAX_NUM_ATOMS_PER_AA];
#endif
	BondStats *bondStatsArray;
	ElementBondStats ****elementBondStatsArray, *elementBondStatsEntry;
	signed short int atomsPositionInBondStatsArray[(NUM_MOLECULES + 1)][MAX_NUM_ATOMS_PER_AA];
	int totalNumValidBonds, numValidAtomBonds, maxArraySize;
	float bondLength, coordDist;
	totalNumValidBonds = 0;
	maxArraySize = 0;
	bondStatsArray = NULL;
#ifdef CHECK_MOLECULES
	printf("> Checking consistency of %d default molecules...\n", NUM_MOLECULES);
	fflush(stdout);
#endif
	for (nm = 1; nm <= NUM_MOLECULES; nm++) { // loop for all molecules
		numMoleculeAtoms = (int)(MoleculesList[nm].numberOfAtoms);
		bondsArray = &(MoleculesList[nm].atomsBonds);
		for (na = 0; na < MAX_NUM_ATOMS_PER_AA; na++) atomsPositionInBondStatsArray[nm][na] = 0; // reset stats array entries
#ifdef CHECK_MOLECULES
		for (na = 0; na < numMoleculeAtoms; na++) { // reset in/out counts
			inBondsCount[na] = 0;
			outBondsCount[na] = 0;
		}
#endif
		for (na = 0; na < numMoleculeAtoms; na++) { // loop for all molecule atoms
#ifdef CHECK_MOLECULES
			atomElement = GET_ELEMENT_LETTER(nm, na);
			atomPlacement = GET_PLACEMENT_ORDER(nm, na);
#endif
			numValidAtomBonds = 0;
			for (nb = 0; nb < MAX_NUM_BONDS_PER_ATOM; nb++) { // loop for all atom bonds
				bondedAtomId = (*bondsArray)[na][nb];
				if (bondedAtomId == SCHAR_MAX) break; // no more bonds
#ifdef CHECK_MOLECULES
				if (bondedAtomId == (-1) && na == 0 && nb == 0) continue; // edge atom "N" bonds to previous molecule
				if (bondedAtomId == na) { // check atom bonded to itself
					printf("> WARNING: Bond #%d of atom \"%s\" of molecule \"%s\" links to itself.\n", nb, MoleculesList[nm].atomsCodes[na], MoleculesList[nm].abbrev);
				}
				if (bondedAtomId >= numMoleculeAtoms) { // check out-of-bounds destination atom
					// TODO: if checking molecules, the bond stats of the last atom will not be processed
					if (bondedAtomId == numMoleculeAtoms && na == 2 && nb == 2) continue; // exception is edge atom "C" which bonds to next molecule
					printf("> WARNING: Bond #%d of atom \"%s\" of molecule \"%s\" is outside molecule (%d>%d).\n", nb, MoleculesList[nm].atomsCodes[na], MoleculesList[nm].abbrev, bondedAtomId, numMoleculeAtoms);
				}
				for (nrb = 0; nrb < MAX_NUM_BONDS_PER_ATOM; nrb++) { // check if reverse bond exists
					if ((*bondsArray)[bondedAtomId][nrb] == na) break;
				}
				if (nrb == MAX_NUM_BONDS_PER_ATOM) {
					printf("> WARNING: Reverse bond \"%s\"<-\"%s\" of molecule \"%s\" does not exist.\n", MoleculesList[nm].atomsCodes[na], MoleculesList[nm].atomsCodes[bondedAtomId], MoleculesList[nm].abbrev);
				}
				outBondsCount[na]++; // update in/out counts
				inBondsCount[bondedAtomId]++;
				bondedAtomElement = GET_ELEMENT_LETTER(nm, bondedAtomId);
				bondedAtomPlacement = GET_PLACEMENT_ORDER(nm, bondedAtomId);
				if (atomElement == 'H' || bondedAtomElement == 'H') { // check placement of bonded H atom
					if (atomElement == 'H' && bondedAtomElement == 'H') { // check H bonded to H
						printf("> WARNING: Bond of Hydrogens \"%s\"->\"%s\" of molecule \"%s\" is invalid.\n", MoleculesList[nm].atomsCodes[na], MoleculesList[nm].atomsCodes[bondedAtomId], MoleculesList[nm].abbrev);
					}
					if (atomPlacement != bondedAtomPlacement) {
						printf("> WARNING: Bond \"%s\"->\"%s\" of molecule \"%s\" does not respect placement.\n", MoleculesList[nm].atomsCodes[na], MoleculesList[nm].atomsCodes[bondedAtomId], MoleculesList[nm].abbrev);
					}
				}
				else { // check placement of bonded non-H atom
					if (((bondedAtomPlacement != atomPlacement) &&
						((bondedAtomId < na && bondedAtomPlacement != (atomPlacement - 1)
							&& (na != 2 || bondedAtomId != 1) // exception are atoms "C" and "CA" which have switched ids/placements
							&& (nm != 15 || (na != 6 || bondedAtomId != 0))) // and also atoms "N" and "CD" of "Proline"
							|| (bondedAtomId > na && bondedAtomPlacement != (atomPlacement + 1)
								&& (na != 1 || bondedAtomId != 2) // exception are atoms "CA" and "C" which have switched ids/placements
								&& (nm != 15 || (na != 0 || bondedAtomId != 6))))) // and also atoms "CD" and "N" of "Proline"
						|| ((bondedAtomPlacement == atomPlacement) &&
							!((atomPlacement == 0) && ((bondedAtomId == 2 && na == 3) || (bondedAtomId == 3 && na == 2))) // exception are atoms "C" and "O" which bond and have the same placement
							&& !((atomPlacement == 5) && (nm == 9) && ((bondedAtomId == 8 && na == 9) || (bondedAtomId == 9 && na == 8))) // and also atoms "CE1" and "NE2" of "Histidine"
							&& !((atomPlacement == 5) && (nm == 18) && ((bondedAtomId == 8 && na == 9) || (bondedAtomId == 9 && na == 8))))) { // and also atoms "NE1" and "CE2" of "Tryptophan"
						printf("> WARNING: Bond \"%s\"->\"%s\" of molecule \"%s\" does not respect placement.\n", MoleculesList[nm].atomsCodes[na], MoleculesList[nm].atomsCodes[bondedAtomId], MoleculesList[nm].abbrev);
					}
				}
#endif
				if (bondedAtomId < (signed)na) continue; // valid atom only if connecting to an atom with a higher id
				if (totalNumValidBonds == maxArraySize) { // allocate more space in array if needed
					maxArraySize += ARRAY_GROW_SIZE;
					bondStatsArray = (BondStats *)realloc(bondStatsArray, maxArraySize * sizeof(BondStats));
				}
				bondStatsArray[totalNumValidBonds].destAtomId = bondedAtomId; // set destination atom and reset entries
				bondStatsArray[totalNumValidBonds].srcAtomId = na;
				bondStatsArray[totalNumValidBonds].moleculeId = nm;
				bondStatsArray[totalNumValidBonds].minLength = FLT_MAX;
				bondStatsArray[totalNumValidBonds].maxLength = FLT_MIN;
				bondStatsArray[totalNumValidBonds].avgLength = 0.0;
				bondStatsArray[totalNumValidBonds].numOccurrences = 0;
				totalNumValidBonds++;
				numValidAtomBonds++;
			} // end of bonds loop
#ifdef CHECK_MOLECULES
			elementId = ElementIdLookupTable[CHAR2ID(atomElement)];
			if (nm == 12 && na == 8) nb--; // exception is atom "NZ" of "Lysine" which has 4 bonds
			if (nb != numBondsPerAtomType[elementId] && nb != (numBondsPerAtomType[elementId] - 1)) { // check number of bonds for this element
				printf("> WARNING: Atom \"%s\" of molecule \"%s\" has an invalid number of bonds (%d).\n", MoleculesList[nm].atomsCodes[na], MoleculesList[nm].abbrev, nb);
			}
#endif
			if (numValidAtomBonds != 0) atomsPositionInBondStatsArray[nm][na] = (totalNumValidBonds - numValidAtomBonds); // start position of atom's bonds in array
			else atomsPositionInBondStatsArray[nm][na] = (-1);
		} // end of atoms loop
#ifdef CHECK_MOLECULES
		for (na = 0; na < numMoleculeAtoms; na++) { // check if the number of in/out bonds is the same
			if (outBondsCount[na] != inBondsCount[na]) {
				printf("> WARNING: Number of inward (%d) and outwards (%d) bonds of atom \"%s\" of molecule \"%s\" does not match.\n", inBondsCount[na], outBondsCount[na], MoleculesList[nm].atomsCodes[na], MoleculesList[nm].abbrev);
			}
		}
#endif
	} // end of molecules loop
	bondStatsArray = (BondStats *)realloc(bondStatsArray, totalNumValidBonds * sizeof(BondStats));
	for (numAtoms = 1; numAtoms <= totalNumAtoms; numAtoms++) {
		nr = allAtomsArray[numAtoms].fromResidueId;
		if (nr == 0) continue;
		nm = allResiduesArray[nr].moleculeId;
		bondsArray = &(MoleculesList[nm].atomsBonds);
		na = allAtomsArray[numAtoms].atomIdInsideResidue;
		for (nb = 0; nb < MAX_NUM_BONDS_PER_ATOM; nb++) { // check all atoms that this atom links to
			bondedAtomId = (*bondsArray)[na][nb];
			if (bondedAtomId == SCHAR_MAX) break;
			if (bondedAtomId < (signed)na) continue;
			nba = (allResiduesArray[nr].startPositionInAtomsArray + bondedAtomId);
			if (bondedAtomId >= MoleculesList[nm].numberOfAtoms) { // if it links to the atom in the next residue, check if that atom is valid
				if (nba >= totalNumAtoms) continue;
				if (allResiduesArray[nr].fromChain != allResiduesArray[(allAtomsArray[nba].fromResidueId)].fromChain) continue;
				if (allResiduesArray[nr].fromModelId != allResiduesArray[(allAtomsArray[nba].fromResidueId)].fromModelId) continue;
			}
			if (allAtomsArray[nba].fromResidueId == 0) continue; // check if atom instance in all atoms array is valid
			nrb = atomsPositionInBondStatsArray[nm][na]; // start of this molecule's atom stats in stats array
			while (bondStatsArray[nrb].destAtomId != bondedAtomId) nrb++; // find location of this bond stats in stats array
			bondLength = 0.0; // calculate 3D square root
			coordDist = (allAtomsArray[numAtoms].coordX - allAtomsArray[nba].coordX);
			bondLength += coordDist*coordDist;
			coordDist = (allAtomsArray[numAtoms].coordY - allAtomsArray[nba].coordY);
			bondLength += coordDist*coordDist;
			coordDist = (allAtomsArray[numAtoms].coordZ - allAtomsArray[nba].coordZ);
			bondLength += coordDist*coordDist;
			bondLength = sqrtf(bondLength);
			if (bondLength < bondStatsArray[nrb].minLength) bondStatsArray[nrb].minLength = bondLength;
			if (bondLength > bondStatsArray[nrb].maxLength) bondStatsArray[nrb].maxLength = bondLength;
			bondStatsArray[nrb].avgLength += (double)bondLength;
			bondStatsArray[nrb].numOccurrences++;
		}
	}
	nm = 0;
	for (nb = 0; nb < (unsigned)totalNumValidBonds; nb++) {
		if (bondStatsArray[nb].moleculeId != nm) {
			nm = bondStatsArray[nb].moleculeId;
			printf("[%02d] %s ('%s'):\n", nm, MoleculesList[nm].name, MoleculesList[nm].abbrev);
		}
		bondStatsArray[nb].avgLength /= (double)bondStatsArray[nb].numOccurrences;
		na = bondStatsArray[nb].srcAtomId;
		nba = bondStatsArray[nb].destAtomId;
		printf("\t%4s -> %4s\t", MoleculesList[nm].atomsCodes[na], (nba >= MoleculesList[nm].numberOfAtoms) ? "*" : MoleculesList[nm].atomsCodes[nba]);
		if (bondStatsArray[nb].numOccurrences != 0) {
			bondLength = (float)bondStatsArray[nb].avgLength;
			printf("[ %+.3f ; %.3f ; %+.3f]\t",
				(bondStatsArray[nb].minLength - bondLength), bondLength, (bondStatsArray[nb].maxLength - bondLength));
		}
		printf("(%u)\n", bondStatsArray[nb].numOccurrences);
	}
	elementBondStatsArray = (ElementBondStats ****)calloc((NUM_ELEMENTS + 1), sizeof(ElementBondStats ***)); // elementBondStatsArray[#SrcElem][#SrcBonds][#DestElem][#DestBonds]
	for (na = 1; na <= NUM_ELEMENTS; na++) {
		elementBondStatsArray[na] = (ElementBondStats ***)calloc((MAX_NUM_BONDS_PER_ATOM + 1), sizeof(ElementBondStats **));
		for (nb = 0; nb <= MAX_NUM_BONDS_PER_ATOM; nb++) {
			elementBondStatsArray[na][nb] = (ElementBondStats **)calloc((NUM_ELEMENTS + 1), sizeof(ElementBondStats *));
			for (nr = 1; nr <= NUM_ELEMENTS; nr++) {
				elementBondStatsArray[na][nb][nr] = (ElementBondStats *)calloc((MAX_NUM_BONDS_PER_ATOM + 1), sizeof(ElementBondStats));
				for (nrb = 0; nrb <= MAX_NUM_BONDS_PER_ATOM; nrb++) {
					elementBondStatsArray[na][nb][nr][nrb].minLength = FLT_MAX;
					elementBondStatsArray[na][nb][nr][nrb].maxLength = FLT_MIN;
					elementBondStatsArray[na][nb][nr][nrb].avgLength = 0.0;
					elementBondStatsArray[na][nb][nr][nrb].numOccurrences = 0;
				}
			}
		}
	}
	for (nba = 0; nba < (unsigned)totalNumValidBonds; nba++) {
		nm = bondStatsArray[nba].moleculeId;
		na = bondStatsArray[nba].srcAtomId;
		nr = bondStatsArray[nba].destAtomId;
		if (nr == (-1) || nr == (MoleculesList[nm].numberOfAtoms) || (bondStatsArray[nba].numOccurrences) == 0) continue;
		bondLength = (float)bondStatsArray[nba].avgLength;
		for (nb = 0; nb < MAX_NUM_BONDS_PER_ATOM; nb++) if (MoleculesList[nm].atomsBonds[na][nb] == SCHAR_MAX) break;
		for (nrb = 0; nrb < MAX_NUM_BONDS_PER_ATOM; nrb++) if (MoleculesList[nm].atomsBonds[nr][nrb] == SCHAR_MAX) break;
		na = GetElementCode(MoleculesList[nm].atomsCodes[na]);
		nr = GetElementCode(MoleculesList[nm].atomsCodes[nr]);
		elementBondStatsEntry = &(elementBondStatsArray[na][nb][nr][nrb]); // first set the stats for the source element
		elementBondStatsEntry->srcElementId = na;
		elementBondStatsEntry->destElementId = nr;
		elementBondStatsEntry->scrElementNumBonds = nb;
		elementBondStatsEntry->destElementNumBonds = nrb;
		elementBondStatsEntry->avgLength += bondLength;
		elementBondStatsEntry->numOccurrences++;
		if (bondLength < elementBondStatsEntry->minLength) {
			elementBondStatsEntry->minLength = bondLength;
			elementBondStatsEntry->minAtMoleculeId = nm;
			elementBondStatsEntry->minAtSrcAtomId = bondStatsArray[nba].srcAtomId;
			elementBondStatsEntry->minAtDestAtomId = bondStatsArray[nba].destAtomId;
		}
		if (bondLength > elementBondStatsEntry->maxLength) {
			elementBondStatsEntry->maxLength = bondLength;
			elementBondStatsEntry->maxAtMoleculeId = nm;
			elementBondStatsEntry->maxAtSrcAtomId = bondStatsArray[nba].srcAtomId;
			elementBondStatsEntry->maxAtDestAtomId = bondStatsArray[nba].destAtomId;
		}
		elementBondStatsEntry = &(elementBondStatsArray[nr][nrb][na][nb]); // now set the stats for the destination element
		elementBondStatsEntry->srcElementId = nr;
		elementBondStatsEntry->destElementId = na;
		elementBondStatsEntry->scrElementNumBonds = nrb;
		elementBondStatsEntry->destElementNumBonds = nb;
		elementBondStatsEntry->avgLength += bondLength;
		elementBondStatsEntry->numOccurrences++;
		if (bondLength < elementBondStatsEntry->minLength) {
			elementBondStatsEntry->minLength = bondLength;
			elementBondStatsEntry->minAtMoleculeId = nm;
			elementBondStatsEntry->minAtSrcAtomId = bondStatsArray[nba].destAtomId;
			elementBondStatsEntry->minAtDestAtomId = bondStatsArray[nba].srcAtomId;
		}
		if (bondLength > elementBondStatsEntry->maxLength) {
			elementBondStatsEntry->maxLength = bondLength;
			elementBondStatsEntry->maxAtMoleculeId = nm;
			elementBondStatsEntry->maxAtSrcAtomId = bondStatsArray[nba].destAtomId;
			elementBondStatsEntry->maxAtDestAtomId = bondStatsArray[nba].srcAtomId;
		}
	}
	for (na = 1; na <= NUM_ELEMENTS; na++) {
		for (nr = 0; nr <= MAX_NUM_BONDS_PER_ATOM; nr++) {
			for (nb = 1; nb <= NUM_ELEMENTS; nb++) {
				for (nrb = 0; nrb <= MAX_NUM_BONDS_PER_ATOM; nrb++) {
					elementBondStatsEntry = &(elementBondStatsArray[na][nr][nb][nrb]);
					if ((elementBondStatsEntry->numOccurrences) == 0) continue;
					printf("%c(%d) <-> %c(%d) : ", ElementsList[(elementBondStatsEntry->srcElementId)].letter , (elementBondStatsEntry->scrElementNumBonds),
						ElementsList[(elementBondStatsEntry->destElementId)].letter, (elementBondStatsEntry->destElementNumBonds));
					(elementBondStatsEntry->avgLength) /= (double)(elementBondStatsEntry->numOccurrences);
					bondLength = (float)(elementBondStatsEntry->avgLength);
					nm = (elementBondStatsEntry->minAtMoleculeId);
					printf("(%3s: %4s - %-4s ) ", MoleculesList[nm].abbrev, MoleculesList[nm].atomsCodes[(elementBondStatsEntry->minAtSrcAtomId)], MoleculesList[nm].atomsCodes[(elementBondStatsEntry->minAtDestAtomId)]);
					printf("[ %+.3f ; %.3f ; %+.3f ]", ((elementBondStatsEntry->minLength) - bondLength), bondLength, ((elementBondStatsEntry->maxLength) - bondLength));
					nm = (elementBondStatsEntry->maxAtMoleculeId);
					printf(" (%3s: %4s - %-4s )\n", MoleculesList[nm].abbrev, MoleculesList[nm].atomsCodes[(elementBondStatsEntry->maxAtSrcAtomId)], MoleculesList[nm].atomsCodes[(elementBondStatsEntry->maxAtDestAtomId)]);
				}
			}
		}
	}
	for (na = 1; na <= NUM_ELEMENTS; na++) {
		for (nb = 0; nb <= MAX_NUM_BONDS_PER_ATOM; nb++) {
			for (nr = 1; nr <= NUM_ELEMENTS; nr++) {
				free(elementBondStatsArray[na][nb][nr]);
			}
			free(elementBondStatsArray[na][nb]);
		}
		free(elementBondStatsArray[na]);
	}
	free(elementBondStatsArray);
	free(bondStatsArray);
}

// TODO: test with charges for carbon atoms too
// TODO: test with and without considering bonded atoms to current atom
// TODO: test with and without considering the charge of the current atom
// TODO: try quadratic distance
// TODO: after initial estimation, refine for all atoms types (e.g. with different number of bonds) with a range around those values 
void GetInteratomicForceStats() {

}

// TODO: support loading from multiple PDB files
int main(int argc, char *argv[]){
	int i;
	printf("[ MDTest v%s ]\n", VERSION);
	if (argc<2){
		printf("\nUsage:\n");
		printf("\t%s <options> <PDB_files>\n", argv[0]);
		printf("\nOptions:\n");
		printf("\t-s\tGet bond length statistics\n");
		printf("\t-f\tGet interatomic force statistics\n");
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
	for (i = 0; i < numPDBs; i++) LoadPDBFile(i, 1);
	if (argv[1][0] == '-') {
		if ((argv[1][1] == 'S') || (argv[1][1] == 's')) GetBondLengthStats();
		if ((argv[1][1] == 'F') || (argv[1][1] == 'f')) GetInteratomicForceStats();
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
