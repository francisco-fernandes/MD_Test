#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "Statistics.h"
#include "AtomsArray.h"
#include "SVGGraphics.h"
#include "ChemicalEntities.h"

//#define CHECK_MOLECULES 1


typedef struct BondStats {
	signed char moleculeId;
	signed char srcAtomId;
	signed char destAtomId;
	float minLength;
	float maxLength;
	double avgLength;
	unsigned int numOccurrences;
} BondStats;

// TODO: implement visualization for this
// TODO: get statistics for the range of distances between a min and a max, and then draw an histogram
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
			bondLength += coordDist * coordDist;
			coordDist = (allAtomsArray[numAtoms].coordY - allAtomsArray[nba].coordY);
			bondLength += coordDist * coordDist;
			coordDist = (allAtomsArray[numAtoms].coordZ - allAtomsArray[nba].coordZ);
			bondLength += coordDist * coordDist;
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
					printf("%c(%d) <-> %c(%d) : ", ElementsList[(elementBondStatsEntry->srcElementId)].letter, (elementBondStatsEntry->scrElementNumBonds),
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
// TODO: 3D scatter plot with transparent points with variable color and radius
// TODO: test "non-spherical" forces, based on direction and strength according to the closeness of the atoms (initial value, drained by nearest atoms in order, when drained, <=0 , no more interactions)
void GetInteratomicForceStats() {
	int elemId, otherElemId;
	unsigned int i, j, atomId, otherAtomId;
	unsigned int numSteps, totalNumSteps, progressMark;
	unsigned int numChargedAtoms;
	ChargedAtom *chargedAtomsArray;
	double globalEnergy, minEnergy, maxEnergy, minAbsEnergy;
	float signValue[(NUM_ELEMENTS + 1)] = { 0.0f , +1.0f , 0.0f , -1.0f , +1.0f , 0.0f , 0.0f }; // -,[H],C,[O],[N],S,P
	float startValue[(NUM_ELEMENTS + 1)] = { 0.0f ,  2.0f , 0.0f ,  2.0f ,  2.0f , 0.0f , 0.0f }; // -,[H],C,[O],[N],S,P
	float endValue[(NUM_ELEMENTS + 1)] = { 0.0f , 10.0f , 0.0f , 10.0f , 10.0f , 0.0f , 0.0f }; // -,[H],C,[O],[N],S,P
	float stepValue[(NUM_ELEMENTS + 1)] = { 0.0f ,  1.0f , 0.0f ,  1.0f ,  1.0f , 0.0f , 0.0f }; // -,[H],C,[O],[N],S,P
	float currValue[(NUM_ELEMENTS + 1)] = { 0.0f ,  0.0f , 0.0f ,  0.0f ,  0.0f , 0.0f , 0.0f };
	float minAbsEnergyValue[(NUM_ELEMENTS + 1)];
	float distance, coordDist, elemChargeValue;
	unsigned int k, *numDistSteps, *numElementAtoms, ***distancesPerElementPair;
	time_t time_start, time_end;
	double time_elapsed;
	ForceFieldStats *plotPoints;
	printf("> Charges to test: ");
	for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) {
		if (stepValue[elemId] == 0.0f) continue;
		printf("[%c]={%+.3f...%+.3f}(%+.3f) ; ", ElementsList[elemId].letter, (signValue[elemId] * startValue[elemId]), (signValue[elemId] * endValue[elemId]), (signValue[elemId] * stepValue[elemId]));
	}
	printf("\n");
	printf("> Identifying charged atoms... ");
	fflush(stdout);
	numChargedAtoms = 0;
	for (atomId = 1; atomId <= totalNumAtoms; atomId++) { // get total number of atoms whose charge will be tested
		if (allAtomsArray[atomId].element == '\0') continue;
		if ((atomId != totalNumAtoms) && (allResiduesArray[(allAtomsArray[atomId].fromResidueId)].fromModelId) != (allResiduesArray[(allAtomsArray[(atomId + 1)].fromResidueId)].fromModelId)) break; // only get first model
		elemId = GetAtomElementCode(&allAtomsArray[atomId]);
		if (stepValue[elemId] == 0.0f) continue; // only H, O and N
		numChargedAtoms++;
	}
	chargedAtomsArray = (ChargedAtom *)malloc(numChargedAtoms * sizeof(ChargedAtom));
	numElementAtoms = (unsigned int *)calloc((NUM_ELEMENTS + 1), sizeof(unsigned int));
	numChargedAtoms = 0;
	for (atomId = 1; atomId <= totalNumAtoms; atomId++) { // initialize arrays with the ids, elements and charges of all the relevant atoms
		if (allAtomsArray[atomId].element == '\0') continue;
		if ((atomId != totalNumAtoms) && (allResiduesArray[(allAtomsArray[atomId].fromResidueId)].fromModelId) != (allResiduesArray[(allAtomsArray[(atomId + 1)].fromResidueId)].fromModelId)) break;
		elemId = GetAtomElementCode(&allAtomsArray[atomId]);
		if (stepValue[elemId] == 0.0f) continue;
		numElementAtoms[elemId]++;
		chargedAtomsArray[numChargedAtoms].id = atomId;
		chargedAtomsArray[numChargedAtoms].elementId = elemId;
		chargedAtomsArray[numChargedAtoms].charge = 0.0f;
		numChargedAtoms++;
	}
	printf("(%u)\n", numChargedAtoms);
	printf("> Collecting distances between pairs of charged atoms ");
	fflush(stdout);
	time_start = time(NULL);
	numDistSteps = (unsigned int *)calloc((NUM_ELEMENTS + 1), sizeof(unsigned int));
	distancesPerElementPair = (unsigned int ***)calloc((NUM_ELEMENTS + 1), sizeof(unsigned int **)); // distancesPerElementPair[elemId][otherElemId][distance]
	for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // set array that, for each pair of elements, stores the number of pairs of atoms at each distance step
		if (stepValue[elemId] == 0.0f) continue;
		numDistSteps[elemId] = (int)roundf((endValue[elemId] - startValue[elemId]) / stepValue[elemId]) + 1;
		distancesPerElementPair[elemId] = (unsigned int **)calloc((NUM_ELEMENTS + 1), sizeof(unsigned int *));
		for (otherElemId = 1; otherElemId <= NUM_ELEMENTS; otherElemId++) {
			if (stepValue[otherElemId] == 0.0f) continue;
			distancesPerElementPair[elemId][otherElemId] = (unsigned int *)calloc(numDistSteps[elemId], sizeof(unsigned int));
		}
	}
	free(numDistSteps);
	totalNumSteps = (numChargedAtoms - 1) * (numChargedAtoms) / 2; // sum = n*(n+1)/2
	progressMark = (totalNumSteps / 20);
	numSteps = 0;
	for (i = 0; i < numChargedAtoms; i++) { // collect the distances between each pair of atoms in the corresponding array entry
		atomId = chargedAtomsArray[i].id;
		elemId = (int)chargedAtomsArray[i].elementId;
		for (j = (i + 1); j < numChargedAtoms; j++) {
			if (numSteps == progressMark) {
				printf(".");
				progressMark += (totalNumSteps / 20);
			}
			numSteps++;
			otherAtomId = chargedAtomsArray[j].id;
			otherElemId = (int)chargedAtomsArray[j].elementId;
			distance = 0.0f; // calculate distance between the two atoms
			coordDist = (allAtomsArray[atomId].coordX - allAtomsArray[otherAtomId].coordX);
			distance += coordDist * coordDist;
			coordDist = (allAtomsArray[atomId].coordY - allAtomsArray[otherAtomId].coordY);
			distance += coordDist * coordDist;
			coordDist = (allAtomsArray[atomId].coordZ - allAtomsArray[otherAtomId].coordZ);
			distance += coordDist * coordDist;
			distance = sqrtf(distance);
			if (distance <= endValue[elemId]) { // not further than maximum range of element
				if (distance <= startValue[elemId]) distancesPerElementPair[elemId][otherElemId][0]++; // closer than minimum range of element
				else {
					k = (int)ceilf((distance - startValue[elemId]) / stepValue[elemId]); // interval: ](step),(step+1)]
					distancesPerElementPair[elemId][otherElemId][k]++; // each entry stores the count of the distances from the previous step (exclusive) to the current step (inclusive)
				}
			}
			if (distance <= endValue[otherElemId]) { // now the other element
				if (distance <= startValue[otherElemId]) distancesPerElementPair[otherElemId][elemId][0]++;
				else {
					k = (int)ceilf((distance - startValue[otherElemId]) / stepValue[otherElemId]);
					distancesPerElementPair[otherElemId][elemId][k]++;
				}
			}
		}
	}
	//free(chargedAtomsArray); // free here when only the optimized code is running
	printf("\n");
	printf("> Calculating optimal element charge values ");
	fflush(stdout);
	totalNumSteps = 1;
	for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // get the number of charge combinations to be tested
		if (stepValue[elemId] == 0.0f) continue;
		currValue[elemId] = endValue[elemId]; // set end values, so they will be changed to the start values on the first iteration of the main loop
		totalNumSteps *= (int)roundf((endValue[elemId] - startValue[elemId]) / stepValue[elemId]) + 1; // count=(max-min+1)
	}
	plotPoints = (ForceFieldStats *)calloc(totalNumSteps, sizeof(ForceFieldStats));
	printf("(%u combinations to test) ", totalNumSteps);
	fflush(stdout);
	progressMark = (totalNumSteps / 20);
	numSteps = 0;
	minEnergy = DBL_MAX;
	maxEnergy = -(DBL_MAX);
	minAbsEnergy = DBL_MAX;
	while (1) { // iterate through all the charge values for all the elements
		if (numSteps == progressMark) {
			printf(".");
			progressMark += (totalNumSteps / 20);
		}
		numSteps++;
		for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // update charge values for this step
			if (stepValue[elemId] == 0.0f) continue;
			if (currValue[elemId] == endValue[elemId]) { // if this element has reached its maximum value
				currValue[elemId] = startValue[elemId]; // reset its counter, and go to the next element
				continue;
			}
			currValue[elemId] += stepValue[elemId]; // otherwise, simply increase its counter
			break;
		}
		if (elemId == (NUM_ELEMENTS + 1) && numSteps != 1) break; // finished processing all values of all elements
		globalEnergy = 0.0f;
		for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // calculate global energy from the sums of distances between each pair of elements
			if (stepValue[elemId] == 0.0f) continue;
			globalEnergy += numElementAtoms[elemId] * (signValue[elemId] * currValue[elemId]); // sum the "self" energy of all the atoms of this element
			for (otherElemId = 1; otherElemId <= NUM_ELEMENTS; otherElemId++) { // add the energy of the atoms of this element to all the other atoms (collected at a distance <= currValue[elemId])
				if (stepValue[otherElemId] == 0.0f) continue;
				k = (int)roundf((currValue[elemId] - startValue[elemId]) / stepValue[elemId]); // get the index of the current energy value inside the array (from 0 to (numSteps-1))
				for (i = 0; i <= k; i++) { // process all stored distances lower or equal than the current energy range
					// number of atoms at that distance times how much is left of the current energy at an average distance between this and the previous step ( = #atoms * energy-(distance-step/2) )
					globalEnergy += distancesPerElementPair[elemId][otherElemId][i] * (signValue[elemId] * (currValue[elemId] - ((startValue[elemId] + i * stepValue[elemId]) - (stepValue[elemId] / 2.0))));
				}
			}
		}
		for (i = 0; i < numChargedAtoms; i++) globalEnergy += (double)chargedAtomsArray[i].charge; // sum the charges of all the atoms
		if (globalEnergy < minEnergy) minEnergy = globalEnergy;
		if (globalEnergy > maxEnergy) maxEnergy = globalEnergy;
		if (fabs(globalEnergy) < minAbsEnergy) { // note down the current charge of each element when the minimum energy was reached
			minAbsEnergy = fabs(globalEnergy);
			for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) minAbsEnergyValue[elemId] = currValue[elemId];
		}
		plotPoints[(numSteps - 1)].globalEnergy = (float)globalEnergy; // fill point information to be drawn on plot
		i = 0;
		for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) {
			if (stepValue[elemId] == 0.0f) continue;
			if (i == 0) plotPoints[(numSteps - 1)].energyX = (signValue[elemId] * currValue[elemId]);
			else if (i == 1) plotPoints[(numSteps - 1)].energyY = (signValue[elemId] * currValue[elemId]);
			else if (i == 2) plotPoints[(numSteps - 1)].energyZ = (signValue[elemId] * currValue[elemId]);
			i++;
		}
	}
	free(numElementAtoms);
	for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) {
		if (stepValue[elemId] == 0.0f) continue;
		for (otherElemId = 1; otherElemId <= NUM_ELEMENTS; otherElemId++) {
			if (stepValue[otherElemId] == 0.0f) continue;
			free(distancesPerElementPair[elemId][otherElemId]);
		}
		free(distancesPerElementPair[elemId]);
	}
	free(distancesPerElementPair);
	time_end = time(NULL);
	time_elapsed = difftime(time_end, time_start);
	printf(" (%.0lfm %.0lfs)\n", (time_elapsed / 60.0), fmod(time_elapsed, 60.0));
	PlotInteratomicForceStats(plotPoints, (numSteps - 1));
	free(plotPoints);
	printf("\tSystem energy   : min = %+.3lf ; max = %+.3lf ; minAbsolute = %+.3lf\n", minEnergy, maxEnergy, minAbsEnergy);
	printf("\tElement charges : ");
	for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) {
		if (stepValue[elemId] == 0.0f) continue;
		printf("[%c]=%+.3f ; ", ElementsList[elemId].letter, (signValue[elemId] * minAbsEnergyValue[elemId]));
	}
	printf("\n");
	/* old slower code bellow */
	/**/
	printf("> Calculating optimal element charge values ");
	totalNumSteps = 1;
	for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // get the number of charge combinations to be tested
		if (stepValue[elemId] == 0.0f) continue;
		currValue[elemId] = endValue[elemId]; // set end values, so they will be changed to the start values on the first iteration of the main loop
		totalNumSteps *= (int)((endValue[elemId] - startValue[elemId] + stepValue[elemId]) / stepValue[elemId]); // +stepValue because count=(max-min+1)
	}
	printf("(%u combinations to test) ", totalNumSteps);
	fflush(stdout);
	time_start = time(NULL);
	numSteps = 0;
	progressMark = (totalNumSteps / 20);
	minEnergy = DBL_MAX;
	maxEnergy = -(DBL_MAX);
	minAbsEnergy = DBL_MAX;
	while (1) { // iterate through all the charge values for all the elements
		if (numSteps == progressMark) {
			printf(".");
			progressMark += (totalNumSteps / 20);
		}
		numSteps++;
		for (elemId = 1; elemId <= NUM_ELEMENTS; elemId++) { // update charge values for this step
			if (stepValue[elemId] == 0.0f) continue;
			if (currValue[elemId] == endValue[elemId]) { // if this element has reached its maximum value
				currValue[elemId] = startValue[elemId]; // reset its counter, and go to the next element
				continue;
			}
			currValue[elemId] += stepValue[elemId]; // otherwise, simply increase its counter
			break;
		}
		if (elemId == (NUM_ELEMENTS + 1) && numSteps != 1) break; // finished processing all values of all elements
		for (i = 0; i < numChargedAtoms; i++) { // reset atoms charges to their initial (current) values
			elemId = (int)chargedAtomsArray[i].elementId;
			chargedAtomsArray[i].charge = signValue[elemId] * currValue[elemId];
		}
		for (i = 0; i < numChargedAtoms; i++) { // process all charged atoms
			atomId = chargedAtomsArray[i].id;
			elemId = (int)chargedAtomsArray[i].elementId;
			for (j = (i + 1); j < numChargedAtoms; j++) { // add the influence of the current atom to all the other atoms in front of him in the array, and vice-versa
				otherAtomId = chargedAtomsArray[j].id;
				otherElemId = (int)chargedAtomsArray[j].elementId;
				distance = 0.0f; // calculate distance between the two atoms
				coordDist = (allAtomsArray[atomId].coordX - allAtomsArray[otherAtomId].coordX);
				distance += coordDist * coordDist;
				coordDist = (allAtomsArray[atomId].coordY - allAtomsArray[otherAtomId].coordY);
				distance += coordDist * coordDist;
				coordDist = (allAtomsArray[atomId].coordZ - allAtomsArray[otherAtomId].coordZ);
				distance += coordDist * coordDist;
				distance = sqrtf(distance);
				elemChargeValue = currValue[otherElemId]; // update the energy of the current atom when it is within the range of the other element's charge (positive value)
				if (distance < elemChargeValue) chargedAtomsArray[i].charge += signValue[otherElemId] * (elemChargeValue - distance);
				elemChargeValue = currValue[elemId]; // update the energy of the other atom when it is within the range of the current element's charge
				if (distance < elemChargeValue) chargedAtomsArray[j].charge += signValue[elemId] * (elemChargeValue - distance);
			}
		}
		globalEnergy = 0.0f;
		for (i = 0; i < numChargedAtoms; i++) globalEnergy += (double)chargedAtomsArray[i].charge; // sum the charges of all the atoms
		if (globalEnergy < minEnergy) minEnergy = globalEnergy;
		if (globalEnergy > maxEnergy) maxEnergy = globalEnergy;
		if (fabs(globalEnergy) < minAbsEnergy) { // note down the current charge of each element when the minimum energy was reached
			minAbsEnergy = fabs(globalEnergy);
			for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) minAbsEnergyValue[elemId] = currValue[elemId];
		}
		//printf("\nE=%+.3lf\t", globalEnergy);
		//for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) if (stepValue[elemId] != 0.0f) printf("[%c]=%+.3f ; ", ElementsList[elemId].letter, (signValue[elemId] * currValue[elemId]));
	}
	free(chargedAtomsArray);
	time_end = time(NULL);
	time_elapsed = difftime(time_end, time_start);
	printf(" (%.0lfm %.0lfs)\n", (time_elapsed / 60.0), fmod(time_elapsed, 60.0));
	printf("\tSystem energy   : min = %+.3lf ; max = %+.3lf ; minAbsolute = %+.3lf\n", minEnergy, maxEnergy, minAbsEnergy);
	printf("\tElement charges : ");
	for (elemId = 0; elemId <= NUM_ELEMENTS; elemId++) {
		if (stepValue[elemId] == 0.0f) continue;
		printf("[%c]=%+.3f ; ", ElementsList[elemId].letter, (signValue[elemId] * minAbsEnergyValue[elemId]));
	}
	printf("\n");
	/**/
}
