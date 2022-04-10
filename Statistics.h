typedef struct ChargedAtom {
	unsigned int id;
	unsigned char elementId;
	float charge;
} ChargedAtom;

typedef struct ForceFieldStats {
	float globalEnergy;
	float energyX;
	float energyY;
	float energyZ;
} ForceFieldStats;


void GetBondLengthStats();
void GetInteratomicForceStats();

#define ARRAY_GROW_SIZE 64
