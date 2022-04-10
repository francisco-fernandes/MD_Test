
typedef struct Atom {
	char letter;
	char name[9];
	char color[7];
} Atom;

const Atom AtomsList[6] = {
		{ '\0', "", "" },				// 0
		{ 'H', "Hydrogen", "white" },	// 1
		{ 'C', "Carbon"  , "black" },	// 2
		{ 'O', "Oxygen"  , "red" },		// 3
		{ 'N', "Nitrogen", "blue" },	// 4
		{ 'S', "Sulfur"  , "yellow" }	// 5
};


#define CHAR2ID(c) ((int)((c) & 0x3F))

const unsigned char AtomIdLookupTable[32] =
//	 _,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,_,_,_,_,_
	{0,0,0,2,0,0,0,0,1,0,0,0,0,0,4,3,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0};

const signed char MoleculeIdLookupTable[4][32] = {
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//    A*     Cys              G*  His Ile         L*  Met         P*         Ser  T*      Val
{ 00, -1, 00, 05, 00, 00, 00, -2,  9, 10, 00, 00, -1, 13, 00, 00, -1, 00, 00, 16, -2, 00, 20, 00, 00, 00, 00, 00, 00, 00, 00, 00 },
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//   Ala             Phe     Arg                         Asn Pro Asp         Lys     Leu
{ 00, 01, 00, 00, 00, 14, 00, 02, 00, 00, 00, 00, 00, 00, 03, 15, 04, 00, 00, 12, 00, 11, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 },
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//                                                       Gln     Trp     T*r         Glu             Gly
{ 00, 00, 00, 00, 00, 00, 00, 06, 00, 00, 00, 00, 00, 00, 06, 00, 18, 00, -3, 00, 00, 07, 00, 00, 00,  8, 00, 00, 00, 00, 00, 00 },
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//                               Thr                                                                 Tyr
{ 00, 00, 00, 00, 00, 00, 00, 00, 17, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 19, 00, 00, 00, 00, 00, 00 }
};

typedef struct Molecule {
	char letter;
	char abbrev[4];
	char name[14];
	unsigned char numberOfAtoms;
	unsigned char atomsCountByAlphabetPosition[8]; // None,(A)lpha,(B)eta,(G)amma,(D)elta,(E)psilon,(Z)eta,(H)Eta
	char atomsLetters[25];
	signed char *atomsBonds[4]; // atomsBonds[][4]
} Molecule;

const Molecule AminoacidsList[21] = {
/*00*/		{ '\0', "", "", 0, NULL, NULL, NULL },
			//	  { 0, A, B, G, D, E, Z, H }
/*01*/		{ 'A', "Ala", "Alanine"				//"0123456789"					 "NCCOCHHHHH" {N,CA,C,O,CB,H,HA,HB1,HB2,HB3} ; (-CH3)
			, 10, { 4, 2, 4, 0, 0, 0, 0, 0 }, "NCOHCHCHHH" },				// {N,C,O,H,CA,HA,CB,HB1,HB2,HB3}
/*02*/		{ 'R', "Arg", "Arginine"			//"012345678901234567890123"	 "NCCOCCCNCNNHHHHHHHHHHHHH" {N,CA,C,O,CB,CG,CD,NE,CZ,NH1,NH2,H,HA,HB2,HB3,HG2,HG3,HD2,HD3,HE,HH11,HH12,HH21,HH22} ; (-CH2-CH2-CH2-NH-C(-NH2)(-NH2))
			, 24, { 4, 2, 3, 3, 3, 2, 1, 6 }, "NCOHCHCHHCHHCHHNHCNNHHHH" },	// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG2,HG3,CD,HD2,HD3,NE,HE,CZ,NH1,NH2,HH11,HH12,HH21,HH22}
/*03*/		{ 'N', "Asn", "Asparagine"			//"01234567890123"				 "NCCOCCONHHHHHH" {N,CA,C,O,CB,CG,OD1,ND2,H,HA,HB2,HB3,HD21,HD22} ; (-CH2-CO-NH2)
			, 14, { 4, 2, 3, 1, 4, 0, 0, 0 }, "NCOHCHCHHCONHH" },			// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,OD1,ND2,HD21,HD22}
/*04*/		{ 'D', "Asp", "Aspartic acid"		//"012345678901"				 "NCCOCCOOHHHH" {N,CA,C,O,CB,CG,OD1,OD2,H,HA,HB2,HB3} ; (-CH2-CO2)
			, 12, { 4, 2, 3, 1, 2, 0, 0, 0 }, "NCOHCHCHHCOO" },				// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,OD1,OD2}
/*05*/		{ 'C', "Cys", "Cysteine"			//"0123456789"					 "NCCOCSHHHH" {N,CA,C,O,CB,SG,H,HA,HB2,HB3} ; (-CH2-SH)
			, 10, { 4, 2, 3, 1, 0, 0, 0, 0 }, "NCOHCHCHHS" },				// {N,C,O,H,CA,HA,CB,HB2,HB3,SG} ... ?HG?
/*06*/		{ 'Q', "Gln", "Glutamine"			//"01234567890123456"			 "NCCOCCCONHHHHHHHH" {N,CA,C,O,CB,CG,CD,OE1,NE2,H,HA,HB2,HB3,HG2,HG3,HE21,HE22} ; (-CH2-CH2-CO-N2)
			, 17, { 4, 2, 3, 3, 1, 4, 0, 0 }, "NCOHCHCHHCHHCONHH" },		// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG2,HG3,CD,OE1,NE2,HE21,HE22}
/*07*/		{ 'E', "Glu", "Glutamic acid"		//"012345678901234"				 "NCCOCCCOOHHHHHH" {N,CA,C,O,CB,CG,CD,OE1,OE2,H,HA,HB2,HB3,HG2,HG3} ; (-CH2-CH2-CO2)
			, 15, { 4, 2, 3, 3, 1, 2, 0, 0 }, "NCOHCHCHHCHHCOO" },			// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG2,HG3,CD,OE1,OE2}
/*08*/		{ 'G', "Gly", "Glycine"				//"0123456"						 "NCCOHHH" {N,CA,C,O,H,HA2,HA3} ; (-H)
			, 07, { 4, 3, 0, 0, 0, 0, 0, 0 }, "NCOHCHH" },					// {N,C,O,H,CA,HA2,HA3}
/*09*/		{ 'H', "His", "Histidine"			//"012345678901234567"			 "NCCOCCNCCNHHHHHHHH" {N,CA,C,O,CB,CG,ND1,CD2,CE1,NE2,H,HA,HB2,HB3,HD1,HD2,HE1,HE2} ; (-CH2-C[-CH-N-CH-NH-])
			, 18, { 4, 2, 3, 1, 4, 4, 0, 0 }, "NCOHCHCHHCNCHHCNHH" },		// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,ND1,CD2,HD1,HD2,CE1,NE2,HE1,HE2}
/*10*/		{ 'I', "Ile", "Isoleucine"			//"0123456789012345678"			 "NCCOCCCCHHHHHHHHHHH" {N,CA,C,O,CB,CG1,CG2,CD1,H,HA,HB,HG12,HG13,HG21,HG22,HG23,HD11,HD12,HD13} ; (-CH(-CH3)(-CH2-CH3))
			, 19, { 4, 2, 2, 7, 4, 0, 0, 0 }, "NCOHCHCHCCHHHHHCHHH" },		// {N,C,O,H,CA,HA,CB,HB,CG1,CG2,HG12,HG13,HG21,HG22,HG23,CD1,HD11,HD12,HD13}
/*11*/		{ 'L', "Leu", "Leucine"				//"0123456789012345678"			 "NCCOCCCCHHHHHHHHHHH" {N,CA,C,O,CB,CG,CD1,CD2,H,HA,HB2,HB3,HG,HD11,HD12,HD13,HD21,HD22,HD23} ; (-CH2-CH(-CH3)(-CH3))
			, 19, { 4, 2, 3, 2, 8, 0, 0, 0 }, "NCOHCHCHHCHCCHHHHHH" },		// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG,CD1,CD2,HD11,HD12,HD13,HD21,HD22,HD23}
/*12*/		{ 'K', "Lys", "Lysine"				//"0123456789012345678901"		 "NCCOCCCCNHHHHHHHHHHHHH" {N,CA,C,O,CB,CG,CD,CE,NZ,H,HA,HB2,HB3,HG2,HG3,HD2,HD3,HE2,HE3,HZ1,HZ2,HZ3} ; (-CH2-CH2-CH2-CH2-NH3)
			, 22, { 4, 2, 3, 3, 3, 3, 4, 0 }, "NCOHCHCHHCHHCHHCHHNHHH" },	// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG2,HG3,CD,HD2,HD3,CE,HE2,HE3,NZ,HZ1,HZ2,HZ3}
/*13*/		{ 'M', "Met", "Methionine"			//"01234567890123456"			 "NCCOCCSCHHHHHHHHH" {N,CA,C,O,CB,CG,SD,CE,H,HA,HB2,HB3,HG2,HG3,HE1,HE2,HE3} ; (-CH2-CH2-S-CH3)
			, 17, { 4, 2, 3, 3, 1, 4, 0, 0 }, "NCOHCHCHHCHHSCHHH" },		// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,HG2,HG3,SD,CE,HE1,HE2,HE3}
/*14*/		{ 'F', "Phe", "Phenylalanine"		//"01234567890123456789"		 "NCCOCCCCCCCHHHHHHHHH" {N,CA,C,O,CB,CG,CD1,CD2,CE1,CE2,CZ,H,HA,HB2,HB3,HD1,HD2,HE1,HE2,HZ} ; (-CH2-C[-CH-CH-CH-CH-CH-])
			, 20, { 4, 2, 3, 1, 4, 4, 2, 0 }, "NCOHCHCHHCCCHHCCHHCH" },		// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,CD1,CD2,HD1,HD2,CE1,CE2,HE1,HE2,CZ,HZ}
/*15*/		{ 'P', "Pro", "Proline"				//"01234567890123"				 "NCCOCCCHHHHHHH" {N,CA,C,O,CB,CG,CD,HA,HB2,HB3,HG2,HG3,HD2,HD3} ; |C|[-CH2-CH2-CH2-|N|-]
			, 14, { 3, 2, 3, 3, 3, 0, 0, 0 }, "NCOCHCHHCHHCHH" },			// {N,C,O,CA,HA,CB,HB2,HB3,CG,HG2,HG3,CD,HD2,HD3}
/*16*/		{ 'S', "Ser", "Serine"				//"01234567890"					 "NCCOCOHHHHH" {N,CA,C,O,CB,OG,H,HA,HB2,HB3,HG} ; (-CH2-OH)
			, 11, { 4, 2, 3, 2, 0, 0, 0, 0 }, "NCOHCHCHHOH" },				// {N,C,O,H,CA,HA,CB,HB2,HB3,OG,HG}
/*17*/		{ 'T', "Thr", "Threonine"			//"01234567890123"				 "NCCOCOCHHHHHHH" {N,CA,C,O,CB,OG1,CG2,H,HA,HB,HG1,HG21,HG22,HG23} ; (-CH(-OH)(-CH3))
			, 14, { 4, 2, 2, 6, 0, 0, 0, 0 }, "NCOHCHCHOCHHHH" },			// {N,C,O,H,CA,HA,CB,HB,OG1,CG2,HG1,HG21,HG22,HG23}
/*18*/		{ 'W', "Trp", "Tryptophan"			//"012345678901234567890123"	 "NCCOCCCCNCCCCCHHHHHHHHHH" {N,CA,C,O,CB,CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2,H,HA,HB2,HB3,HD1,HE1,HE3,HZ2,HZ3,HH2} ; (-CH2-C[-C[-CH-CH-CH-CH]-C-NH-CH-])
			, 24, { 4, 2, 3, 1, 3, 5, 4, 2 }, "NCOHCHCHHCCCHNCCHHCCHHCH" },	// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,CD1,CD2,HD1,NE1,CE2,CE3,HE1,HE3,CZ2,CZ3,HZ2,HZ3,CH2,HH2}
/*19*/		{ 'Y', "Tyr", "Tyrosine"			//"012345678901234567890"		 "NCCOCCCCCCCOHHHHHHHHH" {N,CA,C,O,CB,CG,CD1,CD2,CE1,CE2,CZ,OH,H,HA,HB2,HB3,HD1,HD2,HE1,HE2,HH} ; (-CH2-C[-CH-CH-COH-CH-CH-])
			, 21, { 4, 2, 3, 1, 4, 4, 1, 2 }, "NCOHCHCHHCCCHHCCHHCOH" },	// {N,C,O,H,CA,HA,CB,HB2,HB3,CG,CD1,CD2,HD1,HD2,CE1,CE2,HE1,HE2,CZ,OH,HH}
/*20*/		{ 'V', "Val", "Valine"				//"0123456789012345"			 "NCCOCCCCHHHHHHHHHHH" {N,CA,C,O,CB,CG1,CG2,H,HA,HB,HG11,HG12,HG13,HG21,HG22,HG23} ; (-CH(-CH3)(-CH3))
			, 16, { 4, 2, 2, 8, 0, 0, 0, 0 }, "NCOHCHCHCCHHHHHH" }			// {N,C,O,H,CA,HA,CB,HB,CG1,CG2,HG11,HG12,HG13,HG21,HG22,HG23}
};
