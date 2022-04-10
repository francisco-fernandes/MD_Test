#include <limits.h>

#define NUM_ELEMENTS 6
#define NUM_MOLECULES 23

#define MAX_NUM_ATOMS_PER_AA 24
#define MAX_NUM_BONDS_PER_ATOM 4

typedef struct Element {
	char letter;
	char name[11];
	char color[7];
	char colorHexCode[8];
	char colorHexCodeDarker[8];
	char colorHexCodeLighter[8];
	unsigned char colorHSLCode[3];
	char maxNumBonds;
	float atomicRadius;
} Element;

typedef struct Molecule {
	char letter;
	char abbrev[4];
	char name[14];
	unsigned char numberOfAtoms;
	unsigned char numberOfBonds; // number of distinct bonds among all molecule atoms (including to previous and next residues' atoms)
	char atomsCodes[MAX_NUM_ATOMS_PER_AA][5];
	signed char atomsBonds[MAX_NUM_ATOMS_PER_AA][MAX_NUM_BONDS_PER_ATOM]; // contains ids from (0) to (numberOfAtoms-1)
} Molecule;

static const Element ElementsList[(NUM_ELEMENTS + 1)] = {
		{ '\0', "", "", "", "", "", {0,0,0}, 0 , 0.0 },														// 0
		{ 'H', "Hydrogen"	, "white"	, "#DFDFDF"	, "#000000"	, "#FFFFFF"	, {  0,  0, 80}	, 1	, 0.37f },	// 1
		{ 'C', "Carbon"		, "black"	, "#404040"	, "#000000"	, "#C0C0C0"	, {  0,  0, 20}	, 4	, 0.77f },	// 2
		{ 'O', "Oxygen"		, "red"		, "#FF0000"	, "#7F0000"	, "#FFC080"	, {  0,100, 50}	, 6	, 0.66f },	// 3
		{ 'N', "Nitrogen"	, "blue"	, "#0000FF"	, "#00007F"	, "#C0C0FF"	, {240,100, 50}	, 5	, 0.70f },	// 4
		{ 'S', "Sulfur"		, "yellow"	, "#FFFF00"	, "#7F7F00"	, "#FFFFC0"	, { 60,100, 50}	, 6	, 1.04f },	// 5
		{ 'P', "Phosphorus"	, "orange"	, "#FF7F00"	, "#7F3F00"	, "#FFE0C0"	, { 30,100, 50}	, 5	, 1.10f }	// 6
};

static const unsigned char ElementIdLookupTable[32] =
//	 _,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,_,_,_,_,_
	{0,0,0,2,0,0,0,0,1,0,0,0,0,0,4,3,6,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0};

// [0] ***           : (A)**, (C)ys, (G)**, (H)**, (I)le, (L)**, (M)et, (P)**, (S)er, (T)**, (V)al
// [1] A**, L**, P** : al(A), ph(E), ar(G), as(N), pr(O), as(P), ly(S), le(U)
// [2] G**, H**, T** : ho(H), gl(N), tr(P), t*(R), hi(S), gl(U), gl(Y)
// [3] T*R           : t(H)r, t(Y)r
static const signed char MoleculeIdLookupTable[4][32] = {
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//    A*     Cys              G*  H* Ile          L* Met          P*         Ser  T*     Val
{ 00, -1, 00,  5, 00, 00, 00, -2, -2, 10, 00, 00, -1, 13, 00, 00, -1, 00, 00, 16, -2, 00, 20, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // [0]
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//   Ala             Phe     Arg                         Asn Pro Asp         Lys     Leu
{ 00,  1, 00, 00, 00, 14, 00,  2, 00, 00, 00, 00, 00, 00,  3, 15,  4, 00, 00, 12, 00, 11, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // [1]
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//                               HOH                     Gln     Trp     T*r His     Glu             Gly
{ 00, 00, 00, 00, 00, 00, 00, 00, 21, 00, 00, 00, 00, 00,  6, 00, 18, 00, -3,  9, 00,  7, 00, 00, 00,  8, 00, 00, 00, 00, 00, 00 }, // [2]
//_ , A , B , C , D , E , F , G , H , I , J , K , L , M , N , O , P , Q , R , S , T , U , V , W , X , Y , Z , _ , _ , _ , _ , _
//                               Thr                                                                 Tyr
{ 00, 00, 00, 00, 00, 00, 00, 00, 17, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 19, 00, 00, 00, 00, 00, 00 }  // [3]
};

// TODO: Check if need to add 'N' (from "HN", "HN1" and "HN2")
static const unsigned char AtomPlacementLookupTable[32] = // Alpha,Beta,Gamma,Delta,Epsilon,Zeta,(H)Eta ; 'X' (from "HXT" and "OXT") same as 'A'
//	 _,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,_,_,_,_,_
	{0,1,2,0,4,5,0,3,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,6,0,0,0,0,0};


// TODO: Lookup "atom id" by "atom code" and "molecule id"
// All atoms codes:
// N , NE , NE(1,2) , NH(1,2) , ND(1,2) , NZ
// O , OXT , OD(1,2) , OE(1,2) , OG , OG1 , OH
// C , CA , CB , CG , CG(1,2) , CD , CD(1,2) , CE, CE(1,2,3) , CZ , CZ(2,3) , CH2
// H , HXT , HA , HA(1,2) , HB , HB(1,2,3) , HG , HG(1,2) , HG1 , HG1(1,2,3) , HG2(1,2,3) ,
//		HD(1,2) , HD1(1,2,3) , HD2(1,2,3) , HE(1,2,3) , HE2(1,2) ,
//		HH , HH1(1,2) , HH2, HH2(1,2) , HZ , HZ(1,2,3)
// SG , SD


// TODO: add atom type by their connection, e.g. carbon with 3 Hs, carbon with double bond, carbon with oxygen, etc
// TODO: find way to set bonds for extra "HXT" and "OXT" atoms
// TODO: add pseudo-nucleotides for DNA and RNA
//  00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23
// ___,Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val,HOH,HXT,OXT
static const Molecule MoleculesList[(NUM_MOLECULES+1)] = {
/*00*/	{ '\0', "", "", 0, 0, {0}, {0} },
/*01*/	{ 'A', "Ala", "Alanine"			, 10,	11, // (-CH3)
		//  { 0 ,  1 , 2 , 3 ,  4 , 5 ,  6 ,   7 ,   8 ,   9 }
			{"N","CA","C","O","CB","H","HA","HB1","HB2","HB3"},
			{	{-1,1,5,SCHAR_MAX},					//"N"
				{0,2,4,6},							//"CA"
				{1,3,+10,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,7,8,9},							//"CB"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB1"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HB3"
		},
/*02*/	{ 'R', "Arg", "Arginine"		, 24,	25, // (-CH2-CH2-CH2-NH-C(-NH2)(-NH2))
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,  7 ,  8 ,   9 ,  10 , 11, 12 ,  13 ,  14 ,  15 ,  16 ,  17 ,  18 , 19 ,   20 ,   21 ,   22 ,   23 }
			{"N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2","H","HA","HB2","HB3","HG2","HG3","HD2","HD3","HE","HH11","HH12","HH21","HH22"},
			{	{-1,1,11,SCHAR_MAX},				//"N"
				{0,2,4,12},							//"CA"
				{1,3,+24,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,13,14},						//"CB"
				{4,6,15,16},						//"CG"
				{5,7,17,18},						//"CD"
				{6,8,19,SCHAR_MAX},					//"NE"
				{7,9,10,SCHAR_MAX},					//"CZ"
				{8,20,21,SCHAR_MAX},				//"NH1"
				{8,22,23,SCHAR_MAX},				//"NH2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD3"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE"
				{9,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HH11"
				{9,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HH12"
				{10,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HH21"
				{10,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }//"HH22"
		},
/*03*/	{ 'N', "Asn", "Asparagine"		, 14,	15, // (-CH2-CO-NH2)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 , 8 ,  9 ,  10 ,  11 ,   12 ,   13 }
			{"N","CA","C","O","CB","CG","OD1","ND2","H","HA","HB2","HB3","HD21","HD22"},
			{	{-1,1,8,SCHAR_MAX},					//"N"
				{0,2,4,9},							//"CA"
				{1,3,+14,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,10,11},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OD1"
				{5,12,13,SCHAR_MAX},				//"ND2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD21"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HD22"
		},
/*04*/	{ 'D', "Asp", "Aspartic acid"	, 12,	13, // (-CH2-CO2)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 , 8 ,  9 ,  10 ,  11 }
			{"N","CA","C","O","CB","CG","OD1","OD2","H","HA","HB2","HB3"},
			{	{-1,1,8,SCHAR_MAX},					//"N"
				{0,2,4,9},							//"CA"
				{1,3,+12,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,10,11},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OD1"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OD2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HB3"
		},								
/*05*/	{ 'C', "Cys", "Cysteine"		, 11,	12, // (-CH2-SH) ... NOTE: HG can be absent
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 , 6 ,  7 ,   8 ,   9 , 10 }
			{"N","CA","C","O","CB","SG","H","HA","HB2","HB3","HG"},
			{	{-1,1,6,SCHAR_MAX},					//"N"
				{0,2,4,7},							//"CA"
				{1,3,+10,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,8,9},							//"CB"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"SG"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HG"
		},
/*06*/	{ 'Q', "Gln", "Glutamine"		, 17,	18, // (-CH2-CH2-CO-N2)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,   7 ,   8 , 9 , 10 ,  11 ,  12 ,  13 ,  14 ,   15 ,   16 }
			{"N","CA","C","O","CB","CG","CD","OE1","NE2","H","HA","HB2","HB3","HG2","HG3","HE21","HE22"},
			{	{-1,1,9,SCHAR_MAX},					//"N"
				{0,2,4,10},							//"CA"
				{1,3,+17,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,11,12},						//"CB"
				{4,6,13,14},						//"CG"
				{5,7,8,SCHAR_MAX},					//"CD"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OE1"
				{6,15,16,SCHAR_MAX},				//"NE2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG3"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE21"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HE22"
		},						
/*07*/	{ 'E', "Glu", "Glutamic acid"	, 15,	16, // (-CH2-CH2-CO2)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,   7 ,   8 , 9 , 10 ,  11 ,  12 ,  13 ,  14 }
			{"N","CA","C","O","CB","CG","CD","OE1","OE2","H","HA","HB2","HB3","HG2","HG3"},
			{	{-1,1,9,SCHAR_MAX},					//"N"
				{0,2,4,10},							//"CA"
				{1,3,+15,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,11,12},						//"CB"
				{4,6,13,14},						//"CG"
				{5,7,8,SCHAR_MAX},					//"CD"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OE1"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"OE2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HG3"	
		},
/*08*/	{ 'G', "Gly", "Glycine"			,	7,	8, // (-H)
		//  { 0 ,  1 , 2 , 3 , 4 ,   5 ,   6 }
			{"N","CA","C","O","H","HA2","HA3"},
			{	{-1,1,4,SCHAR_MAX},					//"N"
				{0,2,5,6},							//"CA"
				{1,3,+7,SCHAR_MAX},					//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA2"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HA3"
		},
/*09*/	{ 'H', "His", "Histidine"		, 18,	19, // (-CH2-C[-CH-N-CH-NH-]) ... NOTE: can have both HD1 and HE2, or only one of them
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 ,   8 ,   9 , 10, 11 ,  12 ,  13 ,  14 ,  15 ,  16 ,  17 }
			{"N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2","H","HA","HB2","HB3","HD1","HD2","HE1","HE2"},
			{	{-1,1,10,SCHAR_MAX},				//"N"
				{0,2,4,11},							//"CA"
				{1,3,+18,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,12,13},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,8,14,SCHAR_MAX},					//"ND1"
				{5,9,15,SCHAR_MAX},					//"CD2"
				{6,9,16,SCHAR_MAX},					//"CE1"
				{7,8,17,SCHAR_MAX},					//"NE2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD1"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE1"
				{9,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HE2"
		},
/*10*/	{ 'I', "Ile", "Isoleucine"		, 19,	20, // (-CH(-CH3)(-CH2-CH3))
		//  { 0 ,  1 , 2 , 3 ,  4 ,   5 ,   6 ,   7 , 8 ,  9 , 10 ,   11 ,   12 ,   13 ,   14 ,   15 ,   16 ,   17 ,   18 }
			{"N","CA","C","O","CB","CG1","CG2","CD1","H","HA","HB","HG12","HG13","HG21","HG22","HG23","HD11","HD12","HD13"},
			{	{-1,1,8,SCHAR_MAX},					//"N"
				{0,2,4,9},							//"CA"
				{1,3,+19,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,6,10},							//"CB"
				{4,7,11,12},						//"CG1"
				{4,13,14,15},						//"CG2"
				{5,16,17,18},						//"CD1"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG12"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG13"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG21"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG22"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG23"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD11"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD12"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HD13"
		},
/*11*/	{ 'L', "Leu", "Leucine"			, 19,	20, // (-CH2-CH(-CH3)(-CH3))
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 , 8 ,  9 ,  10 ,  11 , 12 ,   13 ,   14 ,   15 ,   16 ,   17 ,   18 }
			{"N","CA","C","O","CB","CG","CD1","CD2","H","HA","HB2","HB3","HG","HD11","HD12","HD13","HD21","HD22","HD23"},
			{	{-1,1,8,SCHAR_MAX},					//"N"
				{0,2,4,9},							//"CA"
				{1,3,+19,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,10,11},						//"CB"
				{4,6,7,12},							//"CG"
				{5,13,14,15},						//"CD1"
				{5,16,17,18},						//"CD2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD11"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD12"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD13"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD21"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD22"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HD23"
		},
/*12*/	{ 'K', "Lys", "Lysine"			, 22,	22, // (-CH2-CH2-CH2-CH2-NH3) ... NOTE: HZ3 can be absent
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,  7 ,  8 , 9 , 10 ,  11 ,  12 ,  13 ,  14 ,  15 ,  16 ,  17 ,  18 ,  19 ,  20 ,  21 }
			{"N","CA","C","O","CB","CG","CD","CE","NZ","H","HA","HB2","HB3","HG2","HG3","HD2","HD3","HE2","HE3","HZ1","HZ2","HZ3"},
			{	{-1,1,9,SCHAR_MAX},					//"N"
				{0,2,4,10},							//"CA"
				{1,3,+22,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,11,12},						//"CB"
				{4,6,13,14},						//"CG"
				{5,7,15,16},						//"CD"
				{6,8,17,18},						//"CE"
				{7,19,20,21},						//"NZ"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD3"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE2"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE3"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HZ1"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HZ2"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HZ3"
		},
/*13*/	{ 'M', "Met", "Methionine"		, 17,	18, // (-CH2-CH2-S-CH3)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,  7 , 8 ,  9 ,  10 ,  11 ,  12 ,  13 ,  14 ,  15 ,  16 }
			{"N","CA","C","O","CB","CG","SD","CE","H","HA","HB2","HB3","HG2","HG3","HE1","HE2","HE3"},
			{	{-1,1,8,SCHAR_MAX},					//"N"
				{0,2,4,9},							//"CA"
				{1,3,+17,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,10,11},						//"CB"
				{4,6,12,13},						//"CG"
				{5,7,SCHAR_MAX,SCHAR_MAX},			//"SD"
				{6,14,15,16},						//"CE"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG3"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE1"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE2"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HE3"
		},
/*14*/	{ 'F', "Phe", "Phenylalanine"	, 20,	22, // (-CH2-C[-CH-CH-CH-CH-CH-])
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 ,   8 ,   9 , 10 , 11, 12 ,  13 ,  14 ,  15 ,  16 ,  17 ,  18 , 19 }
			{"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","H","HA","HB2","HB3","HD1","HD2","HE1","HE2","HZ"},
			{	{-1,1,11,SCHAR_MAX},				//"N"
				{0,2,4,12},							//"CA"
				{1,3,+20,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,13,14},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,8,15,SCHAR_MAX},					//"CD1"
				{5,9,16,SCHAR_MAX},					//"CD2"
				{6,10,17,SCHAR_MAX},				//"CE1"
				{7,10,18,SCHAR_MAX},				//"CE2"
				{8,9,19,SCHAR_MAX},					//"CZ"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD1"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE1"
				{9,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE2"
				{10,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }//"HZ"
		},
/*15*/	{ 'P', "Pro", "Proline"			, 14,	16, // |C|[-CH2-CH2-CH2-|N|-]
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,  6 ,  7 ,   8 ,   9 ,  10 ,  11 ,  12 ,  13 }
			{"N","CA","C","O","CB","CG","CD","HA","HB2","HB3","HG2","HG3","HD2","HD3"},
			{	{-1,1,6,SCHAR_MAX},					//"N"
				{0,2,4,7},							//"CA"
				{1,3,+14,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,8,9},							//"CB"
				{4,6,10,11},						//"CG"
				{0,5,12,13},						//"CD"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG2"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HD3"
		},
/*16*/	{ 'S', "Ser", "Serine"			, 11,	12, // (-CH2-OH)
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 , 6 ,  7 ,   8 ,   9 , 10 }
			{"N","CA","C","O","CB","OG","H","HA","HB2","HB3","HG"},
			{	{-1,1,6,SCHAR_MAX},					//"N"
				{0,2,4,7},							//"CA"
				{1,3,+11,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,8,9},							//"CB"
				{4,10,SCHAR_MAX,SCHAR_MAX},			//"OG"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HG"
		},
/*17*/	{ 'T', "Thr", "Threonine"		, 14,	15, // (-CH(-OH)(-CH3))
		//  { 0 ,  1 , 2 , 3 ,  4 ,   5 ,   6 , 7 ,  8 ,  9 ,  10 ,   11 ,   12 ,   13 }
			{"N","CA","C","O","CB","OG1","CG2","H","HA","HB","HG1","HG21","HG22","HG23"},
			{	{-1,1,7,SCHAR_MAX},					//"N"
				{0,2,4,8},							//"CA"
				{1,3,+14,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,6,9},							//"CB"
				{4,10,SCHAR_MAX,SCHAR_MAX},			//"OG1"
				{4,11,12,13},						//"CG2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG1"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG21"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG22"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HG23"
		},
/*18*/	{ 'W', "Trp", "Tryptophan"		, 24,	27, // (-CH2-C[-C[-CH-CH-CH-CH]-C-NH-CH-])
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 ,   8 ,   9 ,  10 ,  11 ,  12 ,  13 , 14, 15 ,  16 ,  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 }
			{"N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2","H","HA","HB2","HB3","HD1","HE1","HE3","HZ2","HZ3","HH2"},
			{	{-1,1,14,SCHAR_MAX},				//"N"
				{0,2,4,15},							//"CA"
				{1,3,+24,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,16,17},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,8,18,SCHAR_MAX},					//"CD1"
				{5,9,10,SCHAR_MAX},					//"CD2"
				{6,9,19,SCHAR_MAX},					//"NE1"
				{7,8,11,SCHAR_MAX},					//"CE2"
				{7,12,20,SCHAR_MAX},				//"CE3"
				{9,13,21,SCHAR_MAX},				//"CZ2"
				{10,13,22,SCHAR_MAX},				//"CZ3"
				{11,12,23,SCHAR_MAX},				//"CH2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD1"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE1"
				{10,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE3"
				{11,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HZ2"
				{12,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HZ3"
				{13,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }//"HH2"
		},
/*19*/	{ 'Y', "Tyr", "Tyrosine"		, 21,	23, // (-CH2-C[-CH-CH-COH-CH-CH-])
		//  { 0 ,  1 , 2 , 3 ,  4 ,  5 ,   6 ,   7 ,   8 ,   9 , 10 , 11 , 12, 13 ,  14 ,  15 ,  16 ,  17 ,  18 ,  19 , 20 }
			{"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH","H","HA","HB2","HB3","HD1","HD2","HE1","HE2","HH"},
			{	{-1,1,12,SCHAR_MAX},				//"N"
				{0,2,4,13},							//"CA"
				{1,3,+21,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,14,15},						//"CB"
				{4,6,7,SCHAR_MAX},					//"CG"
				{5,8,16,SCHAR_MAX},					//"CD1"
				{5,9,17,SCHAR_MAX},					//"CD2"
				{6,10,18,SCHAR_MAX},				//"CE1"
				{7,10,19,SCHAR_MAX},				//"CE2"
				{8,9,11,SCHAR_MAX},					//"CZ"
				{10,20,SCHAR_MAX,SCHAR_MAX},		//"OH"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB2"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB3"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD1"
				{7,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HD2"
				{8,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE1"
				{9,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HE2"
				{11,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }//"HH"
		},
/*20*/	{ 'V', "Val", "Valine"			, 16,	17, // (-CH(-CH3)(-CH3))
		//  { 0 ,  1 , 2 , 3 ,  4 ,   5 ,   6 , 7 ,  8 ,  9 ,   10 ,   11 ,   12 ,   13 ,   14 ,   15 }
			{"N","CA","C","O","CB","CG1","CG2","H","HA","HB","HG11","HG12","HG13","HG21","HG22","HG23"},
			{	{-1,1,7,SCHAR_MAX},					//"N"
				{0,2,4,8},							//"CA"
				{1,3,+16,SCHAR_MAX},				//"C"
				{2,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"O"
				{1,5,6,9},							//"CB"
				{4,10,11,12},						//"CG1"
				{4,13,14,15},						//"CG2"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H"
				{1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HA"
				{4,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HB"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG11"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG12"
				{5,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG13"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG21"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"HG22"
				{6,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"HG23"
		},
/*21*/	{ 'w', "HOH", "Water"			,	3,	2,	// (O(-H)(-H))
		//  { 0 ,  1 ,  2 }
			{"O","H1","H2"},
			{	{1,2,SCHAR_MAX,SCHAR_MAX},			//"O"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX},	//"H1"
				{0,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }	//"H2"
		},
/*22*/	{ 'h', "HXT", "Term.hydrogen"	,	1,	1,	// (H)
		//  {  0  }
			{"HXT"},
			{	{-1,SCHAR_MAX,SCHAR_MAX,SCHAR_MAX} }//"HXT" ? should connect to "OXT" (-1) or to "N" (not necessarily -1) ?
		},
/*23*/	{ 'o', "OXT", "Term.oxygen"		,	1,	2,	// (O)
		//  {  0  }
			{"OXT"},
			{	{-1,+1,SCHAR_MAX,SCHAR_MAX} }		//"OXT" ? should connect to "C" (-1) and to "HXT" (+1) ?
		}
};


// 0,...,31
#define CHAR2ID(c) ((int)((c) & 0x1F))
// H,C,O,N,S,P
#define GET_ELEMENT_LETTER(moleculeId,atomId) (MoleculesList[(moleculeId)].atomsCodes[(atomId)][0])
// A,B,G,D,E,Z,H
#define GET_PLACEMENT_ORDER(moleculeId,atomId) (AtomPlacementLookupTable[ CHAR2ID( MoleculesList[(moleculeId)].atomsCodes[(atomId)][1] ) ])
// 0,1,2,3
#define GET_EXTRA_PLACEMENT_ORDER(moleculeId,atomId) ((int) ( (MoleculesList[(moleculeId)].atomsCodes[(atomId)][2]) & 0x0F ) )

int GetElementCode(const char *atomName);
int GetMoleculeCode(char *moleculeName);
