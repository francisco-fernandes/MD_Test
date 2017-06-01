#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "SVGGraphics.h"
#include "AtomsArray.h"
#include "ChemicalEntities.h"

static FILE *svgFile;

static const int svg_pixels_per_angstrom = 50;
static const int svg_text_size = 100;
static const int svg_top_border_length = 125;
static const int svg_bottom_border_length = 75;
static const int svg_horizontal_border_length = 75;

static int svg_width = 0;
static int svg_height = 0;

static float minX = 0.0f;
static float maxY = 0.0f;
static float maxZ = 0.0f;

static const int svg_atom_border_size = 2;
static const int svg_num_depth_levels = 30;
static float svg_angstroms_per_depth_level = 0.0f;
static int *svg_blur_size_at_depth_level = NULL;
static float *svg_opacity_at_depth_level = NULL;

static int svg_current_depth_level = (-1);

int InitializeGraphics(char * svgFilename) {
	int i;
	if ((svgFile = fopen(svgFilename, "w")) == NULL) {
		fprintf(stderr, "> ERROR: Cannot create SVG file\n");
		return (-1);
	}
	fprintf(svgFile,"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 %d %d\">\n", svg_width, svg_height);
	fprintf(svgFile, "<style type=\"text/css\">\n");
	for (i = 0; i < svg_num_depth_levels; i++) {
		fprintf(svgFile, ".d%d { stroke-width: %d; stroke-opacity: %.2f; }\n", i, (svg_atom_border_size + svg_blur_size_at_depth_level[i]), svg_opacity_at_depth_level[i]);
	}
	fprintf(svgFile, "</style>\n");
	fprintf(svgFile, "<defs>\n");
	for (i = 1; i <= NUM_ELEMENTS; i++) {
		fprintf(svgFile, "<radialGradient id='g%c' fx='25%%' fy='25%%' r='50%%' spreadMethod='pad'>\n", ElementsList[i].letter);
		fprintf(svgFile, "\t<stop offset='0%%' stop-color='%s' />\n", (char *)ElementsList[i].colorCode3);
		fprintf(svgFile, "\t<stop offset='100%%' stop-color='%s' />\n", (char *)ElementsList[i].colorCode1);
		fprintf(svgFile, "</radialGradient>\n");
	}
	for (i = 1; i <= NUM_ELEMENTS; i++) {
		//fprintf(svgFile, "<circle id='%c' cx='0' cy='0' r='%d' fill='%s' stroke='%s' />\n",
		//	ElementsList[i].letter, (int)(ElementsList[i].atomicRadius * svg_pixels_per_angstrom), (char *)ElementsList[i].colorCode1, (char *)ElementsList[i].colorCode2);
		fprintf(svgFile, "<circle id='%c' cx='0' cy='0' r='%d' stroke='%s' style='fill:url(#g%c);' />\n",
			ElementsList[i].letter, (int)(ElementsList[i].atomicRadius * svg_pixels_per_angstrom), (char *)ElementsList[i].colorCode2, ElementsList[i].letter);
	}
	for (i = 0; i < svg_num_depth_levels; i++) {
		fprintf(svgFile, "<path id='p%d' d='m0,0 l%d,0 z'/>\n", i, 5*(svg_num_depth_levels-i));
		//fprintf(svgFile, "<path id='p%d' d='m0,0 a50,50 0 0,0 0,0 z'/>\n", i);
	}
	for (i = 0; i < svg_num_depth_levels; i++) {
		fprintf(svgFile, "<animateMotion xlink:href='#d%d' dur='10s' repeatCount='indefinite'><mpath xlink:href='#p%d'/></animateMotion>\n", i, i);
	}
	fprintf(svgFile, "</defs>\n");
	fprintf(svgFile, "<text x=\"%d\" y=\"%d\" font-size=\"%d\" font-weight=\"bolder\" text-anchor=\"middle\" alignment-baseline=\"middle\">%s</text>\n",
		(svg_width / 2), (svg_top_border_length / 2), svg_text_size, "PDB SVG Test");
	fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"grey\" />\n",
		svg_horizontal_border_length, svg_top_border_length, (svg_width - svg_horizontal_border_length), svg_top_border_length);
	fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"grey\" />\n",
		svg_horizontal_border_length, svg_top_border_length, svg_horizontal_border_length, (svg_height - svg_bottom_border_length));
	fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"grey\" />\n",
		svg_horizontal_border_length, (svg_height - svg_bottom_border_length), (svg_width - svg_horizontal_border_length), (svg_height - svg_bottom_border_length));
	fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"grey\" />\n",
		(svg_width - svg_horizontal_border_length), svg_top_border_length, (svg_width - svg_horizontal_border_length), (svg_height - svg_bottom_border_length));
	//<g></g> : group
	return 1;
}

int FinalizeGraphics() {
	fprintf(svgFile, "</svg>\n");
	if (fclose(svgFile) == EOF) {
		fprintf(stderr, "> ERROR: Cannot write SVG file\n");
		return (-1);
	}
	return 1;
}

void DrawSingleAtom(unsigned int atomId) {
	Atom *atom;
	int x, y;
	unsigned char elemId;
	//int radius;
	//char *color;
	//int zfrontdist;
	int depthlevel;
	atom = &(allAtomsArray[atomId]);
	x = svg_horizontal_border_length + (int)(((atom->coordX) - minX)*svg_pixels_per_angstrom);
	y = svg_top_border_length + (int)((maxY - (atom->coordY))*svg_pixels_per_angstrom);
	elemId = GetAtomElementCode(atom);
	//color = (char *)ElementsList[elemId].color;
	//radius = (int)(ElementsList[elemId].atomicRadius * svg_pixels_per_angstrom);
	//zfrontdist = 1 + (int)((maxZ - (atom->coordZ)) / svg_angstroms_per_depth_level); // [ 1 ... svg_atom_blur_max_size ]
	//fprintf(svgFile, "<circle cx='%d' cy='%d' r='%d' style='stroke:black;stroke-width:%d;fill:%s' />\n", x, y, radius, border, color);
	//fprintf(svgFile, "<use xlink:href='#%c' x='%d' y='%d' stroke-width='%d' stroke-opacity='%.2f' />\n",
	//	ElementsList[elemId].letter, x, y, (svg_atom_border_size + zfrontdist), 1.00f - ((float)zfrontdist / (float)svg_num_depth_levels));
	depthlevel = 0 + (int)((maxZ - (atom->coordZ)) / svg_angstroms_per_depth_level); // [ 0 ... (svg_atom_blur_max_size - 1) ]
	if (depthlevel != svg_current_depth_level) {
		if (svg_current_depth_level != (-1)) fprintf(svgFile, "</g>\n");
		fprintf(svgFile, "<g id='d%d'>\n", depthlevel);
		svg_current_depth_level = depthlevel;
	}
	fprintf(svgFile, "<use xlink:href='#%c' x='%d' y='%d' class='d%d' />\n", ElementsList[elemId].letter, x, y, depthlevel);

}

int SortAtomsCoordsByZValue(const void *a, const void *b) {
	float diff = (allAtomsArray[(*(unsigned int *)a)].coordZ - allAtomsArray[(*(unsigned int *)b)].coordZ);
	return ((diff > 0) ? (+1) : (-1));
}

void DrawAtoms() {
	unsigned int *atomsOrder;
	unsigned int numAtoms;
	unsigned int i;
	float minY, minZ;
	float maxX;
	minX = FLT_MAX;
	minY = FLT_MAX;
	minZ = FLT_MAX;
	maxX = FLT_MIN;
	maxY = FLT_MIN;
	maxZ = FLT_MIN;
	numAtoms = 1;
	while (numAtoms != totalNumAtoms) {
		if( (allResiduesArray[(allAtomsArray[numAtoms].fromResidueId)].fromModelId) != (allResiduesArray[(allAtomsArray[(numAtoms + 1)].fromResidueId)].fromModelId)) break;
		numAtoms++;
	}
	atomsOrder = (unsigned int *)calloc(numAtoms, sizeof(unsigned int));
	for (i = 1; i <= numAtoms; i++) {
		atomsOrder[(i-1)] = i;
		if (allAtomsArray[i].coordX < minX) minX = allAtomsArray[i].coordX;
		if (allAtomsArray[i].coordY < minY) minY = allAtomsArray[i].coordY;
		if (allAtomsArray[i].coordZ < minZ) minZ = allAtomsArray[i].coordZ;
		if (allAtomsArray[i].coordX > maxX) maxX = allAtomsArray[i].coordX;
		if (allAtomsArray[i].coordY > maxY) maxY = allAtomsArray[i].coordY;
		if (allAtomsArray[i].coordZ > maxZ) maxZ = allAtomsArray[i].coordZ;
	}
	qsort(atomsOrder, numAtoms, sizeof(unsigned int), SortAtomsCoordsByZValue);
	svg_width = (int)((maxX - minX)*svg_pixels_per_angstrom) + 2 * svg_horizontal_border_length;
	svg_height = (int)((maxY - minY)*svg_pixels_per_angstrom) + svg_top_border_length + svg_bottom_border_length;
	svg_angstroms_per_depth_level = (((maxZ - minZ) + 0.001f) / (float)(svg_num_depth_levels));
	svg_blur_size_at_depth_level = (int *)malloc((svg_num_depth_levels) * sizeof(int));
	svg_opacity_at_depth_level = (float *)malloc((svg_num_depth_levels) * sizeof(float));
	for (i = 0; i < (unsigned int)svg_num_depth_levels; i++) {
		//svg_blur_size_at_depth_level[i] = (i*i) / svg_num_depth_levels; // x^2/N : increases slowly and then speeds up, from 0 to N
		//svg_opacity_at_depth_level[i] = powf(((float)i - (float)svg_num_depth_levels) / (float)svg_num_depth_levels, 2.0f); // (x-N)^2/N : decreases quickly and then slows down, from N to 0 ; /N from 1.0 to 0.0
		svg_blur_size_at_depth_level[i] = i;
		svg_opacity_at_depth_level[i] = 1.0f - i*(1.0f / svg_num_depth_levels);
	}
	InitializeGraphics("test.svg");
	svg_current_depth_level = (-1);
	for (i = 0; i < numAtoms; i++) {
		DrawSingleAtom(atomsOrder[i]);
	}
	if (svg_current_depth_level != (-1)) fprintf(svgFile, "</g>\n");
	FinalizeGraphics();
	free(svg_opacity_at_depth_level);
	free(svg_blur_size_at_depth_level);
	free(atomsOrder);
}
