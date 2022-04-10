#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "SVGGraphics.h"
#include "AtomsArray.h"
#include "Statistics.h"
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
static unsigned int *atomsClosenessToSurface = NULL;

static float svg_shade_level_per_buried_order = 0.0f;
static float svg_fade_level_per_angstrom_distance = 0.0f;

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
	/*
	for (i = 1; i <= NUM_ELEMENTS; i++) {
		fprintf(svgFile, "<radialGradient id='g%c' fx='25%%' fy='25%%' r='50%%' spreadMethod='pad'>\n", ElementsList[i].letter);
		fprintf(svgFile, "\t<stop offset='0%%' stop-color='%s' />\n", (char *)ElementsList[i].colorHexCodeLighter);
		fprintf(svgFile, "\t<stop offset='100%%' stop-color='%s' />\n", (char *)ElementsList[i].colorHexCode);
		fprintf(svgFile, "</radialGradient>\n");
	}
	*/
	for (i = 1; i <= NUM_ELEMENTS; i++) {
		fprintf(svgFile, "<circle id='%c' cx='0' cy='0' r='%d' fill='inherit' stroke='%s' />\n",
			ElementsList[i].letter, (int)(ElementsList[i].atomicRadius * svg_pixels_per_angstrom), (char *)ElementsList[i].colorHexCodeDarker);
		/*
		fprintf(svgFile, "<circle id='%c' cx='0' cy='0' r='%d' stroke='%s' style='fill:url(#g%c);' />\n",
			ElementsList[i].letter, (int)(ElementsList[i].atomicRadius * svg_pixels_per_angstrom), (char *)ElementsList[i].colorHexCodeDarker, ElementsList[i].letter);
		*/
	}
	/*
	for (i = 0; i < svg_num_depth_levels; i++) {
		fprintf(svgFile, "<path id='p%d' d='m0,0 l%d,0 z'/>\n", i, 5*(svg_num_depth_levels-i));
		//fprintf(svgFile, "<path id='p%d' d='m0,0 a50,50 0 0,0 0,0 z'/>\n", i);
	}
	for (i = 0; i < svg_num_depth_levels; i++) {
		fprintf(svgFile, "<animateMotion xlink:href='#d%d' dur='10s' repeatCount='indefinite'><mpath xlink:href='#p%d'/></animateMotion>\n", i, i);
	}
	*/
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
	int depthlevel;
	int satval, lightval;
	atom = &(allAtomsArray[atomId]);
	x = svg_horizontal_border_length + (int)(((atom->coordX) - minX)*svg_pixels_per_angstrom);
	y = svg_top_border_length + (int)((maxY - (atom->coordY))*svg_pixels_per_angstrom);
	elemId = GetAtomElementCode(atom);
	depthlevel = (int)((maxZ - (atom->coordZ)) / svg_angstroms_per_depth_level); // [ 0 ... (svg_atom_blur_max_size - 1) ]
	if (depthlevel != svg_current_depth_level) {
		if (svg_current_depth_level != (-1)) fprintf(svgFile, "</g>\n");
		fprintf(svgFile, "<g id='d%d'>\n", depthlevel);
		svg_current_depth_level = depthlevel;
	}
	/*
	fprintf(svgFile, "<use xlink:href='#%c' x='%d' y='%d' class='d%d' fill='hsl(%d,%d%%,%d%%)' />\n",
		ElementsList[elemId].letter, x, y, depthlevel, ElementsList[elemId].colorHSLCode[0], ElementsList[elemId].colorHSLCode[1], ElementsList[elemId].colorHSLCode[2]);
	*/
	satval = (int)ElementsList[elemId].colorHSLCode[1];
	lightval = (int)ElementsList[elemId].colorHSLCode[2];
	if (satval != 0) {
		satval = (int)(((float)satval) - (svg_shade_level_per_buried_order)*((float)atomsClosenessToSurface[atomId])); // decrease color saturation value
		if (satval < 0) satval = 0;
	}
	else { // if the saturation level is already 0% (H and C), decrease the lightness level instead
		lightval = (int)(((float)lightval) - (svg_shade_level_per_buried_order/2)*((float)atomsClosenessToSurface[atomId]));
		if (lightval < 5) lightval = 5;
	}
	lightval = (int)(((float)lightval) + (maxZ - (atom->coordZ))*svg_fade_level_per_angstrom_distance); // increase color lightness value
	if (lightval > 100) lightval = 100;
	fprintf(svgFile, "<use xlink:href='#%c' x='%d' y='%d' class='d%d' fill='hsl(%d,%d%%,%d%%)' />\n",
		ElementsList[elemId].letter, x, y, depthlevel, ElementsList[elemId].colorHSLCode[0], satval, lightval);
}

int SortAtomsCoordsByXValue(const void *a, const void *b) {
	float diff = (allAtomsArray[(*(unsigned int *)a)].coordX - allAtomsArray[(*(unsigned int *)b)].coordX);
	return ((diff > 0) ? (+1) : (-1));
}

int SortAtomsCoordsByYValue(const void *a, const void *b) {
	float diff = (allAtomsArray[(*(unsigned int *)a)].coordY - allAtomsArray[(*(unsigned int *)b)].coordY);
	return ((diff > 0) ? (+1) : (-1));
}

int SortAtomsCoordsByZValue(const void *a, const void *b) {
	float diff = (allAtomsArray[(*(unsigned int *)a)].coordZ - allAtomsArray[(*(unsigned int *)b)].coordZ);
	return ((diff > 0) ? (+1) : (-1));
}

void DrawAtoms() {
	unsigned int *atomsOrder;
	unsigned int numAtoms;
	unsigned int i, k;
	float minY, minZ;
	float maxX;
	minX = FLT_MAX;
	minY = FLT_MAX;
	minZ = FLT_MAX;
	maxX = -(FLT_MAX);
	maxY = -(FLT_MAX);
	maxZ = -(FLT_MAX);
	numAtoms = 1;
	while (numAtoms != totalNumAtoms) {
		if( (allResiduesArray[(allAtomsArray[numAtoms].fromResidueId)].fromModelId) != (allResiduesArray[(allAtomsArray[(numAtoms + 1)].fromResidueId)].fromModelId)) break;
		numAtoms++;
	}
	// array of the atoms ids sorted by X, Y or Z coordinate
	atomsOrder = (unsigned int *)calloc(numAtoms, sizeof(unsigned int));
	// minimum order among the three X, Y and Z sorting orders for each atom id
	atomsClosenessToSurface = (unsigned int *)malloc((numAtoms + 1) * sizeof(unsigned int));
	for (i = 1; i <= numAtoms; i++) {
		atomsOrder[(i-1)] = i;
		if (allAtomsArray[i].coordX < minX) minX = allAtomsArray[i].coordX;
		if (allAtomsArray[i].coordY < minY) minY = allAtomsArray[i].coordY;
		if (allAtomsArray[i].coordZ < minZ) minZ = allAtomsArray[i].coordZ;
		if (allAtomsArray[i].coordX > maxX) maxX = allAtomsArray[i].coordX;
		if (allAtomsArray[i].coordY > maxY) maxY = allAtomsArray[i].coordY;
		if (allAtomsArray[i].coordZ > maxZ) maxZ = allAtomsArray[i].coordZ;
	}
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
	qsort(atomsOrder, numAtoms, sizeof(unsigned int), SortAtomsCoordsByXValue);
	atomsClosenessToSurface[0] = 0;
	for (i = 0; i <= (numAtoms / 2); i++) {
		atomsClosenessToSurface[(atomsOrder[i])] = i;
		atomsClosenessToSurface[(atomsOrder[(numAtoms - 1) - i])] = i;
	}
	qsort(atomsOrder, numAtoms, sizeof(unsigned int), SortAtomsCoordsByYValue);
	for (i = 0; i <= (numAtoms / 2); i++) {
		k = atomsOrder[i];
		if (i < atomsClosenessToSurface[k]) atomsClosenessToSurface[k] = i;
		k = atomsOrder[(numAtoms - 1) - i];
		if (i < atomsClosenessToSurface[k]) atomsClosenessToSurface[k] = i;
	}
	qsort(atomsOrder, numAtoms, sizeof(unsigned int), SortAtomsCoordsByZValue);
	for (i = 0; i <= (numAtoms / 2); i++) {
		k = atomsOrder[i];
		if (i < atomsClosenessToSurface[k]) atomsClosenessToSurface[k] = i;
		k = atomsOrder[(numAtoms - 1) - i];
		if (i < atomsClosenessToSurface[k]) atomsClosenessToSurface[k] = i;
	}
	k = 0;
	for (i = 1; i <= numAtoms; i++) {
		if (atomsClosenessToSurface[i] > k) k = atomsClosenessToSurface[i];
	}
	svg_shade_level_per_buried_order = (100.0f - 10.0f) / (float)k; // from 100% to 10% color saturation level
	svg_fade_level_per_angstrom_distance = (95.0f - 50.0f) / (maxZ - minZ); // from 50% to 95% color lightness level
	InitializeGraphics("test.svg");
	svg_current_depth_level = (-1);
	for (i = 0; i < numAtoms; i++) {
		DrawSingleAtom(atomsOrder[i]);
	}
	if (svg_current_depth_level != (-1)) fprintf(svgFile, "</g>\n");
	FinalizeGraphics();
	free(svg_opacity_at_depth_level);
	free(svg_blur_size_at_depth_level);
	free(atomsClosenessToSurface);
	free(atomsOrder);
}

// Sort priority (back to front): X, Z, Y
int SortStats(const void *a, const void *b) {
	ForceFieldStats *ffsa = (ForceFieldStats *)a;
	ForceFieldStats *ffsb = (ForceFieldStats *)b;
	float diff = (ffsa->energyY) - (ffsb->energyY);
	if (diff > 0) return (+1);
	if (diff < 0) return (-1);
	diff = (ffsa->energyZ) - (ffsb->energyZ);
	if (diff > 0) return (+1);
	if (diff < 0) return (-1);
	diff = (ffsa->energyX) - (ffsb->energyX);
	if (diff > 0) return (+1);
	else return (-1);
}

// TODO: add #define to convert (x,y,z) pdb coordinates to (X,Y) plot coordinates
int PlotInteratomicForceStats(ForceFieldStats *stats, unsigned int numStats) {
	FILE *imageFile;
	float minX, minY, minZ;
	float maxX, maxY, maxZ;
	float maxAbsX, maxAbsY, maxAbsZ;
	float rangeX, rangeY, rangeZ;
	float minEnergy, maxEnergy, logAbsMinEnergy, logAbsMaxEnergy, rangeEnergy;
	int imageWidth, imageHeight;
	int centerX, centerY;
	float stepX, stepY, stepZ, minStep;
	int k, x, y, size, color, strength, shading;
	float pointValue, opacity;
	unsigned int i;
	int axisWidth, markWidth, markDistance, numMarks, markValue, x2, y2;
	const char backgroundColor[] = "dimgrey";
	const char fontColor[] = "white";
	const int fontWidth = 20;
	const int marginLeft = 100;
	const int marginRight = 150;
	const int marginTop = 100;
	const int marginBottom = 100;
	const int pixelsPerAngstrom = 100;
	const float slopeX = 2.25;
	const float slopeZ = 1.75;
	const int pixelsPerAngstromX = (int)floorf(((2.0f*slopeX - 1.0f)*(float)pixelsPerAngstrom) / (2.0f*slopeX)); // = (ppa-k*(1/slope)*ppa) , with k=(1/2)=0.5
	const int pixelsPerAngstromZ = (int)floorf(((2.0f*slopeZ - 1.0f)*(float)pixelsPerAngstrom) / (2.0f*slopeZ));
	//const int pixelsPerAngstromX = pixelsPerAngstrom;
	//const int pixelsPerAngstromZ = pixelsPerAngstrom;
	minX = FLT_MAX;
	minY = FLT_MAX;
	minZ = FLT_MAX;
	maxX = -(FLT_MAX);
	maxY = -(FLT_MAX);
	maxZ = -(FLT_MAX);
	minEnergy = FLT_MAX;
	maxEnergy = -(FLT_MAX);
	for (i = 0; i < numStats; i++) {
		if (stats[i].energyX < minX) minX = stats[i].energyX;
		if (stats[i].energyX > maxX) maxX = stats[i].energyX;
		if (stats[i].energyY < minY) minY = stats[i].energyY;
		if (stats[i].energyY > maxY) maxY = stats[i].energyY;
		if (stats[i].energyZ < minZ) minZ = stats[i].energyZ;
		if (stats[i].energyZ > maxZ) maxZ = stats[i].energyZ;
		if (stats[i].globalEnergy < minEnergy) minEnergy = stats[i].globalEnergy;
		if (stats[i].globalEnergy > maxEnergy) maxEnergy = stats[i].globalEnergy;
	}
	maxAbsX = maxX;
	if (fabsf(minX) > maxAbsX) maxAbsX = fabsf(minX);
	maxAbsY = maxY;
	if (fabsf(minY) > maxAbsY) maxAbsY = fabsf(minY);
	maxAbsZ = maxZ;
	if (fabsf(minZ) > maxAbsZ) maxAbsZ = fabsf(minZ);
	rangeX = (maxX - minX);
	rangeY = (maxY - minY);
	rangeZ = (maxZ - minZ);
	if (maxEnergy == 0.0) maxEnergy = +1.0; // just to prevent division by zero later on
	if (minEnergy == 0.0) minEnergy = -1.0;
	rangeEnergy = (maxEnergy - minEnergy);
	logAbsMaxEnergy = logf(1.0f + fabsf(maxEnergy)); // +1.0 to prevent negative log values
	logAbsMinEnergy = logf(1.0f + fabsf(minEnergy));
	qsort(stats, numStats, sizeof(ForceFieldStats), SortStats); // sort by drawing order (x, z, y)
	stepX = stats[0].energyX;
	stepY = stats[0].energyY;
	stepZ = stats[0].energyZ;
	i = 1; // get step size of each axis, by going through the sorted list of points and get the value of the next one in that axis
	while ((i < numStats) && (stats[i].energyX == stepX)) i++;
	if (i != numStats) stepX = fabsf(stats[i].energyX - stepX);
	i = 1;
	while ((i < numStats) && (stats[i].energyY == stepY)) i++;
	if (i != numStats) stepY = fabsf(stats[i].energyY - stepY);
	i = 1;
	while ((i < numStats) && (stats[i].energyZ == stepZ)) i++;
	if (i != numStats) stepZ = fabsf(stats[i].energyZ - stepZ);
	minStep = stepX;
	if (stepY < stepX) minStep = stepY;
	if ((stepZ < stepX) && (stepZ < stepY)) minStep = stepZ;
	imageWidth = marginLeft + marginRight + (int)floorf(rangeZ*pixelsPerAngstromZ + rangeX*pixelsPerAngstromX);
	imageHeight = marginTop + marginBottom + (int)floorf((rangeY*pixelsPerAngstrom) + (rangeZ/slopeZ)*pixelsPerAngstromZ + (rangeX/slopeX)*pixelsPerAngstromX);
	centerX = marginLeft + (int)floorf(rangeZ * pixelsPerAngstromZ);
	centerY = marginTop + (int)floorf(rangeY * pixelsPerAngstrom);
	if ((imageFile = fopen("InteratomicForceFieldStats.svg", "w")) == NULL) {
		fprintf(stderr, "> ERROR: Cannot create SVG file\n");
		return (-1);
	}
	fprintf(imageFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %d %d\">\n", imageWidth, imageHeight);
	fprintf(imageFile, "<style type=\"text/css\">\n");
	fprintf(imageFile, "line { stroke: %s; stroke-linecap: round; }\n", fontColor);
	fprintf(imageFile, "text { fill: %s; font-size: %dpx; font-family: sans-serif; text-anchor: middle; dominant-baseline: middle; }\n", fontColor, fontWidth);
	fprintf(imageFile, ".axis { stroke-width: 5; marker-end: url(#arrow); }\n");
	fprintf(imageFile, ".mark { stroke-width: 2; }\n");
	fprintf(imageFile, "circle { transition: all 1s ease; }\n");
	fprintf(imageFile, "circle:hover { stroke: white; stroke-width: 5; opacity: 1; }\n");
	fprintf(imageFile, "</style>\n");
	fprintf(imageFile, "<defs>\n");
	fprintf(imageFile, "<marker id='arrow' markerWidth='20' markerHeight='10' refX='0' refY='5' orient='auto'>\n"); // triangle shaped arrow (pointing to the right)
	fprintf(imageFile, "<path d='M 0 1 L 0 9 L 8 5 Z' fill='%s'/>\n", fontColor);
	fprintf(imageFile, "</marker>\n");
	fprintf(imageFile, "<linearGradient id='grad' x1='0%%' y1='0%%' x2='0%%' y2='100%%'>\n");
	fprintf(imageFile, "<stop offset='0%%' stop-color='blue'/>\n<stop offset='50%%' stop-color='white'/>\n<stop offset='100%%' stop-color='red'/>\n");
	fprintf(imageFile, "</linearGradient>\n");
	fprintf(imageFile, "</defs>\n");
	fprintf(imageFile, "<rect width='100%%' height='100%%' fill='%s'/>\n", backgroundColor);
	// Y axis
	x = centerX;
	y = centerY - (int)floorf(rangeY*pixelsPerAngstrom);
	fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='axis' />\n", centerX, centerY, x, y);
	x -= (fontWidth); // Y axis label
	y -= (fontWidth);
	fprintf(imageFile, "<text x='%d' y='%d'>%s</text>\n", x, y, "Y");
	axisWidth = (int)floorf(rangeY*pixelsPerAngstrom); // size (in pixels) of the full axis
	markWidth = fontWidth * (int)ceilf(log10f(maxAbsY)); // size (in pixels) of each axis mark label
	if (markWidth <= 0) markWidth = fontWidth; // in case the max is 1 or less and the log gives 0 or less
	numMarks = (axisWidth / markWidth);
	markDistance = (int)ceilf(rangeY / (float)numMarks); // size (in angstroms) between consecutive axis marks
	k = (int)powf(10.0, ceilf(log10f((float)markDistance))); // round mark distance to multiple of 10
	if ((k / 5) > markDistance) markDistance = (k / 5); // check if shorter number of type 2*10^k will fit
	else if ((k / 2) > markDistance) markDistance = (k / 2); // check if shorter number of type 5*10^k will fit
	else markDistance = k;
	numMarks = (int)ceilf(minY / (float)markDistance); // to get first mark right next to the minimum value
	while ((float)(markValue = (numMarks*markDistance)) <= maxY) { // draw all marks across the axis
		x = centerX;
		x2 = x + fontWidth; // mark length the same size as a font character
		y = centerY - (int)floorf(pixelsPerAngstrom*((float)markValue - minY));
		y2 = y;
		fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='mark' />\n", x, y, x2, y2);
		x2 += fontWidth; // spacing before the number label
		fprintf(imageFile, "<text x='%d' y='%d'>%d</text>\n", x2, y2, markValue);
		numMarks++;
	}
	// Z axis
	x = centerX - (int)floorf(rangeZ*pixelsPerAngstromZ);
	y = centerY + (int)floorf((rangeZ / slopeZ)*pixelsPerAngstromZ);
	fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='axis' />\n", centerX, centerY, x, y);
	y += (2 * fontWidth); // Z axis label
	x -= (fontWidth);
	fprintf(imageFile, "<text x='%d' y='%d'>%s</text>\n", x, y, "Z");
	axisWidth = (int)floorf(rangeZ*pixelsPerAngstromZ); // Z axis marks
	markWidth = fontWidth * (int)ceilf(log10f(maxAbsZ));
	numMarks = (axisWidth / markWidth);
	markDistance = (int)ceilf(rangeZ / (float)numMarks);
	k = (int)powf(10.0, ceilf(log10f((float)markDistance)));
	if ((k / 5) > markDistance) markDistance = (k / 5);
	else if ((k / 2) > markDistance) markDistance = (k / 2);
	else markDistance = k;
	numMarks = (int)ceilf(minZ / (float)markDistance);
	while ((float)(markValue = (numMarks*markDistance)) <= maxZ) {
		x = centerX - (int)floorf(pixelsPerAngstromZ*((float)markValue - minZ));
		x2 = x - fontWidth;
		y = centerY + (int)floorf(pixelsPerAngstromZ*(((float)markValue - minZ)) / slopeZ);
		y2 = y - (int)floorf((float)fontWidth / slopeX); // opposite slope
		fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='mark' />\n", x, y, x2, y2);
		y2 -= fontWidth;
		fprintf(imageFile, "<text x='%d' y='%d'>%d</text>\n", x2, y2, markValue);
		numMarks++;
	}
	// X axis
	x = centerX + (int)floorf(rangeX*pixelsPerAngstromX);
	y = centerY + (int)floorf((rangeX / slopeX)*pixelsPerAngstromX);
	fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='axis' />\n", centerX, centerY, x, y);
	y += (2 * fontWidth); // X axis label
	x += (fontWidth);
	fprintf(imageFile, "<text x='%d' y='%d'>%s</text>\n", x, y, "X");
	axisWidth = (int)floorf(rangeX*pixelsPerAngstromX); // X axis marks
	markWidth = fontWidth * (int)ceilf(log10f(maxAbsX));
	numMarks = (axisWidth / markWidth);
	markDistance = (int)ceilf(rangeX / (float)numMarks);
	k = (int)powf(10.0, ceilf(log10f((float)markDistance)));
	if ((k / 5) > markDistance) markDistance = (k / 5);
	else if ((k / 2) > markDistance) markDistance = (k / 2);
	else markDistance = k;
	numMarks = (int)ceilf(minX / (float)markDistance);
	while ((float)(markValue = (numMarks*markDistance)) <= maxX) {
		x = centerX + (int)floorf(pixelsPerAngstromX*((float)markValue - minX));
		x2 = x + fontWidth;
		y = centerY + (int)floorf(pixelsPerAngstromX*(((float)markValue - minX)) / slopeX);
		y2 = y - (int)floorf((float)fontWidth / slopeZ); // opposite slope
		fprintf(imageFile, "<line x1='%d' y1='%d' x2='%d' y2='%d' class='mark' />\n", x, y, x2, y2);
		y2 -= fontWidth;
		fprintf(imageFile, "<text x='%d' y='%d'>%d</text>\n", x2, y2, markValue);
		numMarks++;
	}
	// rectangle with color legend
	x = imageWidth - ((marginRight - fontWidth) / 2);
	y = centerY - (10 * fontWidth / 2); // with size (2*fw)x(10*fw)
	fprintf(imageFile, "<rect x='%d' y='%d' width='%d' height='%d' fill='url(#grad)' stroke='%s' stroke-width='1' />\n", x, y, (2 * fontWidth), (10 * fontWidth), fontColor);
	x += ((2*fontWidth) / 2); // centered text of size (fw) above rectangle of width (2*fw)
	y -= (fontWidth / 2);
	fprintf(imageFile, "<text x='%d' y='%d'>%.0f</text>\n", x, y, maxEnergy);
	y += (10 * fontWidth) + (fontWidth + (fontWidth / 2)); // centered text of size (fw) bellow rectangle of height (10*fw)
	fprintf(imageFile, "<text x='%d' y='%d'>%.0f</text>\n", x, y, minEnergy);
	// draw radius (shorter to larger), color (white to red) and opacity (transparent to opaque) according to decreasing energy value (the lower the better)
	for (i = 0; i < numStats; i++) { // draw all points
		pointValue = (fabsf((stats[i].energyX - minX) / rangeX) + fabsf((stats[i].energyY - minY) / rangeY) + fabsf((stats[i].energyZ - minZ) / rangeZ)) / 3.0f; // just to calculate shading ; distance from origin, from 0.0 to 1.0
		shading = 50 + (int)floorf(pointValue * 50); //  darker (50%) near the origin (0.0)
		/*
		pointValue = (maxEnergy - stats[i].globalEnergy) / rangeEnergy; // lower energy is better (1.0), higher energy is worse (0.0)
		color = 100 - (int)floorf(pointValue * 50); // minimum energy (1.0) is red (50), and maximum energy (0.0) is white (100)
		*/
		if (stats[i].globalEnergy >= 0) { // blue for positive values and red for negative values
			pointValue = (logAbsMaxEnergy - logf(1.0f + fabsf(stats[i].globalEnergy))) / logAbsMaxEnergy; // logarithmic values ; lower absolute energy is better (1.0), higher energy is worse (0.0)
			color = 0; // red has hue=0
		}
		else {
			pointValue = (logAbsMinEnergy - logf(1.0f + fabsf(stats[i].globalEnergy))) / logAbsMinEnergy;
			color = 240; // blue has hue=240
		}
		size = (int)floorf((pointValue * minStep / 2) * pixelsPerAngstrom); // maximum radius is half the minimum step among all axis
		strength = 50 + (int)floorf(pointValue * 50); // minimum energy (1.0) is white (100), and maximum energy (0.0) is red/blue (50)
		opacity = 0.1f + 0.9f*pointValue; // minimum energy (1.0) is fully opaque, and maximum energy (0.0) is almost fully transparent
		// X = centerX - (z-minZ)*ppa + (x-minX)*ppa ; (z to the left, x to the right)
		x = centerX + (int)floorf((stats[i].energyX - minX)*pixelsPerAngstromX - (stats[i].energyZ - minZ)*pixelsPerAngstromZ);
		// Y = centerY - (y-minY)*ppa + ((z-minZ)/slopeZ)*ppa + ((x-minX)/slopeX)*ppa ; (y up, z and x down)
		y = centerY + (int)floorf(((stats[i].energyZ - minZ) / slopeZ)*pixelsPerAngstromZ + ((stats[i].energyX - minX) / slopeX)*pixelsPerAngstromX - (stats[i].energyY - minY)*pixelsPerAngstrom);
		// hue=(0/240) (red or blue), saturation=(25%-100%) (from almost grey to full color), lightness=(50%-100%) (from red/blue to white)
		fprintf(imageFile, "<circle cx='%d' cy='%d' r='%d' fill='hsl(%d,%d%%,%d%%)' opacity='%.3f' />\n", x, y, size, color, shading, strength, opacity);
	}
	fprintf(imageFile, "</svg>\n");
	if (fclose(imageFile) == EOF) {
		fprintf(stderr, "> ERROR: Cannot write SVG file\n");
		return (-1);
	}
	return 0;
}
