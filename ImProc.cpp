/*
  First steps in Image Processing
  Creating simple Synthetic Images
  Using BMP Library
*/


#include <stdio.h> // for printf
#include <conio.h> // for getch
#include <string.h>
#include <cstring>
#include <iostream> // for cin cout
#include <fstream>  // For file IO
using namespace std; // Explain some day 

// This needed to do Math calculations
#define _USE_MATH_DEFINES
#include <math.h>

// BMP Library
#include "ImProcInPlainC.h"
#include "PrimeFFTn.h" 

#define MAX_VAL 256
#define WHITE_VAL 255
#define BLACK_VAL 0



tFloat ReIm[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
tFloat ImIm[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
tFloat filter[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];


unsigned char firstDerIm[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
unsigned char secondDerIm[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];

unsigned char RGBImage1[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS][NUMBER_OF_COLORS];
unsigned char RGBImage2[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS][NUMBER_OF_COLORS];

unsigned char imgH[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
unsigned char imgS[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
unsigned char imgL[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];



void RGBtoHSL(unsigned char RGBImage[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS][NUMBER_OF_COLORS], unsigned char imgH[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS]
	, unsigned char imgS[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], unsigned char imgL[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS])
{
	for (int row = 0; row < NUMBER_OF_ROWS; row++)
		for (int col = 0; col < NUMBER_OF_COLUMNS; col++)
		{
			double RHelp, GHelp, BHelp, CMax, CMin, Delta;

			RHelp = (double)(RGBImage[row][col][R] / (double)WHITE_VAL);
			GHelp = (double)(RGBImage[row][col][G] / (double)WHITE_VAL);
			BHelp = (double)(RGBImage[row][col][B] / (double)WHITE_VAL);
			CMax = max(max(RHelp, GHelp), BHelp);
			CMin = min(min(RHelp, GHelp), BHelp);
			Delta = CMax - CMin;
			imgL[row][col] = ((double)(CMax + CMin) / (double)2) * WHITE_VAL;
			if (Delta != 0)
			{
				imgS[row][col] = double(Delta / double(1 - fabs(CMax + CMin - 1)))* WHITE_VAL;
				if (CMax == RHelp)
					imgH[row][col] = double(60 * (fmodl((GHelp - BHelp) / (double)Delta, 6)))* double(WHITE_VAL/360.0);
				else if (CMax == GHelp)
					imgH[row][col] = double(60 * ((BHelp - RHelp) / (double)(Delta)+2)) * double(WHITE_VAL / 360.0);
				else
					imgH[row][col] = double(60 * ((RHelp - GHelp) / (double)(Delta)+4)) * double(WHITE_VAL / 360.0);
			}
			else
			{
				imgS[row][col] = 0;
				imgH[row][col] = 0;
			}
		}
}


void HSLtoRGB(int HSLImage[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS][NUMBER_OF_COLORS], unsigned char RGBImage[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS][NUMBER_OF_COLORS])
{
	for (int row = 0; row < NUMBER_OF_ROWS; row++)
		for (int col = 0; col < NUMBER_OF_COLUMNS; col++)
		{
			double RHelp, GHelp, BHelp, C, X, m;
			C = (1 - abs(2.0 * HSLImage[row][col][2] / 100.0 - 1)) * (HSLImage[row][col][1] / 100.0);
			X = C * (1 - abs(fmodl(HSLImage[row][col][0] / 60.0, 2) - 1));
			m = HSLImage[row][col][2] / 100.0 - C / 2.0;
			if (HSLImage[row][col][0] < 60)
			{
				RHelp = C;
				GHelp = X;
				BHelp = 0;
			}
			else if (HSLImage[row][col][0] < 120)
			{
				RHelp = X;
				GHelp = C;
				BHelp = 0;
			}
			else if (HSLImage[row][col][0] < 180)
			{
				RHelp = 0;
				GHelp = C;
				BHelp = X;
			}
			else if (HSLImage[row][col][0] < 240)
			{
				RHelp = 0;
				GHelp = X;
				BHelp = C;
			}
			else if (HSLImage[row][col][0] < 300)
			{
				RHelp = X;
				GHelp = 0;
				BHelp = C;
			}
			else
			{
				RHelp = C;
				GHelp = 0;
				BHelp = X;
			}
			RGBImage[row][col][R] = (RHelp + m) * 255;
			RGBImage[row][col][G] = (GHelp + m) * 255;
			RGBImage[row][col][B] = (BHelp + m) * 255;
		}
}

void Threshold(unsigned char img[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], unsigned char Threshold)
{
	unsigned char Threshold_LUT[256];      
	for (int i = 0; i < 256; i++)
	{
		Threshold_LUT[i] = 0;
		if (i < Threshold)			 
			Threshold_LUT[i] = 255;	
	}
	unsigned char* ptrtopixel = img[0];
	for (int pixel = 0; pixel < NUMBER_OF_ROWS * NUMBER_OF_COLUMNS; pixel++)
		*ptrtopixel++ = Threshold_LUT[*ptrtopixel]; //applying the LUT on the entire picture
}


void RobertOperator(unsigned char sour[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], unsigned char dest[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS]) {
	for (int row = 1; row < NUMBER_OF_ROWS - 1; row++)
		for (int col = 1; col < NUMBER_OF_COLUMNS - 1; col++)
		{
			int Gx = +(sour[row - 1][col - 1]) - (sour[row + 1][col + 1]);
			int Gy = +(sour[row + 1][col + 1]) - (sour[row - 1][col - 1]);
			int Grad = abs(Gx) + abs(Gy);
			dest[row][col] = (unsigned char)Grad;
		}
}


void ApplyEdgeDetection(unsigned char sour[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], int threshold,char* name)
{
	unsigned char tmp[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
	RobertOperator(sour, tmp);
	Threshold(tmp, threshold);
	StoreGrayImageAsGrayBmpFile(tmp, name);
}


void Convert(tFloat floatImage[][NUMBER_OF_COLUMNS],
	unsigned char byteOriginal[][NUMBER_OF_COLUMNS],
	tFloat minVal,
	tFloat maxVal)
{
	tFloat* ptrToFloatImage = floatImage[0];
	tFloat c, b, temp;
	if (minVal == maxVal)
	{
		c = 1; b = 0;
	}
	else
	{
		c = 255.0 / (maxVal - minVal);
		b = -c * minVal;
	}

	unsigned char* ptrToByteImage = byteOriginal[0];
	ptrToFloatImage = floatImage[0];

	for (int pixel = 0; pixel < NUMBER_OF_COLUMNS * NUMBER_OF_ROWS; pixel++)
	{
		temp = c * (*ptrToFloatImage++) + b;
		if (temp < 0)    temp = 0;
		if (temp > 255)  temp = 255;
		*ptrToByteImage++ = (unsigned char)((int)temp);
	}
}

void Convert(unsigned char byteOriginal[][NUMBER_OF_COLUMNS],
	tFloat floatImage[][NUMBER_OF_COLUMNS])
{
	unsigned char* ptrToByteImage = byteOriginal[0];
	tFloat* ptrToFloatImage = floatImage[0];

	for (int pixel = 0; pixel < NUMBER_OF_COLUMNS * NUMBER_OF_ROWS; pixel++)
	{
		*ptrToFloatImage++ = (tFloat)*ptrToByteImage++;
	}
}

void Clear(tFloat Image[][NUMBER_OF_COLUMNS], unsigned char value)
{
	tFloat* ptr = Image[0];
	for (int i = 0; i < NUMBER_OF_COLUMNS * NUMBER_OF_ROWS; i++) { *(ptr++) = value; }
}

void ShiftHalfSize(tFloat floatImage[][NUMBER_OF_COLUMNS])
{
	for (int row = 0; row < NUMBER_OF_ROWS; row++)
	{
		for (int column = 0; column < NUMBER_OF_COLUMNS; column += 2)
		{
			floatImage[row][column] *= -1;
		}
	}

	for (int row = 0; row < NUMBER_OF_ROWS; row += 2)
	{
		for (int column = 0; column < NUMBER_OF_COLUMNS; column++)
		{
			floatImage[row][column] *= -1;
		}
	}
}

void DoFFT(tFloat imageRe[][NUMBER_OF_COLUMNS],
	tFloat imageIm[][NUMBER_OF_COLUMNS],
	int exponentSign,
	int normalizationScalingType)
{
	tFloat* ptrToRe = imageRe[0];
	tFloat* ptrToIm = imageIm[0];

	tInteger arrayOfDimensions[2];
	arrayOfDimensions[0] = NUMBER_OF_COLUMNS;
	arrayOfDimensions[1] = NUMBER_OF_ROWS;

	PrimeFFTn(2,
		arrayOfDimensions,
		ptrToRe,
		ptrToIm,
		exponentSign,
		normalizationScalingType);
}

void DoFiltrationInFD(tFloat floatRe[][NUMBER_OF_COLUMNS],
	tFloat floatIm[][NUMBER_OF_COLUMNS],
	tFloat floatFilter[][NUMBER_OF_COLUMNS])
{
	tFloat* ptrToRe = floatRe[0];
	tFloat* ptrToIm = floatIm[0];
	tFloat* ptrToFilter = floatFilter[0];

	for (int pixel = 0; pixel < NUMBER_OF_ROWS * NUMBER_OF_COLUMNS; pixel++)
	{
		*ptrToRe++ *= *ptrToFilter;
		*ptrToIm++ *= *ptrToFilter++;
	}

}

void OptimalConvert(tFloat floatImage[][NUMBER_OF_COLUMNS],
	unsigned char byteOriginal[][NUMBER_OF_COLUMNS])
{
	tFloat* ptrToFloatImage = floatImage[0];
	tFloat MinVal, MaxVal, temp;
	// Find MinMax
	MinVal = *ptrToFloatImage;
	MaxVal = MinVal;

	for (int pixel = 0; pixel < NUMBER_OF_COLUMNS * NUMBER_OF_ROWS; pixel++)
	{
		temp = *ptrToFloatImage++;
		if (MinVal > temp) MinVal = temp;
		if (MaxVal < temp) MaxVal = temp;
	}
	Convert(floatImage, byteOriginal, MinVal, MaxVal);
}

void createClippedFirstDerivativeFilter(tFloat floatFilter[][NUMBER_OF_COLUMNS], tFloat sigmaRadius)
{
	tFloat x, y, r, z;
	for (size_t row = 0; row < NUMBER_OF_ROWS; row++)
	{
		for (size_t column = 0; column < NUMBER_OF_COLUMNS; column++)
		{
			y = float(row) - (NUMBER_OF_ROWS / 2);
			x = float(column) - (NUMBER_OF_COLUMNS / 2);
			r = sqrt(y * y + x * x);
			if (r > sigmaRadius)
			{
				z = sigmaRadius;
			}
			else
				z = r;
			floatFilter[row][column] = z;
		}
	}
}


void DoTwoThresholdsCriteria(unsigned char first[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], unsigned char second[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS],
	unsigned char FirstCriteria, unsigned char SecondCriteria, unsigned char dest[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS])
{
	unsigned char* ptrf = first[0];
	unsigned char* ptrs = second[0];
	unsigned char* ptrDest = dest[0];
	for (int i = 0; i < NUMBER_OF_COLUMNS * NUMBER_OF_ROWS; i++)
		if (*(ptrf + i) > FirstCriteria && *(ptrs + i) > SecondCriteria)
			*(ptrDest + i) = WHITE_VAL;
		else 
			*(ptrDest + i) = BLACK_VAL;
}


void TwoCriteriaFirstSecondFilterFD(unsigned char sour[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS], unsigned char FirstCriteria, unsigned char SecondCriteria, char* name)
{
	Convert(sour, ReIm);
	Clear(ImIm, 0);
	ShiftHalfSize(ReIm);
	ShiftHalfSize(ImIm);
	DoFFT(ReIm, ImIm, FORWARD_FFT, NORMALIZE_BY_SQRT);
	createClippedFirstDerivativeFilter(filter, 60);
	DoFiltrationInFD(ReIm, ImIm, filter);
	DoFFT(ReIm, ImIm, REVERSE_FFT, NORMALIZE_BY_SQRT);
	ShiftHalfSize(ReIm);
	ShiftHalfSize(ImIm);
	OptimalConvert(ReIm, firstDerIm);

	Convert(sour, ReIm);
	Clear(ImIm, 0);
	ShiftHalfSize(ReIm);
	ShiftHalfSize(ImIm);
	DoFFT(ReIm, ImIm, FORWARD_FFT, NORMALIZE_BY_SQRT);
	DoFiltrationInFD(ReIm, ImIm, filter);
	DoFiltrationInFD(ReIm, ImIm, filter);
	DoFFT(ReIm, ImIm, REVERSE_FFT, NORMALIZE_BY_SQRT);
	ShiftHalfSize(ReIm);
	ShiftHalfSize(ImIm);
	OptimalConvert(ReIm, secondDerIm);

	unsigned char tmp[NUMBER_OF_ROWS][NUMBER_OF_COLUMNS];
	DoTwoThresholdsCriteria(firstDerIm, secondDerIm, FirstCriteria, SecondCriteria,tmp);
	//StoreGrayImageAsGrayBmpFile(secondDerIm, name);
	StoreGrayImageAsGrayBmpFile(tmp, name);
}


void main()
{

	//------------- img1 ---------------//
	LoadBgrImageFromTrueColorBmpFile(RGBImage1, "Img1.bmp");
	RGBtoHSL(RGBImage1, imgH, imgS, imgL);
	StoreGrayImageAsGrayBmpFile(imgH, "Img1H.bmp");
	StoreGrayImageAsGrayBmpFile(imgS, "Img1S.bmp");
	StoreGrayImageAsGrayBmpFile(imgL, "Img1L.bmp");

	//----- Robert edge detector (alg1) -------//
	ApplyEdgeDetection(imgH, 8, "Img1HAlg1.bmp");
	ApplyEdgeDetection(imgS, 5, "Img1SAlg1.bmp");
	ApplyEdgeDetection(imgL, 10, "Img1LAlg1.bmp");

	//----- First & Second derivative filter (alg2)------//
	TwoCriteriaFirstSecondFilterFD(imgH, 130, 150, "Img1HAlg2.bmp");
	TwoCriteriaFirstSecondFilterFD(imgS, 100, 130, "Img1SAlg2.bmp");
	TwoCriteriaFirstSecondFilterFD(imgL, 115, 148, "Img1LAlg2.bmp");

	//------------- img2 ---------------//
	LoadBgrImageFromTrueColorBmpFile(RGBImage2, "Img2.bmp");
	RGBtoHSL(RGBImage2, imgH, imgS, imgL);
	StoreGrayImageAsGrayBmpFile(imgH, "Img2H.bmp");
	StoreGrayImageAsGrayBmpFile(imgS, "Img2S.bmp");
	StoreGrayImageAsGrayBmpFile(imgL, "Img2L.bmp");

	//----- Robert edge detector (alg1) -------//
	ApplyEdgeDetection(imgH, 50, "Img2HAlg1.bmp");
	ApplyEdgeDetection(imgS, 70, "Img2SAlg1.bmp");
	ApplyEdgeDetection(imgL, 100, "Img2LAlg1.bmp");

	//----- First & Second derivative filter (alg2)------//
	TwoCriteriaFirstSecondFilterFD(imgH, 118, 120, "Img2HAlg2.bmp");
	TwoCriteriaFirstSecondFilterFD(imgS, 95, 113, "Img2SAlg2.bmp");
	TwoCriteriaFirstSecondFilterFD(imgL, 73, 67, "Img2LAlg2.bmp");

	WaitForUserPressKey();
}