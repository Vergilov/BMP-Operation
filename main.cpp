#pragma warning( disable : 4018 4244)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct BitMapHeader
{
    unsigned short Signature;
    unsigned long FileSize;
    unsigned long Reserved;
    unsigned long DataOffSet;
} BitMapHeader;

typedef struct BitMapInfo
{
    unsigned long Size;
    unsigned long Width;
    unsigned long Height;
    unsigned short Planes;
    unsigned short BitsPerPixel;
    unsigned long Compression;
    unsigned long SizeImage;
    unsigned long XPixelsPreMeter;
    unsigned long YPixelsPreMeter;
    unsigned long ColorsUsed;
    unsigned long ColorsImportant;


    unsigned long* ColorTable;
    unsigned char* pixels;
} BitMapInfo;


typedef struct BitMap
{
    BitMapInfo info;
    BitMapHeader header;
} BitMap;

typedef struct Pix
{
    unsigned char R;
    unsigned char G;
    unsigned char B;
    unsigned char L;
    int BW;
} Pix;

double PI = 3.1415927;
double RMAX = RAND_MAX;

// cw1
BitMap readBitMap(FILE *BMPFile);
void writeBitMap(FILE *BMPFile, BitMap bm);
Pix getPix(int x, int y, BitMap bm);
void putPix(int x, int y, BitMap bm, Pix p);

// cw2
int* histogramL(BitMap bm);
BitMap addConst(BitMap bm, int a);
BitMap gray(BitMap bm);
BitMap add(BitMap bm1, BitMap bm2);
BitMap max(BitMap bm1, BitMap bm2);
int* histogramR(BitMap bm);
int* histogramG(BitMap bm);
int* histogramB(BitMap bm);
BitMap multiplyConst(BitMap bm, float a);
BitMap negative(BitMap bm);
BitMap threshold(BitMap bm, int level);
BitMap fixLightness(BitMap bm);
int max(BitMap bm);
int min(BitMap bm);
float mean(BitMap bm);
int median(BitMap bm);
BitMap difference(BitMap bm1, BitMap bm2);
BitMap multiply(BitMap bm1, BitMap bm2);
BitMap min(BitMap bm1, BitMap bm2);
BitMap linearCombination(BitMap bm1, BitMap bm2, float weight1, float weight2);

// cw3
BitMap move(BitMap bm, int wl, int hl, int cycle);
BitMap crop(BitMap bm, int w, int h, int x,int y);
BitMap resize(BitMap bm, int w, int h, const char* method);
BitMap rotate(BitMap bm, int deg);
BitMap flip(BitMap bm, int direct);

// cw4
BitMap filter(BitMap bm, float** matrix, int w, int h);
BitMap dilate(BitMap bm, int** strEl, int w, int h);
BitMap erode(BitMap bm,int ** strEl,int w,int h);
BitMap openMor(BitMap bm, int** strEl,int w, int h);
BitMap closeMor(BitMap bm, int** strEl, int w, int h);
BitMap whiteTH(BitMap bm, int ** strEl, int w, int h);
BitMap blackTH(BitMap bm, int ** strEl, int w, int h);

// cw5
BitMap addNoiseRec(BitMap bm, int level, float density);

BitMap filterMean(BitMap bm);
BitMap filterGauss(BitMap bm);
BitMap filterMeanA(BitMap bm);
BitMap filterMeanB(BitMap bm);

BitMap filterRoberts(BitMap bm);
BitMap filterPrewitt(BitMap bm);
BitMap filterSobel(BitMap bm);
BitMap filterLaplace(BitMap bm);

BitMap filterMedian(BitMap bm, int w, int h);

BitMap addNoiseNor(BitMap bm, int mean, float sigma, float density);
BitMap addNoiseSP(BitMap bm, float density);


// pomocnicze funkcje

int align(int v)
{
    if (v < 0) return 0;
    else if (v > 255) return 255;
    else return v;
}

int min(int a, int b)
{
    return a < b ? a : b;
}
int max(int a, int b)
{
    return a > b ? a : b;
}

int cmp (const void *av, const void *bv)
{
    int a = *(int*)av;
    int b = *(int*)bv;
    if(a < b) return -1;
    else if(a == b) return 0;
    else return 1;
}

double normalDistGenerator(const double mean, const double sigma)
{
	return sigma * sqrt(-2*log((rand()+1.0f)/RAND_MAX))*sin(2*PI*rand()/RMAX)+mean;
}

float** makeMask3x3(float values[])
{
    float** mask = (float**) malloc (3*sizeof(float*));
    for (int i = 0; i < 3; i++)
    {
        mask[i] = (float*)malloc(3 * sizeof(float));
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mask[i][j] = values[3 * i + j];
        }
    }
    return mask;
}

int** makeStrEl3x3(int values[])
{
    int** mask = (int**) malloc (3*sizeof(int*));
    for (int i = 0; i < 3; i++)
    {
        mask[i] = (int*)malloc(3 * sizeof(int));
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mask[i][j] = values[3 * i + j];
        }
    }
    return mask;

}

int randRange(int a, int b)
{
	return a + (rand() % (int)(b - a + 1));
}
double rand01()
{
    return (double)rand() / (double)RAND_MAX ;
}

BitMap makeCopy(BitMap bm)
{
    int i=0;
    BitMap output;
    output.header.DataOffSet = bm.header.DataOffSet;
    output.header.FileSize = bm.header.FileSize;
    output.header.Reserved = bm.header.Reserved;
    output.header.Signature = bm.header.Signature;
    output.info.BitsPerPixel = bm.info.BitsPerPixel;
    output.info.ColorsImportant = bm.info.ColorsImportant;
    output.info.ColorsUsed = bm.info.ColorsUsed;
    output.info.Compression = bm.info.Compression;
    output.info.Height = bm.info.Height;
    output.info.Planes = bm.info.Planes;
    output.info.Size = bm.info.Size;
    output.info.SizeImage = bm.info.SizeImage;
    output.info.Width = bm.info.Width;
    output.info.XPixelsPreMeter = bm.info.XPixelsPreMeter;
    output.info.YPixelsPreMeter = bm.info.YPixelsPreMeter;
    if(bm.header.DataOffSet - 54 > 0)
    {
        output.info.ColorTable = (unsigned long*)malloc(sizeof(unsigned long)*(output.header.DataOffSet-54)/4);
        int colorCount;
        if(output.info.ColorsUsed == 0 && output.info.BitsPerPixel == 8)
        {
            colorCount = 256;
        }
        else if(output.info.ColorsUsed == 0 && output.info.BitsPerPixel == 1)
        {
            colorCount = 2;
        }
        for(i = 0; i<colorCount; i++)
        {
            output.info.ColorTable[i] = bm.info.ColorTable[i];
        }
    }
    output.info.pixels = (unsigned char*)malloc(sizeof(unsigned char)*(bm.info.SizeImage));
    for(i = 0; i<output.info.SizeImage; i++)
    {
        output.info.pixels[i] = bm.info.pixels[i];
    }
    return output;
}

BitMap setSize(int newWidth, int newHeight, BitMap bm)
{
    BitMap output = makeCopy(bm);
    output.info.Height = newHeight;
    output.info.Width = newWidth;
    if(bm.info.BitsPerPixel == 24)
        output.info.SizeImage = newWidth*newHeight*3;
    else if(bm.info.BitsPerPixel == 8)
        output.info.SizeImage = newWidth*newHeight;
    free(output.info.pixels);
    output.info.pixels = (unsigned char*)malloc(sizeof(unsigned char)*(output.info.SizeImage));
    output.info.Size = output.header.DataOffSet - 14;
    return output;
}


///////////////////////////////
// cw1

BitMap readBitMap(FILE *BMPFile)
{
    BitMap bm;
    if(!BMPFile)
    {
        printf("Nie udalo sie otworzyc pliku\n");
        exit(EXIT_FAILURE);
    }

    fseek(BMPFile, 0, SEEK_SET);

    fread(&bm.header.Signature,sizeof(bm.header.Signature), 1, BMPFile);
    fread(&bm.header.FileSize,sizeof(bm.header.FileSize), 1, BMPFile);
    fread(&bm.header.Reserved,sizeof(bm.header.Reserved), 1, BMPFile);
    fread(&bm.header.DataOffSet,sizeof(bm.header.DataOffSet), 1, BMPFile);

    fread(&bm.info.Size,sizeof(bm.info.Size), 1, BMPFile);
    fread(&bm.info.Width,sizeof(bm.info.Width), 1, BMPFile);
    fread(&bm.info.Height,sizeof(bm.info.Height), 1, BMPFile);
    fread(&bm.info.Planes,sizeof(bm.info.Planes), 1, BMPFile);
    fread(&bm.info.BitsPerPixel,sizeof(bm.info.BitsPerPixel), 1, BMPFile);
    fread(&bm.info.Compression,sizeof(bm.info.Compression), 1, BMPFile);
    fread(&bm.info.SizeImage,sizeof(bm.info.SizeImage), 1, BMPFile);
    fread(&bm.info.XPixelsPreMeter,sizeof(bm.info.XPixelsPreMeter), 1, BMPFile);
    fread(&bm.info.YPixelsPreMeter,sizeof(bm.info.YPixelsPreMeter), 1, BMPFile);
    fread(&bm.info.ColorsUsed,sizeof(bm.info.ColorsUsed), 1, BMPFile);
    fread(&bm.info.ColorsImportant,sizeof(bm.info.ColorsImportant), 1, BMPFile);
    if(bm.header.DataOffSet - 54 > 0)
    {
        bm.info.ColorTable = (unsigned long*)malloc(sizeof(unsigned long)*(bm.header.DataOffSet-54)/4);
        fread(bm.info.ColorTable,sizeof(unsigned long), (bm.header.DataOffSet-54)/4, BMPFile);
    }

    bm.info.pixels = (unsigned char*)malloc(sizeof(unsigned char)*(bm.info.SizeImage));
    fread(bm.info.pixels,sizeof(unsigned char), bm.info.SizeImage , BMPFile);
    return bm;
}

void writeBitMap(FILE *BMPFile, BitMap bm)
{

    fseek(BMPFile, 0, SEEK_SET);
    fwrite(&bm.header.Signature,sizeof(bm.header.Signature), 1, BMPFile);
    fwrite(&bm.header.FileSize,sizeof(bm.header.FileSize), 1, BMPFile);
    fwrite(&bm.header.Reserved,sizeof(bm.header.Reserved), 1, BMPFile);
    fwrite(&bm.header.DataOffSet,sizeof(bm.header.DataOffSet), 1, BMPFile);

    fwrite(&bm.info.Size,sizeof(bm.info.Size), 1, BMPFile);
    fwrite(&bm.info.Width,sizeof(bm.info.Width), 1, BMPFile);
    fwrite(&bm.info.Height,sizeof(bm.info.Height), 1, BMPFile);
    fwrite(&bm.info.Planes,sizeof(bm.info.Planes), 1, BMPFile);
    fwrite(&bm.info.BitsPerPixel,sizeof(bm.info.BitsPerPixel), 1, BMPFile);
    fwrite(&bm.info.Compression,sizeof(bm.info.Compression), 1, BMPFile);
    fwrite(&bm.info.SizeImage,sizeof(bm.info.SizeImage), 1, BMPFile);
    fwrite(&bm.info.XPixelsPreMeter,sizeof(bm.info.XPixelsPreMeter), 1, BMPFile);
    fwrite(&bm.info.YPixelsPreMeter,sizeof(bm.info.YPixelsPreMeter), 1, BMPFile);
    fwrite(&bm.info.ColorsUsed,sizeof(bm.info.ColorsUsed), 1, BMPFile);
    fwrite(&bm.info.ColorsImportant,sizeof(bm.info.ColorsImportant), 1, BMPFile);

    if(bm.header.DataOffSet - 54 > 0)
        fwrite(bm.info.ColorTable,sizeof(unsigned long), (bm.header.DataOffSet-54)/4, BMPFile);

    fwrite(bm.info.pixels, sizeof(unsigned char), bm.info.SizeImage , BMPFile);
}

Pix getPix(int x, int y, BitMap bm)
{
    Pix p;
    if(bm.info.BitsPerPixel == 24)
    {
        int Pozycja = (x+y*bm.info.Width)*3;
        if((3*bm.info.Width)%4 > 0)
            Pozycja += (4-(3*bm.info.Width)%4)*y;

        p.B = bm.info.pixels[Pozycja];
        p.G = bm.info.pixels[Pozycja+1];
        p.R = bm.info.pixels[Pozycja+2];
    }
    else if(bm.info.BitsPerPixel == 8)
    {
        int Pozycja = (x+y*bm.info.Width);
        if((3*bm.info.Width)%4 > 0)
            Pozycja += (4-bm.info.Width%4)*y;

        p.BW = 1;
        p.L = bm.info.pixels[Pozycja];
    }
    else if(bm.info.BitsPerPixel == 1)
    {
        int i = 0;
        int Pozycja = x/8+(y*(bm.info.Width/8));

        int Padding = bm.info.Width/8;
        if(Padding < bm.info.Width/8.0)
            Padding++;

        if(Padding%4 > 0)
        {
            Pozycja += y;
            Pozycja += (4-Padding%4)*y;
        }

        int Bit = x%8;
        unsigned char pBit = 128;
        for(i = 0; i < Bit; i++)
            pBit = pBit >> 1;


        p.BW = 1;
        p.L = 0;
        if((bm.info.pixels[Pozycja] & pBit) == pBit)
            p.L = 1;

    }
    return p;
}

void putPix(int x, int y, BitMap bm, Pix p)
{
    if(bm.info.BitsPerPixel == 24)
    {
        int Pozycja = (x+y*bm.info.Width)*3;
        if((3*bm.info.Width)%4 > 0)
            Pozycja += (4-(3*bm.info.Width)%4)*y;

        bm.info.pixels[Pozycja] = p.B;
        bm.info.pixels[Pozycja+1] = p.G;
        bm.info.pixels[Pozycja+2] = p.R;
    }
    else if(bm.info.BitsPerPixel == 8)
    {
        int Pozycja = (x+y*bm.info.Width);
        if(bm.info.Width%4 > 0)
            Pozycja += (4-bm.info.Width%4)*y;

        bm.info.pixels[Pozycja] = p.L;
    }
    else if(bm.info.BitsPerPixel == 1)
    {
        int i = 0;
        int Pozycja = x/8+(y*(bm.info.Width/8));

        int Padding = bm.info.Width/8;
        if(Padding < bm.info.Width/8.0)
            Padding++;

        if(Padding%4 > 0)
        {
            Pozycja += y;
            Pozycja += (4-Padding%4)*y;
        }

        int Bit = x%8; 
        unsigned char pBit = 128;
        for(i = 0; i < Bit; i++)
            pBit = pBit >> 1;

    }
}

///////////////////////////////
// cw2

int* histogramL(BitMap bm)
{
    int* histogram = (int*)malloc(sizeof(int)*256);
    for (int i = 0; i < 256; i++) histogram[i] = 0;
    for (unsigned y = 0; y < bm.info.Height; y++)
    {
        for (unsigned x = 0; x < bm.info.Width; x++)
        {
            histogram[getPix(x, y, bm).L]++;
        }
    }
    return histogram;
}

BitMap addConst(BitMap bm, int a)
{
    BitMap output = makeCopy(bm);
    int x,y;
    switch(bm.info.BitsPerPixel)
    {
    case 8:
    {
        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x, y, bm);
                if(p.L + a > 255)
                    p.L = 255;
                else if(p.R + a < 0)
                    p.L = 0;
                else
                    p.L = p.L + a;
                putPix(x,y,bm,p);
            }
        }
    }
    case 24:
    {
        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x,y,bm);
                if(p.R + a > 255)
                    p.R = 255;
                else if(p.R + a < 0)
                    p.R = 0;
                else
                    p.R = p.R + a;

                if(p.G + a > 255)
                    p.G = 255;
                else if(p.G + a < 0)
                    p.G = 0;
                else
                    p.G = p.G + a;

                if(p.B + a > 255)
                    p.B = 255;
                else if(p.B + a < 0)
                    p.B = 0;
                else
                    p.B = p.B + a;

                putPix(x,y,output,p);
            }
        }
    }
    }

    return output;
}

BitMap gray(BitMap bm)
{
    BitMap output = makeCopy(bm);
    float tmp;
    int x, y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);
            tmp = 0.21f * p.R + 0.72f * p.G + 0.07f * p.B;
            p.R = tmp;
            p.G = tmp;
            p.B = tmp;
            p.L = tmp;
            putPix(x,y,output,p);
        }
    }
    return output;
}

BitMap add(BitMap bm1, BitMap bm2)
{
    BitMap output = makeCopy(bm1);
    int x,y;
    for(y = bm1.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm1.info.Width; x++)
        {
            Pix pix1 = getPix(x,y,bm1);
            Pix pix2 = getPix(x,y,bm2);
            Pix outPix;
            outPix.R = align(pix1.R + pix2.R);
            outPix.G = align(pix1.G + pix2.G);
            outPix.B = align(pix1.B + pix2.B);
            outPix.L = align(pix1.L + pix2.L);
            putPix(x, y, output, outPix);
        }
    }
    return output;
}


int max(BitMap bm)
{
    int max = 0;
    int x,y;
    switch(bm.info.BitsPerPixel)
    {
    case 8:

        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x,y,bm);
                if(p.L > max)
                    max = p.L;
            }
        }
        return max;
        break;

    case 24:

        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x,y,bm);
                if(p.R > max)
                    max = p.R;
                if(p.G > max)
                    max = p.G;
                if(p.B > max)
                    max = p.B;
            }
        }
        return max;
        break;

    }

    return 0;
}

int* histogramR(BitMap bm)
{
    int* histogram = (int*)malloc(sizeof(int)*256);
    for (int i = 0; i < 256; i++) histogram[i] = 0;
    for (unsigned y = 0; y < bm.info.Height; y++)
    {
        for (unsigned x = 0; x < bm.info.Width; x++)
        {
            histogram[getPix(x, y, bm).R]++;
        }
    }
    return histogram;
}

int* histogramG(BitMap bm)
{
    int* histogram = (int*)malloc(sizeof(int)*256);
    for (int i = 0; i < 256; i++) histogram[i] = 0;
    for (unsigned y = 0; y < bm.info.Height; y++)
    {
        for (unsigned x = 0; x < bm.info.Width; x++)
        {
            histogram[getPix(x, y, bm).G]++;
        }
    }
    return histogram;
}

int* histogramB(BitMap bm)
{
    int* histogram = (int*)malloc(sizeof(int)*256);
    for (int i = 0; i < 256; i++) histogram[i] = 0;
    for (unsigned y = 0; y < bm.info.Height; y++)
    {
        for (unsigned x = 0; x < bm.info.Width; x++)
        {
            histogram[getPix(x, y, bm).B]++;
        }
    }
    return histogram;
}

BitMap multiplyConst(BitMap bm, float a)
{
    BitMap output = makeCopy(bm);
    int x,y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x, y, bm);
            p.R = align(p.R * a);
            p.G = align(p.G * a);
            p.B = align(p.B * a);
            putPix(x, y, output, p);
        }
    }
    return output;
}

BitMap negative(BitMap bm)
{
    BitMap output = makeCopy(bm);
    int x,y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);
            p.R =  255 - p.R;
            p.G =  255 - p.G;
            p.B =  255 - p.B;
            putPix(x,y,output,p);
        }
    }
    return output;
}

BitMap threshold(BitMap bm, int level)
{
    BitMap output = gray(bm);
    int x,y;
    for(y = output.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < output.info.Width; x++)
        {
            Pix p = getPix(x,y,output);
            p.R = (p.R > level) ? 255 : 0;
            p.G = (p.G > level) ? 255 : 0;
            p.B = (p.B > level) ? 255 : 0;
            putPix(x, y, output, p);
        }
    }
    return output;
}



BitMap fixLightness(BitMap bm)
{
    BitMap output = makeCopy(bm);
    int x,y;
    int Rmin=255, Gmin=255, Bmin=255, Rmax=1,Gmax=1,Bmax =1;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);
            if(p.R > Rmax)
                Rmax = p.R;
            if(p.G > Gmax)
                Gmax = p.G;
            if(p.B > Bmax)
                Bmax = p.B;
            if(p.R < Rmin)
                Rmin = p.R;
            if(p.G < Gmin)
                Gmin = p.G;
            if(p.B < Bmin)
                Bmin = p.B;
        }
    }
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);
            float tmp =  ((255.0f/((float)Rmax-(float)Rmin))*((float)p.R - (float)Rmin));
            p.R = (int)tmp;
            tmp =  ((255.0f/((float)Gmax-(float)Gmin))*((float)p.G - (float)Gmin));
            p.G = (int)tmp;
            tmp =  ((255.0f/((float)Bmax-(float)Bmin))*((float)p.B - (float)Bmin));
            p.B = (int)tmp;
            putPix(x,y,output, p);
        }
    }
    return output;
}




BitMap max(BitMap bm1, BitMap bm2)
{
    BitMap output = makeCopy(bm1);
    int x,y;
    for(y = bm1.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm1.info.Width; x++)
        {
            Pix pix1 = getPix(x,y,bm1);
            Pix pix2 = getPix(x,y,bm2);
            Pix outPix;
            outPix.R = max(pix1.R, pix2.R);
            outPix.G = max(pix1.G, pix2.G);
            outPix.B = max(pix1.B, pix2.B);
            outPix.L = max(pix1.L, pix2.L);
            putPix(x,y,output, outPix);
        }
    }
    return output;
}

int min(BitMap bm)
{
    int min = 255;
    int x,y;
    switch(bm.info.BitsPerPixel)
    {
    case 8:
    {
        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x,y,bm);
                if(p.L < min)
                    min = p.L;
            }
        }
        return min;
        break;
    }
    case 24:
    {
        for(y = bm.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p = getPix(x,y,bm);
                if(p.R < min)
                    min = p.R;
                if(p.G < min)
                    min = p.G;
                if(p.B < min)
                    min = p.B;
            }
        }
        return min;
        break;
    }
    }

    return 0;
}

float mean(BitMap bm)
{
    int result=0;
    int x,y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);
            if(bm.info.BitsPerPixel == 24)
                result += p.R + p.G + p.B;
            else if(bm.info.BitsPerPixel == 8)
                result += p.L;
        }
    }
    if(bm.info.BitsPerPixel == 24)
        return result/(bm.info.Width*bm.info.Height*3);
    else if(bm.info.BitsPerPixel == 8)
        return result/(bm.info.Width*bm.info.Height);

    return 0;
}

int median(BitMap bm)
{
    int* histL = histogramL(bm);
    int* histR = histogramR(bm);
    int* histG = histogramG(bm);
    int* histB = histogramB(bm);
    int i =0;
    if(bm.info.BitsPerPixel == 24)
    {
        int result = (bm.info.Width*bm.info.Height*3)/2;
        int sum = 0;
        for(i=0; i<256; i++)
        {
            sum+=histR[i]+histG[i]+histB[i];
            if(sum > result)
                return i;
        }
    }
    else
    {
        int result = (bm.info.Width*bm.info.Height)/2;
        int sum = 0;
        for(i=0; i<256; i++)
        {
            sum+=histL[i];
            if(sum > result)
                return i;
        }
    }

    return 0.0f;
}


BitMap difference(BitMap bm1, BitMap bm2)
{
    BitMap output = makeCopy(bm1);
    int x,y;
    for(y = bm1.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm1.info.Width; x++)
        {
            Pix pix1 = getPix(x,y,bm1);
            Pix pix2 = getPix(x,y,bm2);
            Pix outPix;
            outPix.R = ((pix1.R - pix2.R) < 0) ? 0 : pix1.R - pix2.R;
            outPix.G = ((pix1.G - pix2.G) < 0) ? 0 : pix1.G - pix2.G;
            outPix.B = ((pix1.B - pix2.B) < 0) ? 0 : pix1.B - pix2.B;
            outPix.L = ((pix1.L - pix2.L) < 0) ? 0 : pix1.L - pix2.L;
            putPix(x,y,output, outPix);
        }
    }
    return output;
}

BitMap multiply(BitMap bm1, BitMap bm2)
{
    BitMap output = makeCopy(bm1);
    int x,y;
    for(y = bm1.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm1.info.Width; x++)
        {
            Pix pix1 = getPix(x,y,bm1);
            Pix pix2 = getPix(x,y,bm2);
            Pix outPix;
            outPix.R = align(pix1.R * pix2.R / 255);
            outPix.G = align(pix1.G * pix2.G / 255);
            outPix.B = align(pix1.B * pix2.B / 255);
            outPix.L = align(pix1.L * pix2.L / 255);
            putPix(x,y,output, outPix);
        }
    }
    return output;
}

BitMap min(BitMap bm1, BitMap bm2)
{
    BitMap output = makeCopy(bm1);
    int x,y;
    for(y = bm1.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm1.info.Width; x++)
        {
            Pix pix1 = getPix(x,y,bm1);
            Pix pix2 = getPix(x,y,bm2);
            Pix outPix;
            outPix.R = min(pix1.R, pix2.R);
            outPix.G = min(pix1.G, pix2.G);
            outPix.B = min(pix1.B, pix2.B);
            outPix.L = min(pix1.L, pix2.L);
            putPix(x,y,output, outPix);
        }
    }
    return output;
}
BitMap linearCombination(BitMap bm1, BitMap bm2, float weight1, float weight2)
{
    return add(multiplyConst(bm1, weight1), multiplyConst(bm2, weight2));
}


///////////////////////////////
// cw3

BitMap move(BitMap bm,int wl, int hl, int cycle)
{
    BitMap output = makeCopy(bm);
    int x,y;
    for(y = bm.info.Height-(hl+1); y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width-(wl+1); x++)
        {
            Pix p = getPix(x,y,bm);
            putPix(x+wl,y+hl,output, p);
        }
    }
    if(cycle == 0)
    {
        for(y=bm.info.Height-1; y>=0; y--)
        {
            for(x = 0; x < wl+1; x++)
            {
                Pix p;
                p.R = 0;
                p.G = 0;
                p.B = 0;
                p.L = 0;
                putPix(x,y,output,p);
            }
        }
        for(y=hl+1; y>=0; y--)
        {
            for(x = 0; x < bm.info.Width; x++)
            {
                Pix p;
                p.R = 0;
                p.G = 0;
                p.B = 0;
                p.L = 0;
                putPix(x,y,output,p);
            }
        }
    }
    else
    {

        for(y=bm.info.Height-1; y>=hl+1; y--)
        {
            for(x = wl-1; x >= 0; x--)
            {
                Pix p = getPix(bm.info.Width-(wl-x),y,bm);
                putPix(x,y,output,p);
            }
        }

        for(y=0; y<hl; y++)
        {
            for( x = 0; x < bm.info.Width-wl; x++)
            {
                Pix p = getPix(x,bm.info.Height-(hl-1)+y,bm);
                putPix(x+wl,y,output,p);
            }
        }

        for( y=0; y<hl; y++)
        {
            for( x =0; x <wl; x++)
            {
                Pix p = getPix(bm.info.Width-(wl-1)+x,bm.info.Height-(hl-1)+y,bm);
                putPix(x,y,output,p);
            }
        }
    }
    return output;
}
BitMap crop(BitMap bm, int w, int h, int x,int y)
{
    BitMap output = setSize(w,h,bm);
    int i, j;
    for(j = output.info.Height-1; j >= 0 ; j--)
    {
        for(i = 0; i < output.info.Width-1; i++)
        {
            Pix p = getPix(i+x,j+y,bm);
            putPix(i,j,output,p);
        }
    }
    return output;
}

BitMap resize(BitMap bm, int w, int h, const char* method)
{
    BitMap output = setSize(w,h,bm);
    float ratioY = (float)bm.info.Height/(float)output.info.Height;
    float ratioX = (float)bm.info.Width/(float)output.info.Width;
    int x,y;
    if(method == "nearest")
    {
        for(y = output.info.Height-1; y >= 0 ; y--)
        {
            for(x = 0; x < output.info.Width-1; x++)
            {
                Pix p = getPix((int)x*ratioX,(int)y*ratioY,bm);
                putPix(x,y,output,p);
            }
        }
    }
    else if(method == "bilinear")
    {
        for(y = 0; y <output.info.Height-1 ; y++)
        {
            for(x = 0; x < output.info.Width-1; x++)
            {
                float Xratio = x*ratioX;
                float Yratio = y*ratioY;
                float Fa1,Fa0,Fab;
                float a = fabs(Xratio - (int)(x*ratioX));
                float b = fabs(Yratio - (int)(y*ratioY));
                Pix result;
                if(Yratio == output.info.Height && Xratio == output.info.Width)
                {
                    Pix pix1 = getPix((int)Xratio-1,(int)Yratio-1,bm);
                    Pix pix2 = getPix((int)Xratio,(int)Yratio-1,bm);
                    Pix pix3 = getPix((int)Xratio-1,(int)Yratio,bm);
                    Pix pix4 = getPix(Xratio,Yratio,bm);
                    //R
                    Fa0 = (1-a)*pix1.R+a*pix2.R;
                    Fa1 = (1-a)*pix3.R+a*pix4.R;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.R = (Fab > 255) ? 255 : Fab;
                    //B
                    Fa0 = (1-a)*pix1.B+a*pix2.B;
                    Fa1 = (1-a)*pix3.B+a*pix4.B;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.B = (Fab > 255) ? 255 : Fab;
                    //G
                    Fa0 = (1-a)*pix1.G+a*pix2.G;
                    Fa1 = (1-a)*pix3.G+a*pix4.G;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.G = (Fab > 255) ? 255 : Fab;
                }
                else if(Yratio == output.info.Height)
                {
                    Pix pix1 = getPix(Xratio,Yratio-1,bm);
                    Pix pix2 = getPix(Xratio+1,Yratio-1,bm);
                    Pix pix3 = getPix(Xratio,Yratio,bm);
                    Pix pix4 = getPix(Xratio+1,Yratio,bm);
                    //R
                    Fa0 = (1-a)*pix1.R+a*pix2.R;
                    Fa1 = (1-a)*pix3.R+a*pix4.R;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.R = Fab;
                    //B
                    Fa0 = (1-a)*pix1.B+a*pix2.B;
                    Fa1 = (1-a)*pix3.B+a*pix4.B;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.B = Fab;
                    //G
                    Fa0 = (1-a)*pix1.G+a*pix2.G;
                    Fa1 = (1-a)*pix3.G+a*pix4.G;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.G = (Fab > 255) ? 255 : Fab;
                }
                else if(Xratio == output.info.Width)
                {
                    Pix pix1 = getPix(Xratio-1,Yratio,bm);
                    Pix pix2 = getPix(Xratio,Yratio,bm);
                    Pix pix3 = getPix(Xratio-1,Yratio+1,bm);
                    Pix pix4 = getPix(Xratio,Yratio+1,bm);
                    //R
                    Fa0 = (1-a)*pix1.R+a*pix2.R;
                    Fa1 = (1-a)*pix3.R+a*pix4.R;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.R = Fab;
                    //B
                    Fa0 = (1-a)*pix1.B+a*pix2.B;
                    Fa1 = (1-a)*pix3.B+a*pix4.B;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.B = Fab;
                    //G
                    Fa0 = (1-a)*pix1.G+a*pix2.G;
                    Fa1 = (1-a)*pix3.G+a*pix4.G;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.G = Fab;
                }
                else
                {
                    Pix pix1 = getPix(Xratio,Yratio,bm);
                    Pix pix2 = getPix(Xratio+1,Yratio,bm);
                    Pix pix3 = getPix(Xratio,Yratio+1,bm);
                    Pix pix4 = getPix(Xratio+1,Yratio+1,bm);
                    //R
                    Fa0 = (1-a)*pix1.R+a*pix2.R;
                    Fa1 = (1-a)*pix3.R+a*pix4.R;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.R = Fab;
                    //B
                    Fa0 = (1-a)*pix1.B+a*pix2.B;
                    Fa1 = (1-a)*pix3.B+a*pix4.B;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.B = Fab;
                    //G
                    Fa0 = (1-a)*pix1.G+a*pix2.G;
                    Fa1 = (1-a)*pix3.G+a*pix4.G;
                    Fab = (1-b)*Fa0+b*Fa1;
                    result.G = Fab;
                }
                putPix(x,y,output, result);
            }
        }
    }
    return output;
}

BitMap rotate(BitMap bm, int deg)
{
    BitMap output;
    int x,y;
    switch(deg)
    {
    case 90:
    {
        output = setSize(bm.info.Height,bm.info.Width,bm);
        for(y = 0; y < output.info.Height-1 ; y++)
        {
            for( x = 0; x < output.info.Width-1; x++)
            {
                Pix p = getPix(y,bm.info.Height-1-x,bm);
                putPix(x,y,output,p);
            }
        }
        break;
    }
    case 180:
    {
        output = makeCopy(bm);
        for(y = 0; y < output.info.Height; y++)
        {
            for(x = 0; x < output.info.Width; x++)
            {
                Pix p = getPix(bm.info.Width-x,bm.info.Height-y,bm);
                putPix(x,y,output,p);
            }
        }
        break;
    }
    case 270:
    {
        output = setSize(bm.info.Height,bm.info.Width,bm);
        for(y = 0; y < output.info.Height-1 ; y++)
        {
            for( x = 0; x < output.info.Width-1; x++)
            {
                Pix p = getPix(bm.info.Width-1-y,x,bm);
                putPix(x,y,output,p);
            }
        }
        break;
    }
    }
    return output;
}

BitMap flip(BitMap bm, int direct)
{
    BitMap output;
    int x,y;
    if(direct == 1) //odbicie w poziomie
    {
        output = makeCopy(bm);
        for(y = 0; y < bm.info.Height-1 ; y++)
        {
            for(x = 0; x < bm.info.Width-1; x++)
            {
                Pix p = getPix(x,y,bm);
                putPix(x,bm.info.Height-1-y,output,p);
            }
        }
    }
    else
    {
        output = makeCopy(bm);
        for(y = 0; y < output.info.Height-1 ; y++)
        {
            for(x = 0; x < output.info.Width-1; x++)
            {
                Pix p = getPix(bm.info.Width-1-x,y,bm);
                putPix(x,y,output,p);
            }
        }
    }
    return output;
}


///////////////////////////////
// cw4

BitMap filter(BitMap bm, float** matrix, int w, int h)
{
    BitMap output = makeCopy(bm);
    int n = 0;
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			n += matrix[i][j];
		}
	}

	if (n == 0) n = 1;

    for(int y = h; y < bm.info.Height - h ; y++)
    {
        for(int x = w; x < bm.info.Width - w; x++)
        {
            int R = 0;
            int G = 0;
            int B = 0;

            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w; i++)
                {
                    Pix p = getPix(x + i, y + j, bm);

                    R += matrix[i][j] * p.R;
                    G += matrix[i][j] * p.G;
                    B += matrix[i][j] * p.B;
                }
            }

            R /= n;
            G /= n;
            B /= n;

            Pix outPix;
            outPix.R = align(R);
            outPix.G = align(G);
            outPix.B = align(B);

            putPix(x, y, output, outPix);
        }
    }

    return output;
}

BitMap dilate(BitMap bm, int** strEl, int w, int h)
{
    BitMap output = makeCopy(bm);
    int x,y,i,j;
    for(y = h/2; y < bm.info.Height-(h-1); y++)
    {
        for(x = w/2; x < bm.info.Width-(w-1); x++)
        {
            int max=0;
            for(j=0; j<h; j++)
            {
                for(i = 0; i<w; i++)
                {
                    if(strEl[i][j] == 1)
                    {
                        Pix p = getPix(x+(i-w/2),y+(j-h/2),bm);
                        if(p.L > max)
                            max =p.L;
                    }
                }
            }
            Pix result = getPix(x,y,output);
            result.L = max;
            putPix(x,y, output, result);
        }
    }
    return output;
}
BitMap erode(BitMap bm,int ** strEl,int w,int h)
{
    BitMap output = makeCopy(bm);
    int x,y,i,j;
    for(y = 1; y < bm.info.Height-2 ; y++)
    {
        for(x = 1; x < bm.info.Width-2; x++)
        {
            int min;
            if(bm.info.BitsPerPixel == 8)
            {
                min=255;
            }
            else if(bm.info.BitsPerPixel == 1)
            {
                min = 1;
            }
            for( i=0; i<w; i++)
            {
                for(j=0; j<h; j++)
                {
                    if(strEl[i][j] == 1)
                    {
                        Pix p = getPix(x+(i-w/2),y+(j-h/2),bm);
                        if(p.L < min)
                            min =p.L;
                    }
                }
            }
            Pix result = getPix(x,y,output);
            result.L = min;
            putPix(x,y,output, result);
        }
    }
    return output;
}
BitMap openMor(BitMap bm, int** strEl,int w, int h)
{
    return dilate(erode(bm,strEl,w,h),strEl,w,h);
}
BitMap closeMor(BitMap bm, int** strEl,int w, int h)
{
    return erode(dilate(bm,strEl,w,h),strEl,w,h);
}

BitMap whiteTH(BitMap bm, int ** strEl, int w, int h)
{
	return difference(bm, openMor(bm, strEl, w, h));	
}

BitMap blackTH(BitMap bm, int ** strEl, int w, int h)
{
	return difference(closeMor(bm, strEl, w, h), bm);	
}


///////////////////////////////
// cw5

BitMap addNoiseRec(BitMap bm, int level, float density)
{
	BitMap output = makeCopy(bm);
    int x,y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);

			if (density > rand01())
			{
				p.R = align(p.R + randRange(-level, level));
				p.G = align(p.G + randRange(-level, level));
				p.B = align(p.B + randRange(-level, level));
			}

			putPix(x, y, output, p);

		}
	}

	return output;

}

BitMap filterMean(BitMap bm)
{
	float m[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}
BitMap filterGauss(BitMap bm)
{
	float m[] = {1, 2, 1, 2, 4, 2, 1, 2, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}

BitMap filterMeanA(BitMap bm)
{
	float m[] = {1, 1, 1, 1, 2, 1, 1, 1, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}

BitMap filterMeanB(BitMap bm)
{
	float m[] = {1, 1, 1, 1, 0, 1, 1, 1, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}

BitMap filterRoberts(BitMap bm)
{
	float m[] = {0, 0, 0, -1, 0, 0, 0, 1, 0};
	return filter(bm, makeMask3x3(m), 3, 3);
}
BitMap filterPrewitt(BitMap bm)
{
	float m[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}

BitMap filterSobel(BitMap bm)
{
	float m[] = {-1, -2, -1, -1, 0, 0, 1, 2, 1};
	return filter(bm, makeMask3x3(m), 3, 3);
}

BitMap filterLaplace(BitMap bm)
{
	float m[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
	return filter(bm, makeMask3x3(m), 3, 3);
}


BitMap filterMedian(BitMap bm, int w, int h)
{
	BitMap output = makeCopy(bm);

	int n = w * h;

	int* red = new int[n];
	int* green = new int[n];
	int* blue = new int[n];

    for(int y = h; y < bm.info.Height - h ; y++)
    {
        for(int x = w; x < bm.info.Width - w; x++)
        {
            
            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w; i++)
                {
                    Pix p = getPix(x + i, y + j, bm);

                    red[3 * i + j] = p.R;
					green[3 * i + j] = p.G;
					blue[3 * i + j] = p.B;
                }
            }

			qsort(red, n, sizeof(int), cmp);
			qsort(green, n, sizeof(int), cmp);
			qsort(blue, n, sizeof(int), cmp);

            Pix outPix;
            outPix.R = red[n/2];
            outPix.G = green[n/2];
            outPix.B = blue[n/2];

            putPix(x, y, output, outPix);
        }
    }

    return output;
}



BitMap addNoiseNor(BitMap bm, int mean, float sigma, float density)
{
	BitMap output = makeCopy(bm);
    int x,y;
    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);

			if (density > rand01())
			{
				p.R = align(p.R + normalDistGenerator(mean, sigma));
				p.G = align(p.G + normalDistGenerator(mean, sigma));
				p.B = align(p.B + normalDistGenerator(mean, sigma));
			}

			putPix(x, y, output, p);

		}
	}

	return output;
}

BitMap addNoiseSP(BitMap bm, float density)
{
	BitMap output = makeCopy(bm);
    int x,y;

	int r = rand() % 3;

    for(y = bm.info.Height-1; y >= 0 ; y--)
    {
        for(x = 0; x < bm.info.Width; x++)
        {
            Pix p = getPix(x,y,bm);

			if (density > rand01())
			{
				if (rand() % 2 == 0)
				{
					if (r == 0) p.R = 0;
					else if (r == 1) p.G = 0;
					else if (r == 2) p.B = 0;
				}
				else
				{
					if (r == 0) p.R = 255;
					else if (r == 1) p.G = 255;
					else if (r == 2) p.B = 255;
				}
			}

			putPix(x, y, output, p);

		}
	}
	return output;
}



int main(int argc, char* argv[])
{
    BitMap bm1 = readBitMap(fopen("wally1.bmp", "rb"));
    BitMap bm2 = readBitMap(fopen("wally2.bmp", "rb"));
    BitMap bm3 = readBitMap(fopen("SS_filterLaplace.bmp", "rb"));
 	//writeBitMap(fopen("_flip.bmp", "wb"), flip(bm2, 0));
	//writeBitMap(fopen("_difference.bmp", "wb"), difference(bm1, bm3));
	
	
	writeBitMap(fopen("SS_negativeee.bmp", "wb"), negative(bm3));	
	writeBitMap(fopen("SS_treshold.bmp", "wb"), negative(threshold(bm3,50)));
   	/*
    // cw2
    writeBitMap(fopen("_addConst.bmp", "wb"), addConst(bm1, 100));
    writeBitMap(fopen("_gray.bmp", "wb"), gray(bm1));
    writeBitMap(fopen("_add.bmp", "wb"), add(bm1, bm2));
    writeBitMap(fopen("_max.bmp", "wb"), max(bm1, bm2));
    writeBitMap(fopen("_multiplyConst.bmp", "wb"), multiplyConst(bm1, 1.2f));
    writeBitMap(fopen("_negative.bmp", "wb"), negative(bm1));
    writeBitMap(fopen("_threshold.bmp", "wb"), threshold(bm1, 100));
    writeBitMap(fopen("_fixLightness.bmp", "wb"), fixLightness(bm1));

    //printf("max = %d\n", max(bm1));
    //printf("min = %d\n", min(bm1));
    //printf("mean = %f\n", mean(bm1));
    //printf("median = %d\n", median(bm1));

    writeBitMap(fopen("_difference.bmp", "wb"), difference(bm1, bm2));
    writeBitMap(fopen("_multiply.bmp", "wb"), multiply(bm1, bm2));
    writeBitMap(fopen("_min.bmp", "wb"), min(bm1, bm2));
    writeBitMap(fopen("_linearCombination.bmp", "wb"), linearCombination(bm1, bm2, 0.5, 2.0));

    // cw3
    writeBitMap(fopen("_move.bmp", "wb"), move(bm1, 20, 100, 1));
    writeBitMap(fopen("_crop.bmp", "wb"), crop(bm1, 100, 100, 50, 50));
    writeBitMap(fopen("_resize_bilinear.bmp", "wb"), resize(bm1, 200, 200, "bilinear"));
    writeBitMap(fopen("_resize_nearest.bmp", "wb"), resize(bm1, 200, 200, "nearest"));
    writeBitMap(fopen("_rotate.bmp", "wb"), rotate(bm1, 90));
    writeBitMap(fopen("_flip.bmp", "wb"), flip(bm1, 1));

    // cw4
    float m[] = {1, 1, 1, 1, 8, 1, 1, 1, 1};
	float** mask = makeMask3x3(m);
    writeBitMap(fopen("_filter.bmp", "wb"), filter(bm1, mask, 3, 3));

    int s[] = {1, 0, -1, 0, -1, 1, -1, 1, 0};
	int** strEl = makeStrEl3x3(s);
    writeBitMap(fopen("_dilate.bmp", "wb"), dilate(bm1, strEl, 3, 3));
    writeBitMap(fopen("_erode.bmp", "wb"), erode(bm1, strEl, 3, 3));
	writeBitMap(fopen("_openMor.bmp", "wb"), openMor(bm1, strEl, 3, 3));
    writeBitMap(fopen("_closeMor.bmp", "wb"), closeMor(bm1, strEl, 3, 3));

	// cw5
	writeBitMap(fopen("_addNoiseRec.bmp", "wb"), addNoiseRec(bm1, 50, 0.3f));

	writeBitMap(fopen("_filterMean.bmp", "wb"), filterMean(bm1));
	writeBitMap(fopen("_filterGauss.bmp", "wb"), filterGauss(bm1));
	writeBitMap(fopen("_filterMeanA.bmp", "wb"), filterMeanA(bm1));
	writeBitMap(fopen("_filterMeanB.bmp", "wb"), filterMeanB(bm1));

	writeBitMap(fopen("_filterRoberts.bmp", "wb"), filterRoberts(bm1));
	writeBitMap(fopen("_filterPrewitt.bmp", "wb"), filterPrewitt(bm1));
	writeBitMap(fopen("_filterSobel.bmp", "wb"), filterSobel(bm1));
	writeBitMap(fopen("_filterLaplace.bmp", "wb"), filterLaplace(bm1));

	writeBitMap(fopen("_filterMedian.bmp", "wb"), filterMedian(bm1, 3, 3));

	writeBitMap(fopen("_addNoiseNor.bmp", "wb"), addNoiseNor(bm1, 10, 5.0f, 0.3f));
	writeBitMap(fopen("_addNoiseSP.bmp", "wb"), addNoiseSP(bm1, 0.3f));

	*/

    return 0;
}

