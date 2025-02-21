#pragma once

#include <omp.h>
#include <math.h>
#include <stdio.h>

#define TileLevel16_bit 4

typedef unsigned char TileIndex;
typedef unsigned int MatIndex;
typedef double MatValue;
typedef unsigned short ushort;
typedef ushort TileBitmap;

#define TileIndex_size 2  // bits (uint 2)
#define MatIndex_size 32  // bits (uint 32)
#define MatValue_size 64  // bits (double 64)

#define ushort_popcount(x) __builtin_popcount(x)
#define TileBitmap_lowbit(bitmap) ((bitmap) & (-(bitmap)))
#define TileBitmap_popcount(x) __builtin_popcount(x)
#define TileBitmap_ctz(bitmap) __builtin_ctz(bitmap)

#define BSR_BlockSz_bit 4
#define TileBlockSz_bit 2
#define BSR_BlockSz 16
#define TileBlockSz 4
#define TileBlockSz_negetive_one 0xf
#define negetive_one 0xff

// #define max(x, y) ((x) > (y) ? (x) : (y))
// #define min(x, y) ((x) > (y) ? (y) : (x))

#define multiply_5(x) ((x << 2) + x)

typedef enum {
    False,
    True
} Bool;

Bool file_exists(const char* filename)
{
    FILE*f = fopen(filename, "r");
    if (f == NULL)
        return False;
    fclose(f);
    return True;
}