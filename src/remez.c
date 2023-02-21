/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *-------------------------------------------------
 *  Copyright (c) 1995,1998  Jake Janovetz <janovetz@uiuc.edu>
 *  Copyright (c) 2023 Rick Kern <kernrj@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *  Sep 1999 - Paul Kienzle (pkienzle@cs.indiana.edu)
 *      Modified for use in octave as a replacement for the matlab function
 *      remez.mex.  In particular, magnitude responses are required for all
 *      band edges rather than one per band, griddensity is a parameter,
 *      and errors are returned rather than printed directly.
 *  Mar 2000 - Kai Habel (kahacjde@linux.zrz.tu-berlin.de)
 *      Change: ColumnVector x=arg(i).vector_value();
 *      to: ColumnVector x(arg(i).vector_value());
 *  Jan 2023 - Rick Kern (kernrj@gmail.com)
 *      Avoid writing past the end of buffers
 *      Return an error status when an error occurs
 *      Fix memory leak
 *
 *  There appear to be some problems with the routine Search. See comments
 *  therein [search for PAK:].  I haven't looked closely at the rest
 *  of the code---it may also have some problems.
 *************************************************************************/

#include "remez/remez.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
  NEGATIVE,
  POSITIVE,
} Symmetry;

#define Pi 3.1415926535897932
#define PiTimes2 6.2831853071795865

#define RET_ON_ERR(remezCmd__)                 \
  do {                                         \
    RemezStatus retOnErrStatus__ = remezCmd__; \
    if (retOnErrStatus__ != RemezSuccess) {    \
      return retOnErrStatus__;                 \
    }                                          \
  } while (0)

const char* remezStatusToString(RemezStatus status) {
  static char unknownBuffer[64];

  switch (status) {
    case RemezSuccess:
      return "Success";
    case RemezInvalidParameter:
      return "Invalid Parameter";
    case RemezTooManyExtremalFrequencies:
      return "Too many extremal frequencies";
    case RemezTooFewExtremalFrequencies:
      return "Too few extremal frequencies";
    case RemezDidNotConverge:
      return "Did not converge";
    default:
      snprintf(
          unknownBuffer,
          sizeof(unknownBuffer) - 1,
          "Unknown (%d)",
          status);
      unknownBuffer[sizeof(unknownBuffer) - 1] = 0;

      return unknownBuffer;
  }
}

/*******************
 * CreateDenseGrid
 *=================
 * Creates the dense grid of frequencies from the specified bands.
 * Also creates the Desired Frequency Response function (D[]) and
 * the Weight function (W[]) on that dense grid
 *
 *
 * INPUT:
 * ------
 * size_t      r        - 1/2 the number of filter coefficients
 * size_t      numtaps  - Number of taps in the resulting filter
 * double   bands[]  - User-specified band edges [2*numband]
 * double   des[]    - Desired response per band [2*numband]
 * double   weight[] - Weight per band [numband]
 * size_t      symmetry - Symmetry of filter - used for grid check
 * size_t      griddensity
 *
 * OUTPUT:
 * -------
 * double Grid[]     - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]        - Desired response on the dense grid [gridsize]
 * double W[]        - Weight function on the dense grid [gridsize]
 *******************/
static RemezStatus CreateDenseGrid(
    size_t r,
    size_t numtaps,
    const RemezBand* bands,
    size_t bandCount,
    double* Grid,
    size_t gridSize,
    double* D,
    size_t dSize,
    double* W,
    size_t wSize,
    Symmetry symmetry,
    size_t griddensity) {
  if (gridSize != dSize || gridSize != wSize) {
    fprintf(stderr, "Grid buffer size mismatch in CreateDenseGrid()\n");
    return RemezInvalidParameter;
  }

  const double delf = 0.5 / (double)(griddensity * r);

  /*
   * For differentiator, hilbert,
   *   symmetry is odd and Grid[0] = max(delf, bands[0])
   */
  const double grid0 =
      ((symmetry == NEGATIVE) && (delf > bands[0].lowFrequency))
          ? delf
          : bands[0].lowFrequency;

  for (size_t j = 0, bandIndex = 0; bandIndex < bandCount; bandIndex++) {
    const RemezBand* band = &bands[bandIndex];
    double lowf = (bandIndex == 0 ? grid0 : band->lowFrequency);
    const double highf = band->highFrequency;
    const size_t k = lrint((highf - lowf) / delf);

    for (size_t i = 0; i < k; i++) {
      if (j < gridSize) {
        D[j] =
            band->lowFrequencyResponse
            + (double)i
                  * (band->highFrequencyResponse - band->lowFrequencyResponse)
                  / (double)(k - 1);
        W[j] = band->weight;
        Grid[j] = lowf;
        lowf += delf;
        j++;
      }
    }
    if (j - 1 < gridSize) {
      Grid[j - 1] = highf;
    }

    if (j > gridSize) {
      fprintf(stderr, "Grid size is too small. At least [%zd] is needed.\n", k);
      return RemezTooFewExtremalFrequencies;
    }
  }

  /*
   * Similar to above, if odd symmetry, last grid point can't be .5
   *  - but, if there are even taps, leave the last grid point at .5
   */
  if ((symmetry == NEGATIVE) && (Grid[gridSize - 1] > (0.5 - delf))
      && (numtaps % 2)) {
    Grid[gridSize - 1] = 0.5 - delf;
  }

  return RemezSuccess;
}

/********************
 * InitialGuess
 *==============
 * Places Extremal Frequencies evenly throughout the dense grid.
 *
 *
 * INPUT:
 * ------
 * size_t r        - 1/2 the number of filter coefficients
 * size_t gridsize - Number of elements in the dense frequency grid
 *
 * OUTPUT:
 * -------
 * size_t Ext[]    - Extremal indexes to dense frequency grid [r+1]
 ********************/

static RemezStatus InitialGuess(
    size_t r,
    size_t* Ext,
    size_t extSize,
    size_t gridsize) {
  if (gridsize == 0) {
    fprintf(stderr, "Grid size is too small\n");
    return RemezInvalidParameter;
  }

  if (r >= extSize) {
    fprintf(stderr, "r [%zu] is out of range. Max [%zu].\n", r, extSize - 1);
    return RemezInvalidParameter;
  }

  for (size_t i = 0; i < extSize; i++) {
    Ext[i] = i * (gridsize - 1) / r;
  }

  return RemezSuccess;
}

/***********************
 * CalcParms
 *===========
 *
 *
 * INPUT:
 * ------
 * size_t    Ext[]  - Extremal indexes to dense frequency grid [r+1]
 * double Grid[] - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double ad[]   - 'b' in Oppenheim & Schafer [r+1]
 * double x[]    - [r+1]
 * double y[]    - 'C' in Oppenheim & Schafer [r+1]
 ***********************/
static RemezStatus CalcParms(
    size_t r,
    const size_t* Ext,
    size_t extSize,
    const double* Grid,
    size_t gridSize,
    const double* D,
    size_t dSize,
    const double* W,
    size_t wSize,
    double* ad,
    size_t adSize,
    double* x,
    size_t xSize,
    double* y,
    size_t ySize) {
  if (extSize != adSize || extSize != xSize || extSize != ySize) {
    fprintf(stderr, "Ext buffer size mismatch\n");
    return RemezInvalidParameter;
  }

  if (gridSize != dSize || gridSize != wSize) {
    fprintf(stderr, "Grid buffer size mismatch\n");
    return RemezInvalidParameter;
  }

  if (r >= extSize) {
    fprintf(stderr, "r [%zu] is out of range. Max [%zu].\n", r, extSize - 1);
    return RemezInvalidParameter;
  }

  /*
   * Find x[]
   */
  for (size_t i = 0; i <= r; i++) {
    x[i] = cos(PiTimes2 * Grid[Ext[i]]);
  }

  /*
   * Calculate ad[]  - Oppenheim & Schafer eq 7.132
   */
  size_t ld = (r - 1) / 15 + 1; /* Skips around to avoid round errors */
  for (size_t i = 0; i <= r; i++) {
    double denom = 1.0;
    double xi = x[i];
    for (size_t j = 0; j < ld; j++) {
      for (size_t k = j; k <= r; k += ld)
        if (k != i) denom *= 2.0 * (xi - x[k]);
    }
    if (fabs(denom) < 0.00001) denom = 0.00001;
    ad[i] = 1.0 / denom;
  }

  /*
   * Calculate delta  - Oppenheim & Schafer eq 7.131
   */
  double numer = 0.0;
  double denom = 0.0;
  int32_t sign = 1;
  for (size_t i = 0; i <= r; i++) {
    numer += ad[i] * D[Ext[i]];
    denom += sign * ad[i] / W[Ext[i]];
    sign = -sign;
  }
  const double delta = numer / denom;
  sign = 1;

  /*
   * Calculate y[]  - Oppenheim & Schafer eq 7.133b
   */
  for (size_t i = 0; i <= r; i++) {
    y[i] = D[Ext[i]] - sign * delta / W[Ext[i]];
    sign = -sign;
  }

  return RemezSuccess;
}

/*********************
 * ComputeA
 *==========
 * Using values calculated in CalcParms, ComputeA calculates the
 * actual filter response at a given frequency (freq).  Uses
 * eq 7.133a from Oppenheim & Schafer.
 *
 *
 * INPUT:
 * ------
 * double freq - Frequency (0 to 0.5) at which to calculate A
 * size_t    r    - 1/2 the number of filter coefficients
 * double ad[] - 'b' in Oppenheim & Schafer [r+1]
 * double x[]  - [r+1]
 * double y[]  - 'C' in Oppenheim & Schafer [r+1]
 *
 * OUTPUT:
 * -------
 * Returns double value of A[freq]
 *********************/

static RemezStatus ComputeA(
    double freq,
    size_t r,
    const double* ad,
    size_t adSize,
    const double* x,
    size_t xSize,
    const double* y,
    size_t ySize,
    double* outA) {
  if (adSize != xSize || adSize != ySize || outA == NULL) {
    fprintf(stderr, "ad buffer size mismatch in ComputeA()\n");
    return RemezInvalidParameter;
  }

  if (r >= adSize) {
    fprintf(stderr, "r [%zu] is out of range. Max [%zu].\n", r, adSize - 1);
    return RemezInvalidParameter;
  }

  double numer = 0.0;
  double denom = 0.0;
  const double xc = cos(PiTimes2 * freq);
  for (size_t i = 0; i <= r; i++) {
    double c = xc - x[i];
    if (fabs(c) < 1.0e-7) {
      numer = y[i];
      denom = 1;
      break;
    }
    c = ad[i] / c;
    denom += c;
    numer += c * y[i];
  }

  *outA = numer / denom;

  return RemezSuccess;
}

/************************
 * CalcError
 *===========
 * Calculates the Error function from the desired frequency response
 * on the dense grid (D[]), the weight function on the dense grid (W[]),
 * and the present response calculation (A[])
 *
 *
 * INPUT:
 * ------
 * size_t    r      - 1/2 the number of filter coefficients
 * double ad[]   - [r+1]
 * double x[]    - [r+1]
 * double y[]    - [r+1]
 * double Grid[] - Frequencies on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the desnse grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double E[]    - Error function on dense grid [gridsize]
 ************************/
static RemezStatus CalcError(
    size_t r,
    const double* ad,
    size_t adSize,
    const double* x,
    size_t xSize,
    const double* y,
    size_t ySize,
    const double* Grid,
    size_t gridSize,
    const double* D,
    size_t dSize,
    const double* W,
    size_t wSize,
    double* E,
    size_t eSize) {
  if (adSize != xSize || adSize != ySize) {
    fprintf(stderr, "ad buffer size mismatch in CalcError()\n");
    return RemezInvalidParameter;
  }

  if (gridSize != dSize || gridSize != wSize || gridSize != eSize) {
    fprintf(stderr, "Grid buffer size mismatch in CalcError()\n");
    return RemezInvalidParameter;
  }

  for (size_t i = 0; i < gridSize; i++) {
    double A = 0.0;
    RemezStatus status =
        ComputeA(Grid[i], r, ad, adSize, x, xSize, y, ySize, &A);
    if (status != RemezSuccess) {
      return status;
    }

    E[i] = W[i] * (D[i] - A);
  }

  return RemezSuccess;
}

/************************
 * Search
 *========
 * Searches for the maxima/minima of the error curve.  If more than
 * r+1 extrema are found, it uses the following heuristic (thanks
 * Chris Hanson):
 * 1) Adjacent non-alternating extrema deleted first.
 * 2) If there are more than one excess extrema, delete the
 *    one with the smallest error.  This will create a non-alternation
 *    condition that is fixed by 1).
 * 3) If there is exactly one excess extremum, delete the smaller
 *    of the first/last extremum
 *
 *
 * INPUT:
 * ------
 * size_t    r        - 1/2 the number of filter coefficients
 * size_t    Ext[]    - Indexes to Grid[] of extremal frequencies [r+1]
 * size_t    gridsize - Number of elements in the dense frequency grid
 * double E[]      - Array of error values.  [gridsize]
 * OUTPUT:
 * -------
 * size_t    Ext[]    - New indexes to extremal frequencies [r+1]
 ************************/
static RemezStatus Search(
    size_t r,
    size_t* Ext,
    size_t extSize,
    double* E,
    size_t eSize) {
  RemezStatus status = RemezSuccess;
  size_t lastGridIndex, extra;
  size_t k = 0;
  const size_t gridsize = eSize;
  // size_t up, alt;
  const size_t foundExtSize = gridsize;  // originally 2r
  size_t* foundExt =
      calloc(foundExtSize, sizeof(size_t)); /* Array of found extremals */
  if (foundExt == NULL) {
    status = RemezOutOfMemory;
    goto end;
  }

  /*
   * Allocate enough space for found extremals.
   */
  k = 0;

  /*
   * Check for extremum at 0.
   */
  if (((E[0] > 0.0) && (E[0] > E[1])) || ((E[0] < 0.0) && (E[0] < E[1]))) {
    foundExt[k++] = 0;
  }

  /*
   * Check for extrema inside dense grid
   */
  for (size_t i = 1; i < gridsize - 1; i++) {
    if (((E[i] >= E[i - 1]) && (E[i] > E[i + 1]) && (E[i] > 0.0))
        || ((E[i] <= E[i - 1]) && (E[i] < E[i + 1]) && (E[i] < 0.0))) {
      // PAK: we sometimes get too many extremal frequencies
      if (k >= foundExtSize) {
        status = RemezTooManyExtremalFrequencies;
        goto end;
      }
      foundExt[k++] = i;
    }
  }

  /*
   * Check for extremum at 0.5
   */
  lastGridIndex = gridsize - 1;
  if (((E[lastGridIndex] > 0.0) && (E[lastGridIndex] > E[lastGridIndex - 1]))
      || ((E[lastGridIndex] < 0.0)
          && (E[lastGridIndex] < E[lastGridIndex - 1]))) {
    if (k >= foundExtSize) {
      status = RemezTooManyExtremalFrequencies;
      goto end;
    }
    foundExt[k++] = lastGridIndex;
  }

  // PAK: we sometimes get not enough extremal frequencies
  if (k < r + 1) {
    fprintf(
        stderr,
        "Too few extremal frequencies: k (%zd) < r (%zd) + 1\n",
        k,
        r);
    status = RemezTooFewExtremalFrequencies;
    goto end;
  }

  /*
   * Remove extra extremals
   */
  extra = k - (r + 1);
  //   assert(extra >= 0);

  while (extra > 0) {
    size_t up, l, alt;

    if (E[foundExt[0]] > 0.0) {
      up = 1; /* first one is a maxima */
    } else {
      up = 0; /* first one is a minima */
    }

    l = 0;
    alt = 1;
    for (size_t j = 1; j < k; j++) {
      if (fabs(E[foundExt[j]]) < fabs(E[foundExt[l]])) {
        l = j; /* new smallest error. */
      }

      if (up && (E[foundExt[j]] < 0.0)) {
        up = 0; /* switch to a minima */
      } else if ((!up) && (E[foundExt[j]] > 0.0)) {
        up = 1; /* switch to a maxima */
      } else {
        alt = 0;
        // PAK: break now and you will delete the smallest overall
        // extremal.  If you want to delete the smallest of the
        // pair of non-alternating extremals, then you must do:
        //
        // if (fabs(E[foundExt[j]]) < fabs(E[foundExt[j-1]])) l=j;
        // else l=j-1;
        break; /* Ooops, found two non-alternating */
      }        /* extrema.  Delete smallest of them */
    }          /* if the loop finishes, all extrema are alternating */

    /*
     * If there's only one extremal and all are alternating,
     * delete the smallest of the first/last extremals.
     */
    if ((alt) && (extra == 1)) {
      if (fabs(E[foundExt[k - 1]])
          < fabs(E[foundExt[0]])) /* Delete last extremal */
        l = k - 1;
      // PAK: changed from l = foundExt[k-1];
      else
        /* Delete first extremal */
        l = 0;
      // PAK: changed from l = foundExt[0];
    }

    for (size_t j = l; j < k - 1; j++) { /* Loop that does the deletion */
      foundExt[j] = foundExt[j + 1];
      //  assert(foundExt[j]<gridsize);
    }
    k--;
    extra--;
  }

  for (size_t i = 0; i <= r; i++) {
    //      assert(foundExt[i]<gridsize);
    Ext[i] = foundExt[i]; /* Copy found extremals to Ext[] */
  }

end:
  if (foundExt != NULL) {
    free(foundExt);
  }

  return status;
}

/*********************
 * FreqSample
 *============
 * Simple frequency sampling algorithm to determine the impulse
 * response outTaps[] from A's found in ComputeA
 *
 *
 * INPUT:
 * ------
 * double   A[]      - Sample points of desired response [N/2]
 * Symmetry symmetry - Symmetry of desired filter
 *
 * OUTPUT:
 * -------
 * double outTaps[] - Impulse Response of final filter [N]
 *********************/
static RemezStatus FreqSample(
    const double* A,
    size_t aLength,
    double* outTaps,
    size_t outTapsLength,
    Symmetry symmetry) {
  const size_t N = outTapsLength;
  double M = (N - 1.0) / 2.0;

  if (outTapsLength < N) {
    fprintf(stderr, "Taps buffer is too small\n");
    return RemezInternalError;
  }

  if (aLength <= (N / 2 - 1) || aLength <= floor(M)) {
    fprintf(stderr, "input taps length [%zu] is too small.\n", aLength);
    return RemezInternalError;
  }

  if (N == 0) {
    fprintf(stderr, "Size of A must be greater than 0.\n");
    return RemezInvalidParameter;
  }

  if (symmetry == POSITIVE) {
    if (N % 2) {
      for (uint32_t n = 0; n < N; n++) {
        double val = A[0];
        double x = PiTimes2 * (n - M) / N;
        for (uint32_t k = 1; k <= M; k++) {
          val += 2.0 * A[k] * cos(x * k);
        }
        outTaps[n] = val / N;
      }
    } else {
      for (uint32_t n = 0; n < N; n++) {
        double val = A[0];
        double x = PiTimes2 * (n - M) / N;
        for (uint32_t k = 1; k <= (N / 2 - 1); k++) {
          val += 2.0 * A[k] * cos(x * k);
        }
        outTaps[n] = val / N;
      }
    }
  } else {
    if (N % 2) {
      for (uint32_t n = 0; n < N; n++) {
        double val = 0;
        double x = PiTimes2 * (n - M) / N;
        for (uint32_t k = 1; k <= M; k++) {
          val += 2.0 * A[k] * sin(x * k);
        }
        outTaps[n] = val / N;
      }
    } else {
      for (uint32_t n = 0; n < N; n++) {
        double val = A[N / 2] * sin(Pi * (n - M));
        double x = PiTimes2 * (n - M) / N;
        for (uint32_t k = 1; k <= (N / 2 - 1); k++) {
          val += 2.0 * A[k] * sin(x * k);
        }
        outTaps[n] = val / N;
      }
    }
  }

  return RemezSuccess;
}

/*******************
 * isDone
 *========
 * Checks to see if the error function is small enough to consider
 * the result to have converged.
 *
 * INPUT:
 * ------
 * size_t    r     - 1/2 the number of filter coeffiecients
 * size_t    Ext[] - Indexes to extremal frequencies [r+1]
 * double E[]   - Error function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * Returns 1 if the result converged
 * Returns 0 if the result has not converged
 ********************/

static int isDone(
    size_t r,
    const size_t* Ext,
    size_t ExtSize,
    const double* E,
    size_t ESize) {
  double min, max, current;

  if (r >= ExtSize) {
    fprintf(stderr, "r [%zu] is out of range. Max [%zu].\n", r, ExtSize - 1);
    return RemezInvalidParameter;
  }

  min = max = fabs(E[Ext[0]]);
  for (size_t i = 1; i <= r; i++) {
    current = fabs(E[Ext[i]]);

    if (current < min) {
      min = current;
    }

    if (current > max) {
      max = current;
    }
  }

  return (((max - min) / max) < 0.0001);
}
// remez(double h[], int *numtaps,
//      int *numband, const double bands[],
//      const double des[], const double weight[],
//      int *type, int *griddensity)
// assert(remez(15,[0,0.3,0.4,1],[1,1,0,0]),b,1e-14);
RemezStatus remez(
    double* outTaps,
    size_t outTapsLength,
    const RemezBand* bands,
    size_t bandCount,
    RemezFilterType type,
    size_t gridDensity,
    size_t maxIterations) {
  RemezStatus status = RemezSuccess;

  double* Grid = NULL;
  double* D = NULL;
  double* W = NULL;
  double* E = NULL;
  size_t* Ext = NULL;
  double* taps = NULL;
  double* x = NULL;
  double* y = NULL;
  double* ad = NULL;
  double c;

  Symmetry symmetry;

  size_t r, gridSize, eSize, dSize, wSize, extSize, tapsSize, xSize, ySize,
      adSize;

  if (outTapsLength == 0 || bandCount == 0 || gridDensity <= 0
      || maxIterations <= 0) {
    status = RemezInvalidParameter;
    goto end;
  }

  if (type == RemezFilterTypeBandpass) {
    symmetry = POSITIVE;
  } else {
    symmetry = NEGATIVE;
  }

  r = outTapsLength / 2; /* number of extrema */
  if ((outTapsLength % 2) && (symmetry == POSITIVE)) {
    r++;
  }

  outTaps[0] = 32;
  /*
   * Predict dense grid size in advance for memory allocation
   *   .5 is so we round up, not truncate
   */
  gridSize = 0;
  for (size_t i = 0; i < bandCount; i++) {
    const RemezBand* band = &bands[i];
    double bandwidth = band->highFrequency - band->lowFrequency;
    double addToGridDensity = 2.0 * r * gridDensity * bandwidth;
    size_t addToGridDensityLong = lrint(addToGridDensity);

    gridSize += addToGridDensityLong;
  }

  if (symmetry == NEGATIVE) {
    gridSize--;
  }

  dSize = gridSize;
  wSize = gridSize;
  eSize = gridSize;
  extSize = r + 1;
  tapsSize = r + 1;
  xSize = r + 1;
  ySize = r + 1;
  adSize = r + 1;

  /*
   * Dynamically allocate memory for arrays with proper sizes
   */
  Grid = calloc(gridSize, sizeof(double));
  D = calloc(dSize, sizeof(double));
  W = calloc(wSize, sizeof(double));
  E = calloc(eSize, sizeof(double));
  Ext = calloc(extSize, sizeof(size_t));
  taps = calloc(tapsSize, sizeof(double));
  x = calloc(xSize, sizeof(double));
  y = calloc(ySize, sizeof(double));
  ad = calloc(adSize, sizeof(double));

  if (Grid == NULL || D == NULL || W == NULL || E == NULL || Ext == NULL
      || taps == NULL || x == NULL || y == NULL || ad == NULL) {
    status = RemezOutOfMemory;
    goto end;
  }

  /*
   * Create dense frequency grid
   */
  status = CreateDenseGrid(
      r,
      outTapsLength,
      bands,
      bandCount,
      Grid,
      gridSize,
      D,
      dSize,
      W,
      wSize,
      symmetry,
      gridDensity);
  if (status != RemezSuccess) {
    goto end;
  }

  status = InitialGuess(r, Ext, extSize, gridSize);
  if (status != RemezSuccess) {
    goto end;
  }

  /*
   * For Differentiator: (fix grid)
   */
  if (type == RemezFilterTypeDifferentiator) {
    for (size_t i = 0; i < gridSize; i++) {
      /* D[i] = D[i]*Grid[i]; */
      if (D[i] > 0.0001) {
        W[i] = W[i] / Grid[i];
      }
    }
  }

  /*
   * For odd or Negative symmetry filters, alter the
   * D[] and W[] according to Parks McClellan
   */
  if (symmetry == POSITIVE) {
    if (outTapsLength % 2 == 0) {
      for (size_t i = 0; i < gridSize; i++) {
        c = cos(Pi * Grid[i]);
        D[i] /= c;
        W[i] *= c;
      }
    }
  } else {
    if (outTapsLength % 2) {
      for (size_t i = 0; i < gridSize; i++) {
        c = sin(PiTimes2 * Grid[i]);
        D[i] /= c;
        W[i] *= c;
      }
    } else {
      for (size_t i = 0; i < gridSize; i++) {
        c = sin(Pi * Grid[i]);
        D[i] /= c;
        W[i] *= c;
      }
    }
  }

  /*
   * Perform the Remez Exchange algorithm
   */
  for (size_t iter = 0; iter < maxIterations; iter++) {
    status = CalcParms(
        r,
        Ext,
        extSize,
        Grid,
        gridSize,
        D,
        dSize,
        W,
        wSize,
        ad,
        adSize,
        x,
        xSize,
        y,
        ySize);
    if (status != RemezSuccess) {
      goto end;
    }

    status = CalcError(
        r,
        ad,
        adSize,
        x,
        xSize,
        y,
        ySize,
        Grid,
        gridSize,
        D,
        dSize,
        W,
        wSize,
        E,
        eSize);
    if (status != RemezSuccess) {
      goto end;
    }

    status = Search(r, Ext, extSize, E, eSize);
    if (status != RemezSuccess) {
      goto end;
    }

    for (size_t i = 0; i <= r; i++) {
      if (Ext[i] >= gridSize) {
        fprintf(stderr, "Index [%zu] out of range [0, %zu]\n", Ext[i], r);
        return RemezInternalError;
      }
    }

    if (isDone(r, Ext, extSize, E, eSize)) {
      break;
    }
  }

  if (!isDone(r, Ext, extSize, E, eSize)) {
    status = RemezDidNotConverge;
    goto end;
  }

  status = CalcParms(
      r,
      Ext,
      extSize,
      Grid,
      gridSize,
      D,
      dSize,
      W,
      wSize,
      ad,
      adSize,
      x,
      xSize,
      y,
      ySize);
  if (status != RemezSuccess) {
    goto end;
  }

  /*
   * Find the 'taps' of the filter for use with Frequency
   * Sampling.  If odd or Negative symmetry, fix the taps
   * according to Parks McClellan
   */
  for (size_t i = 0; i <= outTapsLength / 2; i++) {
    double A = 0.0;

    if (symmetry == POSITIVE) {
      if (outTapsLength % 2) {
        c = 1;
      } else {
        c = cos(Pi * (double)i / (double)outTapsLength);
      }
    } else {
      if (outTapsLength % 2) {
        c = sin(PiTimes2 * (double)i / (double)outTapsLength);
      } else {
        c = sin(Pi * (double)i / (double)outTapsLength);
      }
    }

    status = ComputeA(
        (double)i / (double)outTapsLength,
        r,
        ad,
        adSize,
        x,
        xSize,
        y,
        ySize,
        &A);
    if (status != RemezSuccess) {
      goto end;
    }

    taps[i] = A * c;
  }

  /*
   * Frequency sampling design with calculated taps
   */
  status = FreqSample(taps, tapsSize, outTaps, outTapsLength, symmetry);
  if (status != RemezSuccess) {
    goto end;
  }

end:
  return status;
}

static RemezStatus checkFrequency(
    float sampleFrequency,
    float frequency,
    float transitionWidth,
    const char* frequencyName) {
  float nyquistFrequency = sampleFrequency / 2.0f;

  if (frequency > nyquistFrequency) {
    fprintf(
        stderr,
        "%s frequency [%f] cannot exceed the Nyquist frequency [%f] for "
        "sample frequency [%f]\n",
        frequencyName,
        frequency,
        nyquistFrequency,
        sampleFrequency);

    return RemezInvalidParameter;
  } else if (frequency + transitionWidth > nyquistFrequency) {
    fprintf(
        stderr,
        "%s frequency [%f] plus the transition width [%f] cannot exceed the "
        "Nyquist frequency [%f] for "
        "sample frequency [%f]\n",
        frequencyName,
        transitionWidth,
        frequency,
        nyquistFrequency,
        sampleFrequency);

    return RemezInvalidParameter;
  }

  return RemezSuccess;
}

RemezStatus remezGenerateLowPassTaps(
    double sampleFrequency,
    double cutoffFrequency,
    double transitionBandwidth,
    double dbAttenuation,
    double* taps,
    size_t numTaps) {
  RemezStatus status = RemezSuccess;
  const size_t numBands = 2;
  RemezBand bands[numBands];

  // kaiserWindowLength(dbAttenuation, transitionBandwidth)
  //  remez() takes normalized frequencies where 0.5 is the Nyquist rate.
  const double nyquistFrequency = sampleFrequency / 2.0;

  if (cutoffFrequency > nyquistFrequency) {
    fprintf(
        stderr,
        "Cutoff frequency [%f] cannot exceed the Nyquist frequency [%f] for "
        "sample frequency [%f]\n",
        cutoffFrequency,
        nyquistFrequency,
        sampleFrequency);

    return RemezInvalidParameter;
  } else if (cutoffFrequency + transitionBandwidth > nyquistFrequency) {
    fprintf(
        stderr,
        "Cutoff frequency [%f] + the transition bandwidth [%f] = [%f] cannot "
        "exceed the Nyquist frequency [%f] for "
        "sample frequency [%f]\n",
        cutoffFrequency,
        transitionBandwidth,
        cutoffFrequency + transitionBandwidth,
        nyquistFrequency,
        sampleFrequency);

    return RemezInvalidParameter;
  }

  const double relativeCutoffFrequency =
      remezNormalizeFrequency(cutoffFrequency, sampleFrequency);
  const double relativeTransitionWidth =
      remezNormalizeFrequency(transitionBandwidth, sampleFrequency);

  memset(bands, 0, sizeof(bands));

  bands[0].lowFrequency = 0.0f;
  bands[0].highFrequency = relativeCutoffFrequency;
  bands[0].lowFrequencyResponse = 1.0f;
  bands[0].highFrequencyResponse = 1.0f;
  bands[0].weight = 1.0f;

  bands[1].lowFrequency = relativeCutoffFrequency + relativeTransitionWidth;
  bands[1].highFrequency = 0.5f;
  bands[1].lowFrequencyResponse = 0.0f;
  bands[1].highFrequencyResponse = 0.0f;
  bands[1].weight = 1.0f;

  const int32_t gridDensity = 16;
  const int32_t maxIterations = 1000;
  status = remez(
      taps,
      numTaps,
      bands,
      numBands,
      RemezFilterTypeBandpass,
      gridDensity,
      maxIterations);
  if (status != RemezSuccess) {
    return status;
  }

  return RemezSuccess;
}

RemezStatus remezGenerateSingleBandPassTaps(
    double sampleFrequency,
    double lowCutoffFrequency,
    double highCutoffFrequency,
    double transitionWidth,
    double dbAttenuation,
    double* taps,
    size_t numTaps) {
  RemezStatus status = RemezSuccess;
  const size_t numBands = 3;
  RemezBand bands[numBands];

  // kaiserWindowLength(dbAttenuation, transitionWidth)
  //  remez() takes normalized frequencies where 0.5 is the Nyquist rate.
  const double nyquistFrequency = sampleFrequency / 2.0;

  RET_ON_ERR(checkFrequency(
      sampleFrequency,
      lowCutoffFrequency,
      transitionWidth,
      "Low cutoff frequency"));

  RET_ON_ERR(checkFrequency(
      sampleFrequency,
      highCutoffFrequency,
      transitionWidth,
      "High cutoff frequency"));

  memset(bands, 0, sizeof(bands));

  bands[0].lowFrequency = 0.0f;
  bands[0].highFrequency = remezNormalizeFrequency(
      lowCutoffFrequency - transitionWidth,
      sampleFrequency);
  bands[0].lowFrequencyResponse = 0.0f;
  bands[0].highFrequencyResponse = 0.0f;
  bands[0].weight = 1.0f;

  bands[1].lowFrequency =
      remezNormalizeFrequency(lowCutoffFrequency, sampleFrequency);
  bands[1].highFrequency =
      remezNormalizeFrequency(highCutoffFrequency, sampleFrequency);
  bands[1].lowFrequencyResponse = 1.0f;
  bands[1].highFrequencyResponse = 1.0f;
  bands[1].weight = 1.0f;

  bands[2].lowFrequency = remezNormalizeFrequency(
      highCutoffFrequency + transitionWidth,
      sampleFrequency);
  bands[2].highFrequency = 0.5f;
  bands[2].lowFrequencyResponse = 0.0f;
  bands[2].highFrequencyResponse = 0.0f;
  bands[2].weight = 1.0f;

  const int32_t gridDensity = 16;
  const int32_t maxIterations = 1000;
  status = remez(
      taps,
      numTaps,
      bands,
      numBands,
      RemezFilterTypeBandpass,
      gridDensity,
      maxIterations);
  if (status != RemezSuccess) {
    return status;
  }

  return RemezSuccess;
}

double remezNormalizeFrequency(
    double frequencyToNormalize,
    double sampleFrequency) {
  return frequencyToNormalize / sampleFrequency;
}
