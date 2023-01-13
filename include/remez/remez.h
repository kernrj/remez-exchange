/*
 * Parks-McClellan algorithm for FIR filter design (C version)
 * -----------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef REMEZ_REMEZ_H_
#define REMEZ_REMEZ_H_

#include <stddef.h>
#include "remez_export.h"

#ifdef __cplusplus
#define REMEZ_C_LINKAGE extern "C"
#else
#define REMEZ_C_LINKAGE
#endif

// type parameter for remez()
typedef enum {
  RemezFilterTypeBandpass,
  RemezFilterTypeDifferentiator
} RemezFilterType;

typedef enum {
  RemezSuccess,
  RemezInvalidParameter,
  RemezDidNotConverge,
  RemezTooManyExtremalFrequencies,
  RemezTooFewExtremalFrequencies,
  RemezInternalError,
  RemezOutOfMemory,
} RemezStatus;

typedef struct RemezBand {
  double lowFrequency;
  double highFrequency;
  double lowFrequencyResponse;
  double highFrequencyResponse;
  double weight;
} RemezBand;

/**
 * Get a human-readable string for a status code.
 */
REMEZ_C_LINKAGE
REMEZ_PUBLIC
const char* remezStatusToString(RemezStatus status);

/**
 * Calculates the optimal (in the Chebyshev/minimax sense)
 * FIR filter impulse response given a set of band edges,
 * the desired reponse on those bands, and the weight given to
 * the error in those bands.
 *
 * INPUT:
 * ------
 * @param outTaps       Array to be populated with FIR taps
 * @param numOutTaps    Number of elements in the outTaps array
 * @param bands         Specifies the frequency bands and their responses
 * @param numBands      Number of elements in the bands array.
 *                      Frequencies are normalized by the sampling rate, and
 *                      must be in the range: 0 <= frequency <= 0.5.
 *                      remezNormalizeFrequency() can be used to compute correct
 *                      values.
 * @param type          Type of filter
 * @param griddensity   Determines how accurately the filter will be
 *                      constructed.
 *                      The minimum value is 16, but higher numbers are slower
 *                      to compute.
 * @param maxIterations The number of attempts to refine error before
 *
 * @return RemezSuccess when the calculation was successful. Other values
 *         indicate an error.
 */
REMEZ_C_LINKAGE
REMEZ_PUBLIC
RemezStatus remez(
    double* outTaps,
    size_t numOutTaps,
    const RemezBand* bands,
    size_t numBands,
    RemezFilterType type,
    size_t gridDensity,
    size_t maxIterations);

/**
 * Convenience method to generate FIR taps for a low-pass filter.
 *
 * Frequencies in the pass-band are passed-through to the output. For a low-pass
 * filter, the pass-band frequencies range from 0 to cutoffFrequency.
 *
 * Frequencies in the stop-band are filtered out. For a low-pass filter, the
 * stop-band frequencies range from (cutoffFrequency + transitionBandwidth) to
 * (samplingFrequency / 2).
 *
 * @param sampleFrequency     The RF sampling frequency.
 * @param cutoffFrequency     The top frequency of the pass-band. This is not
 *                            normalized, so for 100 kHz, pass 100000.
 * @param transitionBandwidth The amount of transition bandwidth between the
 *                            pass-band and stop-band.
 * @param dbAttenuation       The amount of attenuation in the stop-band
 *                            (e.g. -60).
 * @param taps                The calculated taps are written to this buffer.
 * @param tapsLength          The length of the taps buffer.
 *
 * @return RemezSuccess when the calculation was successful. Other values
 *         indicate an error.
 */
REMEZ_C_LINKAGE
REMEZ_PUBLIC
RemezStatus remezGenerateLowPassTaps(
    double sampleFrequency,
    double cutoffFrequency,
    double transitionBandwidth,
    double dbAttenuation,
    double* taps,
    size_t tapsLength);

REMEZ_C_LINKAGE
REMEZ_PUBLIC
double remezNormalizeFrequency(
    double frequencyToNormalize,
    double samplingFrequency);

#endif  // REMEZ_REMEZ_H_
