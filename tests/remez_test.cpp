/*
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

#include <gtest/gtest.h>
#include <remez/remez.h>

using namespace std;

static string doubleVectorToString(const vector<double>& v) {
  ostringstream out;
  out << "{";

  for (size_t i = 0; i < v.size(); i++) {
    if (i > 0) {
      out << ", ";
    }

    out << v[i];
  }

  out << "}";

  return out.str();
}

TEST(WhenLowPassTapsAreCreated, TheyAreCorrect) {
  const vector<double> expectedTaps = {
      0.0415131831103279,
      0.0581639884202646,
      -0.0281579212691008,
      -0.0535575358002337,
      -0.0617245915143180,
      0.0507753178978075,
      0.2079018331396460,
      0.3327160895375440,
      0.3327160895375440,
      0.2079018331396460,
      0.0507753178978075,
      -0.0617245915143180,
      -0.0535575358002337,
      -0.0281579212691008,
      0.0581639884202646,
      0.0415131831103279,
  };

  vector<double> taps(expectedTaps.size());

  RemezBand bands[] = {
      RemezBand {
          .lowFrequency = 0.0,
          .highFrequency = 0.15,
          .lowFrequencyResponse = 1.0,
          .highFrequencyResponse = 1.0,
          .weight = 1.0,
      },
      RemezBand {
          .lowFrequency = 0.2,
          .highFrequency = 0.5,
          .lowFrequencyResponse = 0.0,
          .highFrequencyResponse = 0.0,
          .weight = 1.0,
      },
  };

  const size_t bandCount = sizeof(bands) / sizeof(bands[0]);

  RemezStatus status = remez(
      taps.data(),
      taps.size(),
      bands,
      bandCount,
      RemezFilterTypeBandpass,
      16,
      40);

  ASSERT_EQ(RemezSuccess, status);

  bool passed = true;
  for (size_t i = 0; i < taps.size(); i++) {
    passed &= abs(taps[i] - expectedTaps[i]) < 1e-14;
  }

  if (!passed) {
    EXPECT_TRUE(passed) << "Expected taps" << endl
                        << doubleVectorToString(expectedTaps) << endl
                        << "but got" << endl
                        << doubleVectorToString(taps);
  }
}
