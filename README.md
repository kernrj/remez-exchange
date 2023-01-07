# remez-exchange
A library for low-pass filter and general FIR filter design using the Parks-McClellan algorithm.

Taps for a low-pass FIR filter can be generated using `remezGenerateLowPassTaps()`.

Troubleshooting:

On a Mac, the error `No CMAKE_C_COMPILER could be found` when loading the CMake file can be fixed by adding the following the cmake command `-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++`
