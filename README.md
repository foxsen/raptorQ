# RaptorQ (RFC 6330) — C++ Implementation

A C++ implementation of the [RaptorQ Forward Error Correction (FEC)](https://datatracker.ietf.org/doc/html/rfc6330) fountain code algorithm, originally proposed by Qualcomm.

RaptorQ enables reliable data transmission over lossy channels (wireless links, satellite, UDP-based protocols, etc.) by generating repair symbols from source data. A receiver can reconstruct the original data from *any* sufficient subset of received symbols — it does not matter *which* packets are lost, only how many.

This implementation was originally developed for a wireless communication project and later open-sourced.

## Features

- Full RFC 6330 encoding and decoding with GF(256) arithmetic
- Configurable source block size (K) and symbol size (T)
- Adjustable overhead for tunable reliability vs. bandwidth trade-off
- Both Gaussian elimination and RFC 6330 recommended decoding algorithms
- Optimized GF(256) multiplication via precomputed 256x256 lookup table
- Comprehensive test suite (44 tests) covering correctness and throughput benchmarks
- Static library target for integration into other projects

## How It Works

```
Source Data ──► Encoder ──► Source Symbols + Repair Symbols
                                    │
                            (lossy channel)
                                    │
                Received Symbols ──► Decoder ──► Recovered Data
```

**Encoding:** The source data is split into K symbols of T bytes each. The encoder constructs a matrix with LDPC, Half-rate LDPC, and LT (Luby Transform) components, then generates intermediate symbols via a decoding step. From these, any number of repair symbols can be produced.

**Decoding:** Given at least K received symbols (any mix of source and repair), the decoder reconstructs the encoding matrix, solves for the intermediate symbols, and recovers any missing source symbols.

A key insight of RFC 6330 is that encoding and decoding share the same underlying computation — encoding itself requires an internal decoding step to produce intermediate symbols.

## Building

Requires g++ with C++ standard library support. No external dependencies.

```bash
make            # Build the demo executable (./main)
make libraptorq # Build static library (libraptorq.a)
make tests      # Build and run the full test suite
make clean      # Remove all build artifacts
```

## Running

```bash
./main          # Encode K=8 symbols, simulate 10% packet loss, decode, verify
./test_runner   # Run all 44 tests (correctness + bandwidth benchmarks)
```

## Performance

Measured throughput (K=100, T=1024 bytes, compiled with `-O2`):

| Operation | Throughput |
|-----------|-----------|
| Encoding  | ~50 MB/s  |
| Decoding  | ~25 MB/s  |

Throughput scales with symbol size T — larger symbols amortize the fixed cost of matrix construction.

## Project Structure

| File | Description |
|------|-------------|
| `Symbol.h/.cpp` | Symbol abstraction with GF(256) arithmetic (XOR, multiply, divide, multiply-accumulate) |
| `Tables.h` | Lookup tables — pseudorandom generators (`V0`–`V3`, `f[]`), GF(256) log/exp tables, and precomputed 256x256 multiply table |
| `Helper.h/.cpp` | RFC 6330 parameter calculation and partitioning. Optional — K and T can be specified directly |
| `Generators.h/.cpp` | Core algorithm engine (~1100 lines). Builds the encoding matrix A in stages (GLDPC, GHDPC, GLT), solves via Gaussian elimination or RFC 6330 algorithm. Handles intermediate symbol generation, repair symbol generation, and source symbol recovery |
| `Encoder.h/.cpp` | Encoding wrapper — calls `generate_intermediates()` then `generate_repairs()` |
| `Decoder.h/.cpp` | Decoding wrapper — calls `generate_intermediates()` then `recover()` for lost symbols |
| `Main.cpp` | Demo program — encodes, simulates packet loss, decodes, verifies correctness |
| `test.cpp` | Test suite — GF(256) arithmetic, symbol operations, round-trip encode/decode, and bandwidth benchmarks |

## API Usage

```cpp
#include "Helper.h"
#include "Symbol.h"

int K = 50;        // number of source symbols
int T = 1024;      // symbol size in bytes (must be multiple of 4)
int overhead = 10; // number of repair symbols to generate

// Prepare source data
char **source = new char*[K];
for (int i = 0; i < K; i++) {
    source[i] = new char[T];
    // ... fill with data ...
}

// Encode
Encoder *enc = new Encoder();
enc->init(K, T);
Symbol **repairs = enc->encode(source, overhead);
// repairs[0..overhead-1] are the repair symbols to send alongside source

// Decode (given received symbols and their ESIs)
Decoder *dec = new Decoder();
dec->init(K, T);
dec->decode(received, n_received, esi_list);
Symbol *recovered = dec->recover(lost_index);  // recover a specific lost symbol
```

## License

See [LICENSE](LICENSE) for details.
