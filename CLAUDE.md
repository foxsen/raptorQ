# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C++ implementation of the RaptorQ Forward Error Correction (FEC) algorithm per RFC 6330 (not RFC 5053 Raptor). Provides encoding and decoding for reliable data transmission over lossy channels using Galois field GF(256) arithmetic.

## Build Commands

```bash
make            # Build the test executable (./main)
make libraptorq # Build static library (libraptorq.a)
make clean      # Remove all build artifacts
```

Compiler: g++ with `-fPIC -O2`. No external dependencies beyond standard C++ libraries.

## Running

```bash
./main          # Run the test/demo program (encodes, simulates packet loss, decodes, verifies)
```

## Architecture

**Data flow:** Source symbols → Encoder → (channel with loss) → Decoder → Recovered symbols

### Core layers (bottom-up):

- **Symbol** (`Symbol.h/.cpp`): Symbol abstraction holding payload data. Implements GF(256) arithmetic operations (XOR, multiply, divide, multiply-accumulate) used throughout the codec.
- **Tables** (`Tables.h`): Lookup tables for the RaptorQ algorithm — pseudorandom value tables (`V0[]`, `f[]`) and GF(256) log/exp/multiply tables (`OCT_*`).
- **Helper** (`Helper.h/.cpp`): RFC 6330 parameter calculation (partitioning, symbol counts). Currently bypassed in practice — K and T are specified directly.
- **Generators** (`Generators.h/.cpp`): Core algorithm engine (~1050 lines). Builds the encoding matrix A in stages (`_0_init` → `_1_Tuples` → `_2_Matrix_GLDPC` → `_3_Matrix_GHDPC` → `_4_Matrix_GLT`), then solves via Gaussian elimination. Handles intermediate symbol generation, repair symbol generation, and source symbol recovery.
- **Encoder/Decoder** (`Encoder.h/.cpp`, `Decoder.h/.cpp`): Thin wrappers around `Generators`. Both use the same matrix construction; encoding calls `generate_intermediates()` + `generate_repairs()`, decoding calls `generate_intermediates()` + `recover()`.
- **Main** (`Main.cpp`): Test harness — creates K=8 symbols of T=4 bytes, encodes with overhead, simulates 10% packet loss, decodes, and verifies correctness with timing output.

### Key design points:

- Encoder and Decoder share the same underlying `Generators` class and matrix construction — the difference is which symbols are known vs. to-be-recovered.
- Matrix A has three components: G_LDPC (low-density parity check), G_HDPC (half-rate LDPC), and G_LT (Luby Transform).
- Conditional `#ifdef SPARSE` support exists for sparse matrix representation but is not enabled by default.
- `MAX_OVERHEAD` constant (40) controls pre-calculated repair tuple limit.
- `status` variable in Generators tracks codec state progression (0→3).
