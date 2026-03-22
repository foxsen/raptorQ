# RaptorQ (RFC 6330) Implementation

This is an implementation of the RaptorQ Forward Error Correction (FEC) algorithm specified in [RFC 6330](https://datatracker.ietf.org/doc/html/rfc6330).

RaptorQ is an efficient fountain code FEC algorithm originally proposed by Qualcomm. This implementation was originally developed for a wireless communication project. It sat on a hard drive for a long time and is now shared to avoid it going to waste.

## File Organization

- **Symbol.h / Symbol.cpp** — Abstraction of the "symbol" concept from the RFC document, including initialization, XOR, and assignment operations.
- **Helper.h / Helper.cpp** — Computes RFC 6330 parameters and partitioning. Optional — you can bypass it and specify K, T, etc. directly.
- **Generators.h / Generators.cpp** — Core algorithm. Initializes the encoder/decoder based on parameters such as K and T, generates the encoding matrix, intermediate symbols, and repair symbols. Provides both Gaussian elimination and the RFC 6330 recommended decoding algorithm. Also includes a sparse matrix representation intended to reduce memory usage and improve efficiency, but in practice it increased computational overhead and is not actively used. The code contains some comments.
- **Encoder.\* / Decoder.\*** — Wrappers around the core. An interesting aspect of RFC 6330 is that encoding and decoding are essentially the same computation — encoding requires a decoding step first.
- **Main.cpp** — Test program. Invokes the encoder, simulates transmission with packet loss, receives and decodes, and verifies correctness.

## Building

```bash
make            # Build the test executable (./main)
make libraptorq # Build static library (libraptorq.a)
make tests      # Build and run the test suite
make clean      # Remove all build artifacts
```

## Running

```bash
./main          # Run the demo (encode, simulate packet loss, decode, verify)
```
