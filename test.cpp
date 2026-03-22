#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "Helper.h"
#include "Symbol.h"
#include "Generators.h"

using namespace std;

/* ── Minimal test framework ──────────────────────────────────── */

static int tests_run    = 0;
static int tests_passed = 0;
static int tests_failed = 0;
static const char *current_test = "";

#define TEST(name)                                      \
    do {                                                \
        current_test = #name;                           \
        tests_run++;                                    \
        cout << "  " << #name << " ... ";               \
        name();                                         \
        tests_passed++;                                 \
        cout << "OK" << endl;                           \
    } while(0)

/* ASSERT prints file:line, records failure, and returns from the
   calling test function so later asserts don't crash on bad state. */
#define ASSERT(cond)                                                    \
    do {                                                                \
        if (!(cond)) {                                                  \
            cout << "FAIL" << endl;                                     \
            cerr << "    " << __FILE__ << ":" << __LINE__               \
                 << ": ASSERT(" << #cond << ") failed" << endl;         \
            tests_failed++;                                             \
            tests_passed--;                                             \
            return;                                                     \
        }                                                               \
    } while(0)

#define ASSERT_EQ(a, b)                                                 \
    do {                                                                \
        if ((a) != (b)) {                                               \
            cout << "FAIL" << endl;                                     \
            cerr << "    " << __FILE__ << ":" << __LINE__               \
                 << ": ASSERT_EQ(" << #a << ", " << #b << ") failed: "  \
                 << (a) << " != " << (b) << endl;                       \
            tests_failed++;                                             \
            tests_passed--;                                             \
            return;                                                     \
        }                                                               \
    } while(0)

/* ── GF(256) arithmetic tests ──────────────────────────────────── */

void test_octmul_identity() {
    for (int a = 0; a < 256; a++) {
        ASSERT_EQ(octmul(a, 1), a);
        ASSERT_EQ(octmul(1, a), a);
    }
}

void test_octmul_zero() {
    for (int a = 0; a < 256; a++) {
        ASSERT_EQ(octmul(a, 0), 0);
        ASSERT_EQ(octmul(0, a), 0);
    }
}

void test_octmul_commutative() {
    for (int a = 1; a < 256; a += 17)
        for (int b = 1; b < 256; b += 13)
            ASSERT_EQ(octmul(a, b), octmul(b, a));
}

void test_octdiv_identity() {
    for (int a = 0; a < 256; a++)
        ASSERT_EQ(octdiv(a, 1), a);
}

void test_octdiv_zero_numerator() {
    for (int b = 1; b < 256; b++)
        ASSERT_EQ(octdiv(0, b), 0);
}

void test_octmul_octdiv_inverse() {
    /* (a * b) / b == a  for all non-zero a,b */
    for (int a = 1; a < 256; a += 7)
        for (int b = 1; b < 256; b += 11) {
            unsigned char prod = octmul(a, b);
            ASSERT_EQ(octdiv(prod, b), a);
        }
}

void test_octdiv_self() {
    /* a / a == 1  for all non-zero a */
    for (int a = 1; a < 256; a++)
        ASSERT_EQ(octdiv(a, a), 1);
}

/* ── Deg() tests ───────────────────────────────────────────── */

/*
 * RFC 6330 Section 5.3.5.2: d = min(j, W-2) where j is the minimum
 * value such that v < f[j].  f[31] = {0, 5243, 529531, ..., 1048576}.
 * Deg() needs a Generators object with W initialized.
 */

void test_deg_at_boundary() {
    /* Use K=8 which gives W=17, so min(j, W-2) = min(j, 15) */
    Generators g;
    g.gen(8, 8, 4);

    /* f[0]=0: v=0 must give j=1 (not 0), since we need v < f[j] */
    ASSERT_EQ((int)g.Deg(0), 1);

    /* f[1]=5243: v=5242 gives j=1, v=5243 must give j=2 (not 1) */
    ASSERT_EQ((int)g.Deg(5242), 1);
    ASSERT_EQ((int)g.Deg(5243), 2);

    /* f[2]=529531: v=529530 gives j=2, v=529531 must give j=3 */
    ASSERT_EQ((int)g.Deg(529530), 2);
    ASSERT_EQ((int)g.Deg(529531), 3);
}

void test_deg_interior_values() {
    Generators g;
    g.gen(8, 8, 4);  /* W=17, so cap is min(j, 15) */

    /* v=1 is in [f[0]=0, f[1]=5243), so j=1 */
    ASSERT_EQ((int)g.Deg(1), 1);

    /* v=100000 is in [f[1]=5243, f[2]=529531), so j=2 */
    ASSERT_EQ((int)g.Deg(100000), 2);

    /* v=1048575 is in [f[29]=1017662, f[30]=1048576), so j=30 -> min(30,15)=15 */
    ASSERT_EQ((int)g.Deg(1048575), 15);
}

/* ── Symbol tests ──────────────────────────────────────────── */

void test_symbol_default_constructor() {
    Symbol s;
    ASSERT_EQ(s.nbytes, 0);
    ASSERT(s.data == NULL);
    ASSERT_EQ(s.sbn, -1);
    ASSERT_EQ(s.esi, -1);
}

void test_symbol_sized_constructor() {
    Symbol s(8);
    ASSERT_EQ(s.nbytes, 8);
    ASSERT(s.data != NULL);
    /* should be zeroed */
    for (int i = 0; i < 8; i++)
        ASSERT_EQ(((unsigned char*)s.data)[i], 0);
}

void test_symbol_init_and_fill() {
    Symbol s;
    s.init(4);
    ASSERT_EQ(s.nbytes, 4);

    char buf[4] = {0x11, 0x22, 0x33, 0x44};
    s.fillData(buf, 4);
    ASSERT(memcmp(s.data, buf, 4) == 0);
}

void test_symbol_assign() {
    Symbol a(4), b(4);
    char buf[4] = {(char)0xAA, (char)0xBB, (char)0xCC, (char)0xDD};
    a.fillData(buf, 4);

    b = a;
    ASSERT_EQ(b.nbytes, 4);
    ASSERT(memcmp(b.data, buf, 4) == 0);
    /* verify it's a copy, not shared pointer */
    ASSERT(a.data != b.data);
}

void test_symbol_xor() {
    Symbol a(4), b(4);
    char abuf[4] = {(char)0xFF, 0x00, (char)0xAA, 0x55};
    char bbuf[4] = {0x0F, (char)0xF0, 0x55, (char)0xAA};
    a.fillData(abuf, 4);
    b.fillData(bbuf, 4);

    a.xxor(&b);

    unsigned char *p = (unsigned char *)a.data;
    ASSERT_EQ(p[0], 0xF0);
    ASSERT_EQ(p[1], 0xF0);
    ASSERT_EQ(p[2], 0xFF);
    ASSERT_EQ(p[3], 0xFF);
}

void test_symbol_mul_div() {
    Symbol s(4);
    char buf[4] = {0x00, 0x01, 0x53, (char)0xFF};
    s.fillData(buf, 4);

    /* save original */
    char orig[4];
    memcpy(orig, s.data, 4);

    unsigned char factor = 0x37;
    s.mul(factor);
    s.div(factor);

    /* mul then div by same value should round-trip */
    ASSERT(memcmp(s.data, orig, 4) == 0);
}

void test_symbol_muladd() {
    /* muladd(s, u) computes: this ^= s * u  (in GF(256), byte-wise) */
    Symbol a(4), b(4);
    char abuf[4] = {0, 0, 0, 0};
    char bbuf[4] = {0x01, 0x02, 0x03, 0x04};
    a.fillData(abuf, 4);
    b.fillData(bbuf, 4);

    /* a ^= b * 1  should give a == b */
    a.muladd(&b, 1);
    ASSERT(memcmp(a.data, b.data, 4) == 0);

    /* a ^= b * 1 again should give zero (a ^ b ^ b == 0) */
    a.muladd(&b, 1);
    for (int i = 0; i < 4; i++)
        ASSERT_EQ(((unsigned char*)a.data)[i], 0);
}

/* ── Helper / Combin tests ─────────────────────────────────── */

void test_combin_base_cases() {
    ASSERT_EQ(Combin(5, 0), 1);
    ASSERT_EQ(Combin(5, 1), 5);
    ASSERT_EQ(Combin(5, 5), 1);  /* C(n,n)=1 via Combin(5,0) */
}

void test_combin_known_values() {
    ASSERT_EQ(Combin(10, 2), 45);
    ASSERT_EQ(Combin(10, 3), 120);
    ASSERT_EQ(Combin(6, 3),  20);
    ASSERT_EQ(Combin(7, 4),  35);
}

void test_combin_symmetry() {
    /* C(n,k) == C(n, n-k) */
    ASSERT_EQ(Combin(8, 3), Combin(8, 5));
    ASSERT_EQ(Combin(10, 4), Combin(10, 6));
}

void test_helper_partition() {
    Helper h;
    PartitionS p = h.partition(13, 4);
    /* 13 / 4: ceil=4, floor=3, JL = 13 - 3*4 = 1, JS = 4 - 1 = 3 */
    ASSERT_EQ(p.IL, 4);
    ASSERT_EQ(p.IS, 3);
    ASSERT_EQ(p.JL, 1);
    ASSERT_EQ(p.JS, 3);
}

void test_helper_partition_even() {
    Helper h;
    PartitionS p = h.partition(12, 4);
    ASSERT_EQ(p.IL, 3);
    ASSERT_EQ(p.IS, 3);
    ASSERT_EQ(p.JL, 0);
    ASSERT_EQ(p.JS, 4);
}

/* ── Encoder-only tests ────────────────────────────────────── */

void test_encoder_init_valid() {
    Encoder *enc = new Encoder();
    ASSERT(enc->init(8, 4));
    delete enc;
}

void test_encoder_init_invalid_K() {
    Encoder *enc = new Encoder();
    ASSERT(!enc->init(0, 4));
    delete enc;
}

void test_encoder_produces_repairs() {
    int K = 8, T = 4, overhead = 4;

    char **source = new char*[K];
    for (int i = 0; i < K; i++) {
        source[i] = new char[T];
        for (int j = 0; j < T; j++)
            source[i][j] = (char)(i * 10 + j);
    }

    Encoder *enc = new Encoder();
    enc->init(K, T);
    Symbol **repairs = enc->encode(source, overhead);
    ASSERT(repairs != NULL);

    /* each repair symbol should have the right size */
    for (int i = 0; i < overhead; i++) {
        ASSERT(repairs[i] != NULL);
        ASSERT_EQ(repairs[i]->nbytes, T);
        delete repairs[i];
    }
    delete[] repairs;
    delete enc;

    for (int i = 0; i < K; i++)
        delete[] source[i];
    delete[] source;
}

/* ── Full encode→decode round-trip helpers ─────────────────── */

/*
 * Encode K symbols of size T, drop specific source symbols listed in
 * lost_indices[], supply 'overhead' extra repair symbols, decode, and
 * verify every lost symbol is recovered exactly.
 * Returns true on success.
 */
static bool roundtrip(int K, int T, int overhead,
                      int *lost_indices, int nlost)
{
    int i, j;

    /* build deterministic source data */
    char **source = new char*[K];
    for (i = 0; i < K; i++) {
        source[i] = new char[T];
        for (j = 0; j < T; j++)
            source[i][j] = (char)((i * 7 + j * 3) & 0xFF);
    }

    /* encode */
    Encoder *enc = new Encoder();
    if (!enc->init(K, T)) { delete enc; return false; }
    Symbol **repairs = enc->encode(source, overhead);
    if (!repairs) { delete enc; return false; }

    /* build received set: source symbols minus lost ones, then all repairs */
    int max_received = K + overhead;
    char **received = new char*[max_received];
    int *esi = new int[max_received];
    int n = 0;

    for (i = 0; i < K; i++) {
        bool is_lost = false;
        for (j = 0; j < nlost; j++)
            if (lost_indices[j] == i) { is_lost = true; break; }
        if (!is_lost) {
            received[n] = new char[T];
            memcpy(received[n], source[i], T);
            esi[n] = i;
            n++;
        }
    }

    for (i = 0; i < overhead; i++) {
        received[n] = new char[T];
        memcpy(received[n], repairs[i]->data, T);
        esi[n] = K + i;
        n++;
        delete repairs[i];
    }
    delete[] repairs;
    delete enc;

    /* decode */
    bool ok = true;
    if (n < K) {
        ok = false;
    } else {
        Decoder *dec = new Decoder();
        if (!dec->init(K, T)) { ok = false; }
        else if (!dec->decode(received, n, esi)) { ok = false; }
        else {
            for (i = 0; i < nlost; i++) {
                Symbol *s = dec->recover(lost_indices[i]);
                if (!s || memcmp(s->data, source[lost_indices[i]], T) != 0)
                    ok = false;
                delete s;
                if (!ok) break;
            }
        }
        delete dec;
    }

    /* cleanup */
    for (i = 0; i < n; i++)
        delete[] received[i];
    delete[] received;
    delete[] esi;
    for (i = 0; i < K; i++)
        delete[] source[i];
    delete[] source;

    return ok;
}

/* ── Round-trip test cases ─────────────────────────────────── */

void test_roundtrip_no_loss() {
    /* K=8, T=4, 0 lost, 4 overhead — trivial case */
    ASSERT(roundtrip(8, 4, 4, NULL, 0));
}

void test_roundtrip_one_lost() {
    int lost[] = {3};
    ASSERT(roundtrip(8, 4, 4, lost, 1));
}

void test_roundtrip_two_lost() {
    int lost[] = {0, 7};  /* first and last */
    ASSERT(roundtrip(8, 4, 6, lost, 2));
}

void test_roundtrip_three_lost() {
    int lost[] = {1, 4, 6};
    ASSERT(roundtrip(8, 4, 8, lost, 3));
}

void test_roundtrip_larger_symbol() {
    /* T=16 bytes (4 ints) */
    int lost[] = {2, 5};
    ASSERT(roundtrip(8, 16, 6, lost, 2));
}

void test_roundtrip_K4() {
    int lost[] = {1};
    ASSERT(roundtrip(4, 4, 4, lost, 1));
}

void test_roundtrip_K16() {
    int lost[] = {0, 5, 10};
    ASSERT(roundtrip(16, 4, 8, lost, 3));
}

void test_roundtrip_K20() {
    int lost[] = {3, 7, 12, 18};
    ASSERT(roundtrip(20, 8, 10, lost, 4));
}

void test_roundtrip_K50() {
    int lost[] = {0, 12, 25, 37, 49};
    ASSERT(roundtrip(50, 4, 12, lost, 5));
}

void test_roundtrip_K100() {
    int lost[] = {0, 33, 66, 99};
    ASSERT(roundtrip(100, 4, 10, lost, 4));
}

/* Regression: K > 466 caused OCT_EXP out-of-bounds in _3_Matrix_GHDPC.
   OCT_EXP[510] was indexed by (k-j) which can reach K1+S-1 > 509. */
void test_roundtrip_K512() {
    int lost[] = {0, 100, 255, 400, 511};
    ASSERT(roundtrip(512, 4, 12, lost, 5));
}

/* ── Main ──────────────────────────────────────────────────── */

int main()
{
    cout << "=== GF(256) arithmetic ===" << endl;
    TEST(test_octmul_identity);
    TEST(test_octmul_zero);
    TEST(test_octmul_commutative);
    TEST(test_octdiv_identity);
    TEST(test_octdiv_zero_numerator);
    TEST(test_octmul_octdiv_inverse);
    TEST(test_octdiv_self);

    cout << endl << "=== Deg (degree distribution) ===" << endl;
    TEST(test_deg_at_boundary);
    TEST(test_deg_interior_values);

    cout << endl << "=== Symbol ===" << endl;
    TEST(test_symbol_default_constructor);
    TEST(test_symbol_sized_constructor);
    TEST(test_symbol_init_and_fill);
    TEST(test_symbol_assign);
    TEST(test_symbol_xor);
    TEST(test_symbol_mul_div);
    TEST(test_symbol_muladd);

    cout << endl << "=== Helper / Combin ===" << endl;
    TEST(test_combin_base_cases);
    TEST(test_combin_known_values);
    TEST(test_combin_symmetry);
    TEST(test_helper_partition);
    TEST(test_helper_partition_even);

    cout << endl << "=== Encoder ===" << endl;
    TEST(test_encoder_init_valid);
    TEST(test_encoder_init_invalid_K);
    TEST(test_encoder_produces_repairs);

    cout << endl << "=== Encode/Decode round-trip ===" << endl;
    TEST(test_roundtrip_no_loss);
    TEST(test_roundtrip_one_lost);
    TEST(test_roundtrip_two_lost);
    TEST(test_roundtrip_three_lost);
    TEST(test_roundtrip_larger_symbol);
    TEST(test_roundtrip_K4);
    TEST(test_roundtrip_K16);
    TEST(test_roundtrip_K20);
    TEST(test_roundtrip_K50);
    TEST(test_roundtrip_K100);
    TEST(test_roundtrip_K512);

    cout << endl << "================================" << endl;
    cout << tests_run << " tests: "
         << tests_passed << " passed, "
         << tests_failed << " failed" << endl;

    return tests_failed > 0 ? 1 : 0;
}
