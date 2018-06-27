
#include <string.h>
#include <stdio.h>
#include <endian.h>

#include <stdint.h>
#if defined(_MSC_VER)
#include <immintrin.h>
#elif defined(__GNUC__)
#include <x86intrin.h>
#endif

__m128i MASK;
__m128i INIT0;
__m128i INIT1;

static const size_t OUTPUT_SIZE = 32;

struct sha256
{
  uint32_t s[8];
  unsigned char buf[64];
  uint64_t bytes;
};

static struct sha256 sha256_g;


inline static uint16_t ReadLE16(const unsigned char* ptr)
{
    uint16_t x;
    memcpy((char*)&x, ptr, 2);
    return le16toh(x);
}

inline static uint32_t ReadLE32(const unsigned char* ptr)
{
    uint32_t x;
    memcpy((char*)&x, ptr, 4);
    return le32toh(x);
}

static inline uint64_t ReadLE64(const unsigned char* ptr)
{
    uint64_t x;
    memcpy((char*)&x, ptr, 8);
    return le64toh(x);
}

static inline void WriteLE16(unsigned char* ptr, uint16_t x)
{
    uint16_t v = htole16(x);
    memcpy(ptr, (char*)&v, 2);
}

static inline void WriteLE32(unsigned char* ptr, uint32_t x)
{
    uint32_t v = htole32(x);
    memcpy(ptr, (char*)&v, 4);
}

static inline void WriteLE64(unsigned char* ptr, uint64_t x)
{
    uint64_t v = htole64(x);
    memcpy(ptr, (char*)&v, 8);
}

static inline uint32_t ReadBE32(const unsigned char* ptr)
{
    uint32_t x;
    memcpy((char*)&x, ptr, 4);
    return be32toh(x);
}

static inline uint32_t ReadBE64(const unsigned char* ptr)
{
    uint64_t x;
    memcpy((char*)&x, ptr, 8);
    return be64toh(x);
}

static inline void WriteBE32(unsigned char* ptr, uint32_t x)
{
    uint32_t v = htobe32(x);
    memcpy(ptr, (char*)&v, 4);
}

static inline void WriteBE64(unsigned char* ptr, uint64_t x)
{
    uint64_t v = htobe64(x);
    memcpy(ptr, (char*)&v, 8);
}

/** Return the smallest number n such that (x >> n) == 0 (or 64 if the highest bit in x is set. */
static inline uint64_t CountBits(uint64_t x)
{
#if HAVE_DECL___BUILTIN_CLZL
    if (sizeof(unsigned long) >= sizeof(uint64_t)) {
        return x ? 8 * sizeof(unsigned long) - __builtin_clzl(x) : 0;
    }
#endif
#if HAVE_DECL___BUILTIN_CLZLL
    if (sizeof(unsigned long long) >= sizeof(uint64_t)) {
        return x ? 8 * sizeof(unsigned long long) - __builtin_clzll(x) : 0;
    }
#endif
    int ret = 0;
    while (x) {
        x >>= 1;
        ++ret;
    }
    return ret;
}

void Transform(uint32_t* s, const unsigned char* chunk, size_t blocks);

void GlobalInit() {
  MASK = _mm_set_epi64x(0x0c0d0e0f08090a0bULL, 0x0405060700010203ULL);
  INIT0 = _mm_set_epi64x(0x6a09e667bb67ae85ull, 0x510e527f9b05688cull);
  INIT1 = _mm_set_epi64x(0x3c6ef372a54ff53aull, 0x1f83d9ab5be0cd19ull);
}

/** Initialize SHA-256 state. */
void Initialize(uint32_t* s)
{
  s[0] = 0x6a09e667ul;
  s[1] = 0xbb67ae85ul;
  s[2] = 0x3c6ef372ul;
  s[3] = 0xa54ff53aul;
  s[4] = 0x510e527ful;
  s[5] = 0x9b05688cul;
  s[6] = 0x1f83d9abul;
  s[7] = 0x5be0cd19ul;
}


void ResetG()
{
  sha256_g.bytes = 0;
  Initialize(sha256_g.s);
}



void Write(struct sha256 *sha256, const unsigned char* data, size_t len)
{
    const unsigned char* end = data + len;
    size_t bufsize = sha256->bytes % 64;
    if (bufsize && bufsize + len >= 64) {
        // Fill the buffer, and process it.
        memcpy(sha256->buf + bufsize, data, 64 - bufsize);
        sha256->bytes += 64 - bufsize;
        data += 64 - bufsize;
        Transform(sha256->s, sha256->buf, 1);
        bufsize = 0;
    }
    if (end - data >= 64) {
        size_t blocks = (end - data) / 64;
        Transform(sha256->s, data, blocks);
        data += 64 * blocks;
        sha256->bytes += 64 * blocks;
    }
    if (end > data) {
        // Fill the buffer with what remains.
        memcpy(sha256->buf + bufsize, data, end - data);
        sha256->bytes += end - data;
    }

}

void WriteG(const unsigned char *data, size_t len) {
  Write(&sha256_g, data, len);
}


void Finalize(struct sha256 *sha256, unsigned char hash[OUTPUT_SIZE])
{
    static const unsigned char pad[64] = {0x80};
    unsigned char sizedesc[8];
    WriteBE64(sizedesc, sha256->bytes << 3);
    Write(sha256, pad, 1 + ((119 - (sha256->bytes % 64)) % 64));
    Write(sha256, sizedesc, 8);
    WriteBE32(hash, sha256->s[0]);
    WriteBE32(hash + 4, sha256->s[1]);
    WriteBE32(hash + 8, sha256->s[2]);
    WriteBE32(hash + 12, sha256->s[3]);
    WriteBE32(hash + 16, sha256->s[4]);
    WriteBE32(hash + 20, sha256->s[5]);
    WriteBE32(hash + 24, sha256->s[6]);
    WriteBE32(hash + 28, sha256->s[7]);
}

void FinalizeG(unsigned char hash[OUTPUT_SIZE]) {
  Finalize(&sha256_g, hash);
}


inline void __attribute__((always_inline)) QuadRound(__m128i* state0, __m128i* state1, uint64_t k1, uint64_t k0)
{
    const __m128i msg = _mm_set_epi64x(k1, k0);
    *state1 = _mm_sha256rnds2_epu32(*state1, *state0, msg);
    *state0 = _mm_sha256rnds2_epu32(*state0, *state1, _mm_shuffle_epi32(msg, 0x0e));
}

inline void __attribute__((always_inline)) QuadRoundm(__m128i* state0, __m128i* state1, __m128i m, uint64_t k1, uint64_t k0)
{
    const __m128i msg = _mm_add_epi32(m, _mm_set_epi64x(k1, k0));
    *state1 = _mm_sha256rnds2_epu32(*state1, *state0, msg);
    *state0 = _mm_sha256rnds2_epu32(*state0, *state1, _mm_shuffle_epi32(msg, 0x0e));
}

inline void __attribute__((always_inline)) ShiftMessageA(__m128i* m0, __m128i m1)
{
    *m0 = _mm_sha256msg1_epu32(*m0, m1);
}

inline void __attribute__((always_inline)) ShiftMessageC(__m128i* m0, __m128i m1, __m128i* m2)
{
    *m2 = _mm_sha256msg2_epu32(_mm_add_epi32(*m2, _mm_alignr_epi8(m1, *m0, 4)), m1);
}

inline void __attribute__((always_inline)) ShiftMessageB(__m128i* m0, __m128i m1, __m128i* m2)
{
    ShiftMessageC(m0, m1, m2);
    ShiftMessageA(m0, m1);
}

inline void __attribute__((always_inline)) Shuffle(__m128i* s0, __m128i* s1)
{
    const __m128i t1 = _mm_shuffle_epi32(*s0, 0xB1);
    const __m128i t2 = _mm_shuffle_epi32(*s1, 0x1B);
    *s0 = _mm_alignr_epi8(t1, t2, 0x08);
    *s1 = _mm_blend_epi16(t2, t1, 0xF0);
}

inline void __attribute__((always_inline)) Unshuffle(__m128i* s0, __m128i* s1)
{
    const __m128i t1 = _mm_shuffle_epi32(*s0, 0x1B);
    const __m128i t2 = _mm_shuffle_epi32(*s1, 0xB1);
    *s0 = _mm_blend_epi16(t1, t2, 0xF0);
    *s1 = _mm_alignr_epi8(t2, t1, 0x08);
}

inline __m128i __attribute__((always_inline)) Load(const unsigned char* in)
{
    return _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)in), MASK);
}

inline void __attribute__((always_inline)) Save(unsigned char* out, __m128i s)
{
    _mm_storeu_si128((__m128i*)out, _mm_shuffle_epi8(s, MASK));
}

typedef void (*TransformType)(uint32_t*, const unsigned char*, size_t);

int SelfTest(TransformType tr) {
    static const unsigned char in1[65] = {0, 0x80};
    static const unsigned char in2[129] = {
        0,
        32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
        32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
        0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0
    };
    static const uint32_t init[8] = {0x6a09e667ul, 0xbb67ae85ul, 0x3c6ef372ul, 0xa54ff53aul, 0x510e527ful, 0x9b05688cul, 0x1f83d9abul, 0x5be0cd19ul};
    static const uint32_t out1[8] = {0xe3b0c442ul, 0x98fc1c14ul, 0x9afbf4c8ul, 0x996fb924ul, 0x27ae41e4ul, 0x649b934cul, 0xa495991bul, 0x7852b855ul};
    static const uint32_t out2[8] = {0xce4153b0ul, 0x147c2a86ul, 0x3ed4298eul, 0xe0676bc8ul, 0x79fc77a1ul, 0x2abe1f49ul, 0xb2b055dful, 0x1069523eul};
    uint32_t buf[8];
    memcpy(buf, init, sizeof(buf));
    // Process nothing, and check we remain in the initial state.
    tr(buf, NULL, 0);
    if (memcmp(buf, init, sizeof(buf))) return 0;
    // Process the padded empty string (unaligned)
    tr(buf, in1 + 1, 1);
    if (memcmp(buf, out1, sizeof(buf))) return 0;
    // Process 64 spaces (unaligned)
    memcpy(buf, init, sizeof(buf));
    tr(buf, in2 + 1, 2);
    if (memcmp(buf, out2, sizeof(buf))) return 0;
    return 1;
}

void Transform(uint32_t* s, const unsigned char* chunk, size_t blocks)
{
    __m128i m0, m1, m2, m3, s0, s1, so0, so1;

    /* Load state */
    s0 = _mm_loadu_si128((const __m128i*)s);
    s1 = _mm_loadu_si128((const __m128i*)(s + 4));
    Shuffle(&s0, &s1);

    while (blocks--) {
        /* Remember old state */
        so0 = s0;
        so1 = s1;

        /* Load data and transform */
        m0 = Load(chunk);
        QuadRoundm(&s0, &s1, m0, 0xe9b5dba5b5c0fbcfull, 0x71374491428a2f98ull);
        m1 = Load(chunk + 16);
        QuadRoundm(&s0, &s1, m1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
        ShiftMessageA(&m0, m1);
        m2 = Load(chunk + 32);
        QuadRoundm(&s0, &s1, m2, 0x550c7dc3243185beull, 0x12835b01d807aa98ull);
        ShiftMessageA(&m1, m2);
        m3 = Load(chunk + 48);
        QuadRoundm(&s0, &s1, m3, 0xc19bf1749bdc06a7ull, 0x80deb1fe72be5d74ull);
        ShiftMessageB(&m2, m3, &m0);
        QuadRoundm(&s0, &s1, m0, 0x240ca1cc0fc19dc6ull, 0xefbe4786E49b69c1ull);
        ShiftMessageB(&m3, m0, &m1);
        QuadRoundm(&s0, &s1, m1, 0x76f988da5cb0a9dcull, 0x4a7484aa2de92c6full);
        ShiftMessageB(&m0, m1, &m2);
        QuadRoundm(&s0, &s1, m2, 0xbf597fc7b00327c8ull, 0xa831c66d983e5152ull);
        ShiftMessageB(&m1, m2, &m3);
        QuadRoundm(&s0, &s1, m3, 0x1429296706ca6351ull, 0xd5a79147c6e00bf3ull);
        ShiftMessageB(&m2, m3, &m0);
        QuadRoundm(&s0, &s1, m0, 0x53380d134d2c6dfcull, 0x2e1b213827b70a85ull);
        ShiftMessageB(&m3, m0, &m1);
        QuadRoundm(&s0, &s1, m1, 0x92722c8581c2c92eull, 0x766a0abb650a7354ull);
        ShiftMessageB(&m0, m1, &m2);
        QuadRoundm(&s0, &s1, m2, 0xc76c51A3c24b8b70ull, 0xa81a664ba2bfe8a1ull);
        ShiftMessageB(&m1, m2, &m3);
        QuadRoundm(&s0, &s1, m3, 0x106aa070f40e3585ull, 0xd6990624d192e819ull);
        ShiftMessageB(&m2, m3, &m0);
        QuadRoundm(&s0, &s1, m0, 0x34b0bcb52748774cull, 0x1e376c0819a4c116ull);
        ShiftMessageB(&m3, m0, &m1);
        QuadRoundm(&s0, &s1, m1, 0x682e6ff35b9cca4full, 0x4ed8aa4a391c0cb3ull);
        ShiftMessageC(&m0, m1, &m2);
        QuadRoundm(&s0, &s1, m2, 0x8cc7020884c87814ull, 0x78a5636f748f82eeull);
        ShiftMessageC(&m1, m2, &m3);
        QuadRoundm(&s0, &s1, m3, 0xc67178f2bef9A3f7ull, 0xa4506ceb90befffaull);

        /* Combine with old state */
        s0 = _mm_add_epi32(s0, so0);
        s1 = _mm_add_epi32(s1, so1);

        /* Advance */
        chunk += 64;
    }

    Unshuffle(&s0, &s1);
    _mm_storeu_si128((__m128i*)s, s0);
    _mm_storeu_si128((__m128i*)(s + 4), s1);
}

void Transform_2way(unsigned char* out, const unsigned char* in)
{
    __m128i am0, am1, am2, am3, as0, as1, aso0, aso1;
    __m128i bm0, bm1, bm2, bm3, bs0, bs1, bso0, bso1;

    /* Transform 1 */
    bs0 = as0 = INIT0;
    bs1 = as1 = INIT1;
    am0 = Load(in);
    bm0 = Load(in + 64);
    QuadRoundm(&as0, &as1, am0, 0xe9b5dba5b5c0fbcfull, 0x71374491428a2f98ull);
    QuadRoundm(&bs0, &bs1, bm0, 0xe9b5dba5b5c0fbcfull, 0x71374491428a2f98ull);
    am1 = Load(in + 16);
    bm1 = Load(in + 80);
    QuadRoundm(&as0, &as1, am1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    QuadRoundm(&bs0, &bs1, bm1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    ShiftMessageA(&am0, am1);
    ShiftMessageA(&bm0, bm1);
    am2 = Load(in + 32);
    bm2 = Load(in + 96);
    QuadRoundm(&as0, &as1, am2, 0x550c7dc3243185beull, 0x12835b01d807aa98ull);
    QuadRoundm(&bs0, &bs1, bm2, 0x550c7dc3243185beull, 0x12835b01d807aa98ull);
    ShiftMessageA(&am1, am2);
    ShiftMessageA(&bm1, bm2);
    am3 = Load(in + 48);
    bm3 = Load(in + 112);
    QuadRoundm(&as0, &as1, am3, 0xc19bf1749bdc06a7ull, 0x80deb1fe72be5d74ull);
    QuadRoundm(&bs0, &bs1, bm3, 0xc19bf1749bdc06a7ull, 0x80deb1fe72be5d74ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x240ca1cc0fc19dc6ull, 0xefbe4786E49b69c1ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x240ca1cc0fc19dc6ull, 0xefbe4786E49b69c1ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x76f988da5cb0a9dcull, 0x4a7484aa2de92c6full);
    QuadRoundm(&bs0, &bs1, bm1, 0x76f988da5cb0a9dcull, 0x4a7484aa2de92c6full);
    ShiftMessageB(&am0, am1, &am2);
    ShiftMessageB(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0xbf597fc7b00327c8ull, 0xa831c66d983e5152ull);
    QuadRoundm(&bs0, &bs1, bm2, 0xbf597fc7b00327c8ull, 0xa831c66d983e5152ull);
    ShiftMessageB(&am1, am2, &am3);
    ShiftMessageB(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0x1429296706ca6351ull, 0xd5a79147c6e00bf3ull);
    QuadRoundm(&bs0, &bs1, bm3, 0x1429296706ca6351ull, 0xd5a79147c6e00bf3ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x53380d134d2c6dfcull, 0x2e1b213827b70a85ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x53380d134d2c6dfcull, 0x2e1b213827b70a85ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x92722c8581c2c92eull, 0x766a0abb650a7354ull);
    QuadRoundm(&bs0, &bs1, bm1, 0x92722c8581c2c92eull, 0x766a0abb650a7354ull);
    ShiftMessageB(&am0, am1, &am2);
    ShiftMessageB(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0xc76c51A3c24b8b70ull, 0xa81a664ba2bfe8a1ull);
    QuadRoundm(&bs0, &bs1, bm2, 0xc76c51A3c24b8b70ull, 0xa81a664ba2bfe8a1ull);
    ShiftMessageB(&am1, am2, &am3);
    ShiftMessageB(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0x106aa070f40e3585ull, 0xd6990624d192e819ull);
    QuadRoundm(&bs0, &bs1, bm3, 0x106aa070f40e3585ull, 0xd6990624d192e819ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x34b0bcb52748774cull, 0x1e376c0819a4c116ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x34b0bcb52748774cull, 0x1e376c0819a4c116ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x682e6ff35b9cca4full, 0x4ed8aa4a391c0cb3ull);
    QuadRoundm(&bs0, &bs1, bm1, 0x682e6ff35b9cca4full, 0x4ed8aa4a391c0cb3ull);
    ShiftMessageC(&am0, am1, &am2);
    ShiftMessageC(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0x8cc7020884c87814ull, 0x78a5636f748f82eeull);
    QuadRoundm(&bs0, &bs1, bm2, 0x8cc7020884c87814ull, 0x78a5636f748f82eeull);
    ShiftMessageC(&am1, am2, &am3);
    ShiftMessageC(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0xc67178f2bef9A3f7ull, 0xa4506ceb90befffaull);
    QuadRoundm(&bs0, &bs1, bm3, 0xc67178f2bef9A3f7ull, 0xa4506ceb90befffaull);
    as0 = _mm_add_epi32(as0, INIT0);
    bs0 = _mm_add_epi32(bs0, INIT0);
    as1 = _mm_add_epi32(as1, INIT1);
    bs1 = _mm_add_epi32(bs1, INIT1);

    /* Transform 2 */
    aso0 = as0;
    bso0 = bs0;
    aso1 = as1;
    bso1 = bs1;
    QuadRound(&as0, &as1, 0xe9b5dba5b5c0fbcfull, 0x71374491c28a2f98ull);
    QuadRound(&bs0, &bs1, 0xe9b5dba5b5c0fbcfull, 0x71374491c28a2f98ull);
    QuadRound(&as0, &as1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    QuadRound(&bs0, &bs1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    QuadRound(&as0, &as1, 0x550c7dc3243185beull, 0x12835b01d807aa98ull);
    QuadRound(&bs0, &bs1, 0x550c7dc3243185beull, 0x12835b01d807aa98ull);
    QuadRound(&as0, &as1, 0xc19bf3749bdc06a7ull, 0x80deb1fe72be5d74ull);
    QuadRound(&bs0, &bs1, 0xc19bf3749bdc06a7ull, 0x80deb1fe72be5d74ull);
    QuadRound(&as0, &as1, 0x240cf2540fe1edc6ull, 0xf0fe4786649b69c1ull);
    QuadRound(&bs0, &bs1, 0x240cf2540fe1edc6ull, 0xf0fe4786649b69c1ull);
    QuadRound(&as0, &as1, 0x16f988fa61b9411eull, 0x6cc984be4fe9346full);
    QuadRound(&bs0, &bs1, 0x16f988fa61b9411eull, 0x6cc984be4fe9346full);
    QuadRound(&as0, &as1, 0xb9d99ec7b019fc65ull, 0xa88e5a6df2c65152ull);
    QuadRound(&bs0, &bs1, 0xb9d99ec7b019fc65ull, 0xa88e5a6df2c65152ull);
    QuadRound(&as0, &as1, 0xc7353eb0fdb1232bull, 0xe70eeaa09a1231c3ull);
    QuadRound(&bs0, &bs1, 0xc7353eb0fdb1232bull, 0xe70eeaa09a1231c3ull);
    QuadRound(&as0, &as1, 0xdc1eeefd5a0f118full, 0xcb976d5f3069bad5ull);
    QuadRound(&bs0, &bs1, 0xdc1eeefd5a0f118full, 0xcb976d5f3069bad5ull);
    QuadRound(&as0, &as1, 0xe15d5b1658f4ca9dull, 0xde0b7a040a35b689ull);
    QuadRound(&bs0, &bs1, 0xe15d5b1658f4ca9dull, 0xde0b7a040a35b689ull);
    QuadRound(&as0, &as1, 0x6fab9537a507ea32ull, 0x37088980007f3e86ull);
    QuadRound(&bs0, &bs1, 0x6fab9537a507ea32ull, 0x37088980007f3e86ull);
    QuadRound(&as0, &as1, 0xc0bbbe37cdaa3b6dull, 0x0d8cd6f117406110ull);
    QuadRound(&bs0, &bs1, 0xc0bbbe37cdaa3b6dull, 0x0d8cd6f117406110ull);
    QuadRound(&as0, &as1, 0x6fd15ca70b02e931ull, 0xdb48a36383613bdaull);
    QuadRound(&bs0, &bs1, 0x6fd15ca70b02e931ull, 0xdb48a36383613bdaull);
    QuadRound(&as0, &as1, 0x6d4378906ed41a95ull, 0x31338431521afacaull);
    QuadRound(&bs0, &bs1, 0x6d4378906ed41a95ull, 0x31338431521afacaull);
    QuadRound(&as0, &as1, 0x532fb63cb5c9a0e6ull, 0x9eccabbdc39c91f2ull);
    QuadRound(&bs0, &bs1, 0x532fb63cb5c9a0e6ull, 0x9eccabbdc39c91f2ull);
    QuadRound(&as0, &as1, 0x4c191d76a4954b68ull, 0x07237ea3d2c741c6ull);
    QuadRound(&bs0, &bs1, 0x4c191d76a4954b68ull, 0x07237ea3d2c741c6ull);
    as0 = _mm_add_epi32(as0, aso0);
    bs0 = _mm_add_epi32(bs0, bso0);
    as1 = _mm_add_epi32(as1, aso1);
    bs1 = _mm_add_epi32(bs1, bso1);

    /* Extract hash */
    Unshuffle(&as0, &as1);
    Unshuffle(&bs0, &bs1);
    am0 = as0;
    bm0 = bs0;
    am1 = as1;
    bm1 = bs1;

    /* Transform 3 */
    bs0 = as0 = INIT0;
    bs1 = as1 = INIT1;
    QuadRoundm(&as0, &as1, am0, 0xe9b5dba5B5c0fbcfull, 0x71374491428a2f98ull);
    QuadRoundm(&bs0, &bs1, bm0, 0xe9b5dba5B5c0fbcfull, 0x71374491428a2f98ull);
    QuadRoundm(&as0, &as1, am1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    QuadRoundm(&bs0, &bs1, bm1, 0xab1c5ed5923f82a4ull, 0x59f111f13956c25bull);
    ShiftMessageA(&am0, am1);
    ShiftMessageA(&bm0, bm1);
    bm2 = am2 = _mm_set_epi64x(0x0ull, 0x80000000ull);
    QuadRound(&as0, &as1, 0x550c7dc3243185beull, 0x12835b015807aa98ull);
    QuadRound(&bs0, &bs1, 0x550c7dc3243185beull, 0x12835b015807aa98ull);
    ShiftMessageA(&am1, am2);
    ShiftMessageA(&bm1, bm2);
    bm3 = am3 = _mm_set_epi64x(0x10000000000ull, 0x0ull);
    QuadRound(&as0, &as1, 0xc19bf2749bdc06a7ull, 0x80deb1fe72be5d74ull);
    QuadRound(&bs0, &bs1, 0xc19bf2749bdc06a7ull, 0x80deb1fe72be5d74ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x240ca1cc0fc19dc6ull, 0xefbe4786e49b69c1ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x240ca1cc0fc19dc6ull, 0xefbe4786e49b69c1ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x76f988da5cb0a9dcull, 0x4a7484aa2de92c6full);
    QuadRoundm(&bs0, &bs1, bm1, 0x76f988da5cb0a9dcull, 0x4a7484aa2de92c6full);
    ShiftMessageB(&am0, am1, &am2);
    ShiftMessageB(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0xbf597fc7b00327c8ull, 0xa831c66d983e5152ull);
    QuadRoundm(&bs0, &bs1, bm2, 0xbf597fc7b00327c8ull, 0xa831c66d983e5152ull);
    ShiftMessageB(&am1, am2, &am3);
    ShiftMessageB(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0x1429296706ca6351ull, 0xd5a79147c6e00bf3ull);
    QuadRoundm(&bs0, &bs1, bm3, 0x1429296706ca6351ull, 0xd5a79147c6e00bf3ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x53380d134d2c6dfcull, 0x2e1b213827b70a85ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x53380d134d2c6dfcull, 0x2e1b213827b70a85ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x92722c8581c2c92eull, 0x766a0abb650a7354ull);
    QuadRoundm(&bs0, &bs1, bm1, 0x92722c8581c2c92eull, 0x766a0abb650a7354ull);
    ShiftMessageB(&am0, am1, &am2);
    ShiftMessageB(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0xc76c51a3c24b8b70ull, 0xa81a664ba2bfe8A1ull);
    QuadRoundm(&bs0, &bs1, bm2, 0xc76c51a3c24b8b70ull, 0xa81a664ba2bfe8A1ull);
    ShiftMessageB(&am1, am2, &am3);
    ShiftMessageB(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0x106aa070f40e3585ull, 0xd6990624d192e819ull);
    QuadRoundm(&bs0, &bs1, bm3, 0x106aa070f40e3585ull, 0xd6990624d192e819ull);
    ShiftMessageB(&am2, am3, &am0);
    ShiftMessageB(&bm2, bm3, &bm0);
    QuadRoundm(&as0, &as1, am0, 0x34b0bcb52748774cull, 0x1e376c0819a4c116ull);
    QuadRoundm(&bs0, &bs1, bm0, 0x34b0bcb52748774cull, 0x1e376c0819a4c116ull);
    ShiftMessageB(&am3, am0, &am1);
    ShiftMessageB(&bm3, bm0, &bm1);
    QuadRoundm(&as0, &as1, am1, 0x682e6ff35b9cca4full, 0x4ed8aa4a391c0cb3ull);
    QuadRoundm(&bs0, &bs1, bm1, 0x682e6ff35b9cca4full, 0x4ed8aa4a391c0cb3ull);
    ShiftMessageC(&am0, am1, &am2);
    ShiftMessageC(&bm0, bm1, &bm2);
    QuadRoundm(&as0, &as1, am2, 0x8cc7020884c87814ull, 0x78a5636f748f82eeull);
    QuadRoundm(&bs0, &bs1, bm2, 0x8cc7020884c87814ull, 0x78a5636f748f82eeull);
    ShiftMessageC(&am1, am2, &am3);
    ShiftMessageC(&bm1, bm2, &bm3);
    QuadRoundm(&as0, &as1, am3, 0xc67178f2bef9a3f7ull, 0xa4506ceb90befffaull);
    QuadRoundm(&bs0, &bs1, bm3, 0xc67178f2bef9a3f7ull, 0xa4506ceb90befffaull);
    as0 = _mm_add_epi32(as0, INIT0);
    bs0 = _mm_add_epi32(bs0, INIT0);
    as1 = _mm_add_epi32(as1, INIT1);
    bs1 = _mm_add_epi32(bs1, INIT1);

    /* Extract hash into out */
    Unshuffle(&as0, &as1);
    Unshuffle(&bs0, &bs1);

    Save(out, as0);
    Save(out + 16, as1);
    Save(out + 32, bs0);
    Save(out + 48, bs1);
}

void InitG() {
  Initialize(sha256_g.s);
}

/* int main(int argc, char *argv[]) */
/* { */
/*   (void)argc; */
/*   (void)argv; */
/*   GlobalInit(); */
/*   struct sha256 sha256; */
/*   Initialize(sha256.s); */
/*   printf("sizeof %d\n", (int)sizeof(struct sha256)); */
/*   return 0; */
/* } */
