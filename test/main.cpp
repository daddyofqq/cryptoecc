// SPDX-License-Identifier: GPL-3.0-only
/*
 *  main.cpp
 *
 *  Copyright (C) 2020 Daniel Hu <daddy.of.qq@gmail.com>
 */

#include "debug.h"
#include "default_rng.h"
#include "ecc.h"
#include "mbedtls/bignum.h"
#include "mbedtls/ecdh.h"
#include "mbedtls/ecdsa.h"
#include "mbedtls/ecp.h"
#include "mbedtls/sha512.h"
#include "mpi.h"
#include "utility.h"
#include <iostream>
#include <string>
#include <vector>

using namespace mpi_ns;
using namespace utility;
using namespace ecc_ns;

static default_rng _default_rng;
mpi_ns::randomizer* rng = &_default_rng;

extern "C" int c_get_random(void* ctx, unsigned char* p, size_t size)
{
    (*rng)(p, size);
    return 0;
}

void print(ecc_ns::ops_stat stat)
{
    pr_debug("[I ", stat.I);
    pr_debug(",M ", stat.M);
    pr_debug(",A ", stat.A);
    pr_debug(",FA ", stat.FA);
    pr_debug(",R ", stat.R);
    pr_debug(",F ", stat.F, "]\n");
};

template <unsigned N>
void print(mpi<N> const& m)
{
    std::vector<uint8_t> o(m.num_of_bytes);
    m.output(o.data(), o.size());
    if (m.sign < 0) {
        pr_debug("-");
    }
    pr_debug(o, "\n");
}

template <unsigned S>
inline void print(FixedSizeBuf<S> const& buf)
{
    std::vector<uint8_t> o(buf.begin(), buf.end());
    pr_debug(o, "\n");
}

void Logger::log(const std::string& msg)
{
    std::cout << msg;
};

void Logger::log(const char* str)
{
    std::cout << std::string(str);
};

Logger logger{};

bool test_assert(std::string const& msg, bool cond)
{
    if (cond) {
        std::cout << "<" << msg << "> passed!\n";
        return true;
    } else {
        std::cout << "<" << msg << "> failed!\n";
        exit(1);
        return false;
    }
}

extern "C" void print_mpi(const mbedtls_mpi* X)
{
    std::vector<uint8_t> v((mbedtls_mpi_bitlen(X) + 7) / 8);
    mbedtls_mpi_write_binary(X, v.data(), v.size());
    pr_debug(v, "\n");
}

void test_mul()
{
    constexpr int N = 256;
    uint8_t a[N / 8];
    (*rng)(a, sizeof(a));

    mbedtls_mpi A, X;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&X);

    mbedtls_mpi_read_binary(&A, a, sizeof(a));
    measure(std::string("mbedtls multiplication"), [&]() { mbedtls_mpi_mul_mpi(&X, &A, &A); });

    std::vector<uint8_t> r(sizeof(a) * 2);
    mbedtls_mpi_write_binary(&X, r.data(), r.size());

    mpi<N> m{a, sizeof(a)};
    mpi<N * 2> p1;
    measure(std::string("my multiplication"), [&m, &p1]() { mul(m, m, p1); });
    test_assert("multiplication", p1.compare(mpi<p1.num_of_bits>(r.data(), r.size())) == 0);

    measure(std::string("square"), [&m, &p1]() { square(m, p1); });
    test_assert("square", p1.compare(mpi<p1.num_of_bits>(r.data(), r.size())) == 0);

    mbedtls_mpi_free(&A);
    mbedtls_mpi_free(&X);
}

void test_perf()
{
    constexpr int N = 256;
    uint8_t a[N / 8];
    (*rng)(a, sizeof(a));

    uint8_t b[N / 8];
    (*rng)(b, sizeof(a));

    mpi<N> A(a, sizeof(a));
    mpi<N> B(b, sizeof(b));

    mpi<N * 2> C;
    measure(std::string("mul"), [&]() { mul(A, B, C); });
    measure(std::string("shift"), [&]() { C.lshift(1); });
    mpi<N> R;
    measure(std::string("add"), [&]() { R.add(B); });
    measure(std::string("sub"), [&]() { R.sub(A); });
}

std::vector<uint8_t> gen_prime_number(int bits)
{
    mbedtls_mpi X;
    mbedtls_mpi_init(&X);
    mbedtls_mpi_gen_prime(&X, bits, 0,
                          c_get_random,
                          NULL);
    std::vector<uint8_t> r((bits + 7) / 8);
    mbedtls_mpi_write_binary(&X, r.data(), r.size());
    mbedtls_mpi_free(&X);
    return r;
};

void test_gcd()
{
    constexpr int N = 256;
    auto p1 = gen_prime_number(N);
    auto p2 = gen_prime_number(N);
    auto p3 = gen_prime_number(N);

    auto a = mpi<N>(p1.data(), p1.size()) * mpi<N>(p2.data(), p2.size());
    auto b = mpi<N>(p2.data(), p2.size()) * mpi<N>(p3.data(), p3.size());
    mpi<N * 2> g;
    measure(std::string("my gcd"), [&]() { gcd(a, b, g); });
    test_assert("gcd", g.compare(mpi<N>(p2.data(), p2.size())) == 0);
}

void test_div()
{
    constexpr int N = 64;
    std::vector<uint8_t> b(N / 8);
    (*rng)(b.data(), b.size());

    std::vector<uint8_t> a(N * 4 / 8);
    (*rng)(a.data(), a.size());

    mbedtls_mpi A, B, Q, R;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&Q);
    mbedtls_mpi_init(&R);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&B, b.data(), b.size());

    measure(std::string("mbed div"), [&]() { mbedtls_mpi_div_mpi(&Q, &R, &A, &B); });

    std::vector<uint8_t> r(N / 8);
    mbedtls_mpi_write_binary(&R, r.data(), r.size());

    std::vector<uint8_t> q(a.size());
    mbedtls_mpi_write_binary(&Q, q.data(), q.size());

    mpi<N> Y(b.data(), b.size());
    mpi<N * 4> X(a.data(), a.size());
    mpi<X.num_of_bits> _q;
    mpi<N> _r;
    measure(std::string("my div"), [&]() { div(X, Y, _q, _r); });

    test_assert("div q", _q.compare(mpi<N * 4>(q.data(), q.size())) == 0);
    test_assert("div r", _r.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_free(&A);
    mbedtls_mpi_free(&B);
    mbedtls_mpi_free(&Q);
    mbedtls_mpi_free(&R);
}

bool test_barret_mod()
{
    constexpr unsigned N = 256;
    auto m = gen_prime_number(N);
    std::vector<uint8_t> a(N * 2 / 8);
    (*rng)(a.data(), a.size());

    mbedtls_mpi A, B, R;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&R);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&B, m.data(), m.size());
    mbedtls_mpi_mod_mpi(&R, &A, &B);

    std::vector<uint8_t> r(m.size());
    mbedtls_mpi_write_binary(&R, r.data(), r.size());

    mpi<N> mm(m.data(), m.size());
    mpi<N * 2> ma(a.data(), a.size());
    mpi<N> mr;
    barrett_helper bh(mm);
    mpi<N> br;
    bh.reduce(ma, br);
    auto ret = (br.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&R);
    return ret;
}

void test_mod()
{
    constexpr unsigned N = 256;
    auto m = gen_prime_number(N);
    std::vector<uint8_t> a(N * 2 / 8);
    (*rng)(a.data(), a.size());

    mbedtls_mpi A, B, R;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&R);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&B, m.data(), m.size());
    measure(std::string("mbedtls reduction"), [&]() { mbedtls_mpi_mod_mpi(&R, &A, &B); });

    std::vector<uint8_t> r(m.size());
    mbedtls_mpi_write_binary(&R, r.data(), r.size());
    //pr_debug("mbedtls mod\n", r, "\n");

    mpi<N> mm(m.data(), m.size());
    mpi<N * 2> ma(a.data(), a.size());
    mpi<N> mr;
    measure(std::string("my reduction"), [&]() { fp_reduce(ma, mm, mr); });
    //pr_debug("modulus:\n");
    //print(mr);

    test_assert("modulus", mr.compare(mpi<N>(r.data(), r.size())) == 0);

    barrett_helper bh(mm);
    mpi<N> br;
    measure("barret modlus", [&]() { bh.reduce(ma, br); });
    test_assert("barret modulus", br.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&R);

    test_run(std::string("test barret modulus"),
             [&]() {
                 return test_barret_mod();
             });
}

void test_inverse()
{
    constexpr unsigned N = 256;
    auto m = gen_prime_number(N);
    std::vector<uint8_t> a(N / 8);
    (*rng)(a.data(), a.size());

    mbedtls_mpi X, A, M;
    mbedtls_mpi_init(&X);
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&M);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&M, m.data(), m.size());
    measure_slow(std::string("mbedtls inverse"), [&]() { mbedtls_mpi_inv_mod(&X, &A, &M); });

    std::vector<uint8_t> r(m.size());
    mbedtls_mpi_write_binary(&X, r.data(), r.size());
    //pr_debug("mbedtls inverse\n", "a: ", a, "\nm: ", m, "\nr: ", r, "\n");

    mpi<N> mm(m.data(), m.size());
    mpi<N> ma(a.data(), a.size());
    mpi<N> mr;

    // test montgomery inverse
    montg_inverse_helper mih(mm);
    measure_slow(std::string("montgomery inverse"), [&]() { mih(ma, mr); });
    test_assert("montgomery inverse", mr.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_init(&X);
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&M);
}

void test_naf()
{
    limb_t x = 1122334455;
    static int NAF2[] = {1, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 1, 0, -1, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, -1};

    mpi<bits_per_limb> m(x);
    NAF naf2(m);
    auto it1 = naf2.get_iterator();

    unsigned i = 0;
    while (!it1.end()) {
        int v = it1;
        if (v != NAF2[i])
            break;
        it1.next();
        i++;
    }
    test_assert(std::string("NAF2"), i == sizeof(NAF2) / sizeof(int));

    NAF<bits_per_limb, 6> naf6(m);
    auto it2 = naf6.get_iterator();

    static int NAF6[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 23, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, -9, 0, 0, 0, 0, 0, 0, 0, -9};
    i = 0;
    while (!it2.end()) {
        int v = it2;
        if (v != NAF6[i])
            break;
        it2.next();
        i++;
    }
    test_assert(std::string("NAF6"), i == sizeof(NAF6) / sizeof(int));
}

void test_mod_exp()
{
    constexpr unsigned N = 128;
    auto m = gen_prime_number(N);
    std::vector<uint8_t> a(N / 8);
    (*rng)(a.data(), a.size());
    auto e = gen_prime_number(N);

    mbedtls_mpi X, A, E, M;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&X);
    mbedtls_mpi_init(&E);
    mbedtls_mpi_init(&M);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&E, e.data(), e.size());
    mbedtls_mpi_read_binary(&M, m.data(), m.size());

    measure(std::string("mbedtls mod exp"), [&]() { mbedtls_mpi_exp_mod(&X, &A, &E, &M, NULL); });

    std::vector<uint8_t> r(m.size());
    mbedtls_mpi_write_binary(&X, r.data(), r.size());
    //pr_debug("a: ", a, "\nm: ", m, "\ne: ", e, "\nr: ", r, "\n");

    mpi<N> mm(m.data(), m.size());
    mpi<N> ma(a.data(), a.size());
    mpi<N> me(e.data(), e.size());
    mpi<N> mr;
    measure(std::string("my mod exp"), [&]() { fp_exp_mod(ma, me, mm, mr); });

    //print(mr);
    test_assert("mod exponent", mr.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_free(&A);
    mbedtls_mpi_free(&X);
    mbedtls_mpi_free(&E);
    mbedtls_mpi_free(&M);
}

extern "C" void mbedtls_mont_init(mbedtls_mpi_uint* mm, const mbedtls_mpi* N);
extern "C" void mbedtls_mont_mul(mbedtls_mpi* A, const mbedtls_mpi* B,
                                 const mbedtls_mpi* N, mbedtls_mpi_uint mm,
                                 const mbedtls_mpi* T);
void test_montg()
{
    constexpr unsigned N = 256;
    auto m = gen_prime_number(N);
    std::vector<uint8_t> a(N / 8);
    (*rng)(a.data(), a.size());
    std::vector<uint8_t> b(N / 8);
    (*rng)(b.data(), b.size());
    std::vector<uint8_t> t(N / 8 * 2);

    mbedtls_mpi A, B, M, T;
    mbedtls_mpi_init(&A);
    mbedtls_mpi_init(&B);
    mbedtls_mpi_init(&M);
    mbedtls_mpi_init(&T);

    mbedtls_mpi_read_binary(&A, a.data(), a.size());
    mbedtls_mpi_read_binary(&B, b.data(), b.size());
    mbedtls_mpi_read_binary(&M, m.data(), m.size());
    mbedtls_mpi_grow(&T, (m.size() + 1) * 2);

    mbedtls_mpi_uint mm;
    mbedtls_mont_init(&mm, &M);
    //printf("mbedtls mm %08lx\n", (unsigned long)mm);
    mbedtls_mont_mul(&A, &B, &M, mm, &T);

    std::vector<uint8_t> r(m.size());
    mbedtls_mpi_write_binary(&A, r.data(), r.size());
    //pr_debug("mbedtls montg mul\n", r, "\n");

    mpi<N> module(m.data(), m.size());
    mpi<N> ma(a.data(), a.size());
    mpi<N> mb(b.data(), b.size());
    mpi<N> mr;

    montg_helper mh(module);
    //printf("my mm %08lx\n", (unsigned long)mh.get_mm());

    mh.montred(ma, mb, mr);
    //pr_debug("montg mul resut -> ");
    //print(mr);

    test_assert(std::string("montgomery multiplication"), mr.compare(mpi<N>(r.data(), r.size())) == 0);

    mbedtls_mpi_free(&A);
    mbedtls_mpi_free(&B);
    mbedtls_mpi_free(&M);
    mbedtls_mpi_free(&T);
}

static void mbedtls_point_print(const mbedtls_ecp_point* p)
{
    pr_debug("mbedtls point:\n");
    pr_debug("X -> ");
    print_mpi(&p->X);
    pr_debug("Y -> ");
    print_mpi(&p->Y);
    pr_debug("Z -> ");
    print_mpi(&p->Z);
}

bool test_curve25519_scalar_mul()
{
    curve25519 c(rng);

    curve25519::keypair key;
    if (!c.genkeypair(key))
        return false;

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_CURVE25519);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    std::vector<uint8_t> prv(32);
    key.d.output(prv.data(), prv.size());
    mbedtls_mpi_read_binary(&m, prv.data(), prv.size());
    mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL);
    std::vector<uint8_t> x(32);
    mbedtls_mpi_write_binary(&R.X, x.data(), x.size());
    mpi<256> mpix(x.data(), x.size());

    ecp_aff<256> QX;
    c.scalar_mul(key.d, c.getG(), QX);
    print(c.get_statis());
    auto ret = (QX.X.compare(mpix) == 0);

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

bool test_curve448_scalar_mul()
{
    curve448 c(rng);

    curve448::keypair key;
    if (!c.genkeypair(key))
        return false;

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_CURVE448);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    std::vector<uint8_t> prv(56);
    key.d.output(prv.data(), prv.size());
    mbedtls_mpi_read_binary(&m, prv.data(), prv.size());
    mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL);
    std::vector<uint8_t> x(56);
    mbedtls_mpi_write_binary(&R.X, x.data(), x.size());
    mpi<448> mpix(x.data(), x.size());

    ecp_aff<448> QX;
    c.scalar_mul(key.d, c.getG(), QX);
    print(c.get_statis());
    auto ret = (QX.X.compare(mpix) == 0);

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

void test_curve448()
{
    std::vector<uint8_t> a(896 / 8);
    (*rng)(a.data(), a.size());
    curve448 c(rng);
    mpi<896> ma(a.data(), a.size());
    mpi<448> r1;
    measure(std::string("benchmark raw reduction"), [&]() { fp_reduce(ma, c.dp_p(), r1); });

    mpi<448> r2;
    measure(std::string("curve448 reduction"), [&]() { c.reduce(ma, r2); });
    test_assert("curve448 reduction", r1.compare(r2) == 0);

    curve448::keypair key;
    auto ret = c.genkeypair(key);
    test_assert(std::string("curve448 key generation"), ret);

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_CURVE448);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    std::vector<uint8_t> prv(56);
    key.d.output(prv.data(), prv.size());
    mbedtls_mpi_read_binary(&m, prv.data(), prv.size());

    measure_very_slow(std::string("mbedtls curve448 scalar"), [&]() mutable { mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL); });

    ecp_aff<448> QX;
    measure_very_slow(std::string("curve448 mul"), [&]() { c.scalar_mul(key.d, c.getG(), QX); });

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    test_run("curve448 scalar multiplication",
             [&]() {
                 return test_curve448_scalar_mul();
             });
};

void test_curve25519()
{
    std::vector<uint8_t> a(512 / 8);
    (*rng)(a.data(), a.size());
    curve25519 c(rng);
    mpi<512> ma(a.data(), a.size());
    mpi<256> r1;
    measure(std::string("benchmark raw reduction"), [&]() { fp_reduce(ma, c.dp_p(), r1); });

    mpi<256> r2;
    measure(std::string("curve25519 p-255 reduction"), [&]() { c.reduce(ma, r2); });
    test_assert("p-255 reduction", r1.compare(r2) == 0);

    curve25519::keypair key;
    auto ret = c.genkeypair(key);
    test_assert(std::string("curve25519 key generation"), ret);

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_CURVE25519);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    std::vector<uint8_t> prv(32);
    key.d.output(prv.data(), prv.size());
    mbedtls_mpi_read_binary(&m, prv.data(), prv.size());

    measure_very_slow(std::string("mbedtls curve25519 scalar"), [&]() { mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL); });
    ecp_aff<256> QX;
    measure_very_slow(std::string("curve25519 mul"), [&]() { c.scalar_mul(key.d, c.getG(), QX); });

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    test_run("curve25519 scalar multiplication",
             [&]() {
                 return test_curve25519_scalar_mul();
             });
};

bool test_short_weierstrass_scalar_mul_fixed_base()
{
    const int N = 256;
    std::vector<uint8_t> scalar(N / 8);
    (*rng)(scalar.data(), scalar.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    mbedtls_mpi_read_binary(&m, scalar.data(), scalar.size());
    mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL);

    std::vector<uint8_t> x(32);
    mbedtls_mpi_write_binary(&R.X, x.data(), x.size());
    std::vector<uint8_t> y(32);
    mbedtls_mpi_write_binary(&R.Y, y.data(), y.size());
    mpi<256> mpix(x.data(), x.size());
    mpi<256> mpiy(y.data(), y.size());

    secp256r1 c(rng);
    mpi<N> k(scalar.data(), scalar.size());
    ecp_aff<256> affp;
    c.scalar_mul(k, c.getG(), affp);
    print(c.get_statis());
    auto ret = (affp.X.compare(mpix) == 0 &&
                affp.Y.compare(mpiy) == 0);

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

bool test_short_weierstrass_scalar_mul_unknown_base()
{
    const int N = 256;
    std::vector<uint8_t> scalar(N / 8);
    (*rng)(scalar.data(), scalar.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    mbedtls_mpi_read_binary(&m, scalar.data(), scalar.size());
    mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL);

    std::vector<uint8_t> x(32);
    mbedtls_mpi_write_binary(&R.X, x.data(), x.size());
    std::vector<uint8_t> y(32);
    mbedtls_mpi_write_binary(&R.Y, y.data(), y.size());
    mpi<256> mpix(x.data(), x.size());
    mpi<256> mpiy(y.data(), y.size());

    secp256r1 c(rng);
    mpi<N> k(scalar.data(), scalar.size());

    // test the multiplication for unknown base
    auto affg = c.getG();
    ecp_aff<256> affp;
    c.scalar_mul(k, affg, affp);
    print(c.get_statis());
    auto ret = (affp.X.compare(mpix) == 0 &&
                affp.Y.compare(mpiy) == 0);

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

void test_short_weierstrass()
{
    const int N = 256;
    std::vector<uint8_t> scalar(N / 8);
    (*rng)(scalar.data(), scalar.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point R;
    mbedtls_ecp_point_init(&R);
    mbedtls_mpi m;
    mbedtls_mpi_init(&m);
    mbedtls_mpi_read_binary(&m, scalar.data(), scalar.size());
    measure_very_slow(std::string("mbedtls scalar multiplication"), [&]() { mbedtls_ecp_mul(&grp, &R, &m, &grp.G, NULL, NULL); });

    secp256r1 c(rng);
    mpi<N> k(scalar.data(), scalar.size());
    ecp_aff<256> affp;
    measure_very_slow(std::string("fixed comb point multiplication"), [&]() { c.scalar_mul(k, affp); });

    auto affg = c.getG();
    ecp_aff<256> affp2;
    measure_very_slow(std::string("unknown base multiplication"), [&]() { c.scalar_mul(k, affg, affp2); });

    mbedtls_mpi_free(&m);
    mbedtls_ecp_point_free(&R);
    mbedtls_ecp_group_free(&grp);

    test_run(std::string("short weierstrass scalar multiplication - fixed base"),
             [&]() {
                 return test_short_weierstrass_scalar_mul_fixed_base();
             });

    test_run(std::string("short weierstrass scalar multiplication - unknown base"),
             [&]() {
                 return test_short_weierstrass_scalar_mul_unknown_base();
             });
}

bool test_ecdsa_verify()
{
    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point Q;
    mbedtls_ecp_point_init(&Q);
    mbedtls_mpi r, s, d;
    mbedtls_mpi_init(&r);
    mbedtls_mpi_init(&s);
    mbedtls_mpi_init(&d);
    std::vector<uint8_t> buf(32);
    (*rng)(buf.data(), buf.size());

    mbedtls_ecp_gen_keypair(&grp, &d, &Q, c_get_random, NULL);
    std::vector<uint8_t> vx(32);
    std::vector<uint8_t> vy(32);
    mbedtls_mpi_write_binary(&Q.X, vx.data(), vx.size());
    mbedtls_mpi_write_binary(&Q.Y, vy.data(), vy.size());

    mbedtls_ecdsa_sign(&grp, &r, &s,
                       &d, buf.data(), buf.size(),
                       c_get_random, NULL);

    std::vector<uint8_t> vr(32);
    std::vector<uint8_t> vs(32);
    mbedtls_mpi_write_binary(&r, vr.data(), vr.size());
    mbedtls_mpi_write_binary(&s, vs.data(), vs.size());

    mpi<256> mm(buf.data(), buf.size());
    ecp_aff<256> mq(mpi<256>(vx.data(), vx.size()), mpi<256>(vy.data(), vy.size()));
    mpi<256> mr(vr.data(), vr.size());
    mpi<256> ms(vs.data(), vs.size());
    secp256r1 c(rng);

    auto ret = c.ecdsa_verify(mq, mm, secp256r1::signature_t{mr, ms});

    mbedtls_mpi_free(&r);
    mbedtls_mpi_free(&s);
    mbedtls_mpi_free(&d);
    mbedtls_ecp_point_free(&Q);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

void measure_ecdsa_verify()
{
    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point Q;
    mbedtls_ecp_point_init(&Q);
    mbedtls_mpi r, s, d;
    mbedtls_mpi_init(&r);
    mbedtls_mpi_init(&s);
    mbedtls_mpi_init(&d);
    std::vector<uint8_t> buf(32);
    (*rng)(buf.data(), buf.size());

    mbedtls_ecp_gen_keypair(&grp, &d, &Q, c_get_random, NULL);
    std::vector<uint8_t> vx(32);
    std::vector<uint8_t> vy(32);
    mbedtls_mpi_write_binary(&Q.X, vx.data(), vx.size());
    mbedtls_mpi_write_binary(&Q.Y, vy.data(), vy.size());

    measure_very_slow(std::string("mbedtls sign"),
                      [&]() {
                          mbedtls_ecdsa_sign(&grp, &r, &s,
                                             &d, buf.data(), buf.size(),
                                             c_get_random, NULL);
                      });

    std::vector<uint8_t> vr(32);
    std::vector<uint8_t> vs(32);
    mbedtls_mpi_write_binary(&r, vr.data(), vr.size());
    mbedtls_mpi_write_binary(&s, vs.data(), vs.size());

    mpi<256> mm(buf.data(), buf.size());
    ecp_aff<256> mq(mpi<256>(vx.data(), vx.size()), mpi<256>(vy.data(), vy.size()));
    mpi<256> mr(vr.data(), vr.size());
    mpi<256> ms(vs.data(), vs.size());
    secp256r1 c(rng);

    measure_very_slow(std::string("ecdsa verify"),
                      [&]() {
                          c.ecdsa_verify(mq, mm, secp256r1::signature_t{mr, ms});
                      });

    mbedtls_mpi_free(&r);
    mbedtls_mpi_free(&s);
    mbedtls_mpi_free(&d);
    mbedtls_ecp_point_free(&Q);
    mbedtls_ecp_group_free(&grp);
}

bool test_ecdsa_sign()
{
    std::vector<uint8_t> buf(32);
    (*rng)(buf.data(), buf.size());

    secp256r1 c(rng);
    secp256r1::keypair key;
    c.genkeypair(key);
    secp256r1::signature_t sig;
    mpi<256> mm(buf.data(), buf.size());

    c.ecdsa_sign(key.d, mm, sig);

    std::vector<uint8_t> vx(32);
    std::vector<uint8_t> vy(32);
    std::vector<uint8_t> vr(32);
    std::vector<uint8_t> vs(32);
    std::vector<uint8_t> vd(32);

    sig.s.output(vs.data(), vs.size());
    sig.r.output(vr.data(), vr.size());
    key.Q.X.output(vx.data(), vx.size());
    key.Q.Y.output(vy.data(), vy.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point Q;
    mbedtls_ecp_point_init(&Q);
    mbedtls_mpi r, s;
    mbedtls_mpi_init(&r);
    mbedtls_mpi_init(&s);

    mbedtls_mpi_read_binary(&Q.X, vx.data(), vx.size());
    mbedtls_mpi_read_binary(&Q.Y, vy.data(), vy.size());
    mbedtls_mpi_lset(&Q.Z, 1);
    mbedtls_mpi_read_binary(&r, vr.data(), vr.size());
    mbedtls_mpi_read_binary(&s, vs.data(), vs.size());

    bool ret = (mbedtls_ecdsa_verify(&grp,
                                     buf.data(), buf.size(),
                                     &Q, &r, &s) == 0);

    mbedtls_mpi_free(&r);
    mbedtls_mpi_free(&s);
    mbedtls_ecp_point_free(&Q);
    mbedtls_ecp_group_free(&grp);

    return ret;
}

void measure_ecdsa_sign()
{
    std::vector<uint8_t> buf(32);
    (*rng)(buf.data(), buf.size());

    secp256r1 c(rng);
    secp256r1::keypair key;
    c.genkeypair(key);
    secp256r1::signature_t sig;
    mpi<256> mm(buf.data(), buf.size());

    measure_very_slow(std::string("ecdsa sign"),
                      [&]() {
                          c.ecdsa_sign(key.d, mm, sig);
                      });

    std::vector<uint8_t> vx(32);
    std::vector<uint8_t> vy(32);
    std::vector<uint8_t> vr(32);
    std::vector<uint8_t> vs(32);
    std::vector<uint8_t> vd(32);

    sig.s.output(vs.data(), vs.size());
    sig.r.output(vr.data(), vr.size());
    key.Q.X.output(vx.data(), vx.size());
    key.Q.Y.output(vy.data(), vy.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point Q;
    mbedtls_ecp_point_init(&Q);
    mbedtls_mpi r, s;
    mbedtls_mpi_init(&r);
    mbedtls_mpi_init(&s);

    mbedtls_mpi_read_binary(&Q.X, vx.data(), vx.size());
    mbedtls_mpi_read_binary(&Q.Y, vy.data(), vy.size());
    mbedtls_mpi_lset(&Q.Z, 1);
    mbedtls_mpi_read_binary(&r, vr.data(), vr.size());
    mbedtls_mpi_read_binary(&s, vs.data(), vs.size());

    measure_very_slow(std::string("mbedtls ecdsa verify"),
                      [&]() {
                          mbedtls_ecdsa_verify(&grp,
                                               buf.data(), buf.size(),
                                               &Q, &r, &s);
                      });

    mbedtls_mpi_free(&r);
    mbedtls_mpi_free(&s);
    mbedtls_ecp_point_free(&Q);
    mbedtls_ecp_group_free(&grp);
}

void test_ecdsa()
{
    test_run(std::string("ecdsa verify"),
             [&]() {
                 return test_ecdsa_verify();
             });

    measure_ecdsa_verify();

    test_run(std::string("ecdsa sign"),
             [&]() {
                 return test_ecdsa_sign();
             });

    measure_ecdsa_sign();
}

bool test_ecdh_compute_secret()
{
    secp256r1 c(rng);
    secp256r1::keypair mkey;
    c.genkeypair(mkey);

    std::vector<uint8_t> vx(32);
    std::vector<uint8_t> vy(32);
    mkey.Q.X.output(vx.data(), vx.size());
    mkey.Q.Y.output(vy.data(), vy.size());

    mbedtls_ecp_group grp;
    mbedtls_ecp_keypair key;
    mbedtls_ecp_keypair_init(&key);
    mbedtls_ecp_group_init(&grp);
    mbedtls_ecp_group_load(&grp, MBEDTLS_ECP_DP_SECP256R1);
    mbedtls_ecp_point Q;
    mbedtls_ecp_point Qp;
    mbedtls_ecp_point_init(&Q);
    mbedtls_ecp_point_init(&Qp);
    mbedtls_mpi d, secret;
    mbedtls_mpi_init(&d);
    mbedtls_mpi_init(&secret);

    mbedtls_ecp_gen_keypair(&grp, &d, &Q, c_get_random, NULL);
    mbedtls_mpi_read_binary(&Qp.X, vx.data(), vx.size());
    mbedtls_mpi_read_binary(&Qp.Y, vy.data(), vy.size());
    mbedtls_mpi_lset(&Qp.Z, 1);

    mbedtls_ecdh_compute_shared(&grp, &secret,
                                &Qp, &d,
                                c_get_random, NULL);

    std::vector<uint8_t> vs(32);
    mbedtls_mpi_write_binary(&secret, vs.data(), vs.size());
    mbedtls_mpi_write_binary(&Q.X, vx.data(), vx.size());
    mbedtls_mpi_write_binary(&Q.Y, vy.data(), vy.size());

    mpi<256> ms;
    auto ok = c.ecdh_compute_shared(mkey, ecp_aff<256>(mpi<256>(vx.data(), vx.size()), mpi<256>(vy.data(), vy.size())), ms);

    mbedtls_mpi_free(&d);
    mbedtls_mpi_free(&secret);
    mbedtls_ecp_point_free(&Q);
    mbedtls_ecp_point_free(&Qp);
    mbedtls_ecp_group_free(&grp);

    return ok && ms.compare(mpi<256>(vs.data(), vs.size())) == 0;
}

void test_ecdh()
{
    test_run(std::string("ecdh"),
             [&]() {
                 return test_ecdh_compute_secret();
             });
}

static void shuffle(uint8_t* p, size_t len)
{
    uint8_t* q = &p[len - 1];
    while (p < q) {
        auto tmp = *p;
        *p = *q;
        *q = tmp;
        p++;
        q--;
    }
}

void test_sha512()
{
    uint8_t msg[] = {0xaf, 0x82};
    std::vector<uint8_t> hash1(64);
    mbedtls_sha512(msg, sizeof(msg), hash1.data(), false);

    std::vector<uint8_t> hash2(64);
    default_sha512 sha{};
    sha.update(msg, sizeof(msg));
    sha.finish(hash2.data());

    test_assert("sha512", hash1 == hash2);
}

bool test_ed25519_kp()
{
    ed25519 c(rng);
    ed25519::keypair kp;
    c.genkeypair(kp);
    return c.verify_pk(kp.Q);
}

void test_ed25519()
{
    test_run(std::string("ed25519 keypair"),
             [&]() {
                 return test_ed25519_kp();
             });

    ed25519 c(rng);
    curve25519 cc(rng);
    auto G = c.getG();
    test_assert("edward verify G", c.verify_point(G));

    auto m = c.to_montgomery(G);
    test_assert("edward G conversion", cc.verify_point(m));

    ed25519::signature_t sig;
    // test vector from RFC8032 (Test 3)
    std::array<uint8_t, 32> dx{0xc5, 0xaa, 0x8d, 0xf4, 0x3f, 0x9f, 0x83, 0x7b,
                               0xed, 0xb7, 0x44, 0x2f, 0x31, 0xdc, 0xb7, 0xb1,
                               0x66, 0xd3, 0x85, 0x35, 0x07, 0x6f, 0x09, 0x4b,
                               0x85, 0xce, 0x3a, 0x2e, 0x0b, 0x44, 0x58, 0xf7};

    shuffle(dx.data(), dx.size());
    mpi<256> d(dx.data(), dx.size());
    std::array<uint8_t, 2> msg{0xaf, 0x82};
    measure_very_slow(std::string("eddsa sign"),
                      [&]() {
                          c.eddsa_sign(d, msg.data(), msg.size(), sig);
                      });

    std::vector<uint8_t> expected_sig{0x62, 0x91, 0xd6, 0x57, 0xde, 0xec, 0x24, 0x02,
                                      0x48, 0x27, 0xe6, 0x9c, 0x3a, 0xbe, 0x01, 0xa3,
                                      0x0c, 0xe5, 0x48, 0xa2, 0x84, 0x74, 0x3a, 0x44,
                                      0x5e, 0x36, 0x80, 0xd7, 0xdb, 0x5a, 0xc3, 0xac,
                                      0x18, 0xff, 0x9b, 0x53, 0x8d, 0x16, 0xf2, 0x90,
                                      0xae, 0x67, 0xf7, 0x60, 0x98, 0x4d, 0xc6, 0x59,
                                      0x4a, 0x7c, 0x15, 0xe9, 0x71, 0x6e, 0xd2, 0x8d,
                                      0xc0, 0x27, 0xbe, 0xce, 0xea, 0x1e, 0xc4, 0x0a};
    test_assert("ed25519 signing", std::vector<uint8_t>(sig.begin(), sig.end()) == expected_sig);

    auto pk = c.genpk(d);
    test_assert("ed25519 verify", c.eddsa_verify(pk, msg.data(), msg.size(), sig));
    measure_very_slow(std::string("eddsa verify"),
                      [&]() {
                          c.eddsa_verify(pk, msg.data(), msg.size(), sig);
                      });
}


void test_ed448()
{
    constexpr unsigned N = 448;
    std::vector<uint8_t> a(N * 2 / 8);
    (*rng)(a.data(), a.size());
    ed448 c(rng);
    mpi<N * 2> ma(a.data(), a.size());
    mpi<N> r1;
    measure(std::string("448 raw reduction"), [&]() { fp_reduce(ma, c.dp_p(), r1); });

    mpi<N> r2;
    measure(std::string("ed448 fast reduction"), [&]() { c.reduce(ma, r2); });

    ed448::signature_t sig;
    // test vector from RFC8032
    std::array<uint8_t, 57> dx{0xc4, 0xea, 0xb0, 0x5d, 0x35, 0x70, 0x07, 0xc6,
                               0x32, 0xf3, 0xdb, 0xb4, 0x84, 0x89, 0x92, 0x4d,
                               0x55, 0x2b, 0x08, 0xfe, 0x0c, 0x35, 0x3a, 0x0d,
                               0x4a, 0x1f, 0x00, 0xac, 0xda, 0x2c, 0x46, 0x3a,
                               0xfb, 0xea, 0x67, 0xc5, 0xe8, 0xd2, 0x87, 0x7c,
                               0x5e, 0x3b, 0xc3, 0x97, 0xa6, 0x59, 0x94, 0x9e,
                               0xf8, 0x02, 0x1e, 0x95, 0x4e, 0x0a, 0x12, 0x27,
                               0x4e};

    std::vector<uint8_t> pk_expected{0x43, 0xba, 0x28, 0xf4, 0x30, 0xcd, 0xff, 0x45,
                                     0x6a, 0xe5, 0x31, 0x54, 0x5f, 0x7e, 0xcd, 0x0a,
                                     0xc8, 0x34, 0xa5, 0x5d, 0x93, 0x58, 0xc0, 0x37,
                                     0x2b, 0xfa, 0x0c, 0x6c, 0x67, 0x98, 0xc0, 0x86,
                                     0x6a, 0xea, 0x01, 0xeb, 0x00, 0x74, 0x28, 0x02,
                                     0xb8, 0x43, 0x8e, 0xa4, 0xcb, 0x82, 0x16, 0x9c,
                                     0x23, 0x51, 0x60, 0x62, 0x7b, 0x4c, 0x3a, 0x94,
                                     0x80};

    shuffle(dx.data(), dx.size());
    mpi<57 * 8> d(dx.data(), dx.size());
    auto pk = c.genpk(d);
    test_assert(std::string("ed448 pk generation"), pk_expected == std::vector<uint8_t>(pk.begin(), pk.end()));

    std::vector<uint8_t> msg{0x03};
    std::vector<uint8_t> ctx{0x66, 0x6f, 0x6f};

    measure_very_slow(std::string("ed448 eddsa sign"),
                      [&]() {
                          c.eddsa_sign(d, msg.data(), msg.size(), sig, ctx.data(), ctx.size());
                      });

    std::vector<uint8_t> expected_sig{
        0xd4, 0xf8, 0xf6, 0x13, 0x17, 0x70, 0xdd, 0x46,
        0xf4, 0x08, 0x67, 0xd6, 0xfd, 0x5d, 0x50, 0x55,
        0xde, 0x43, 0x54, 0x1f, 0x8c, 0x5e, 0x35, 0xab,
        0xbc, 0xd0, 0x01, 0xb3, 0x2a, 0x89, 0xf7, 0xd2,
        0x15, 0x1f, 0x76, 0x47, 0xf1, 0x1d, 0x8c, 0xa2,
        0xae, 0x27, 0x9f, 0xb8, 0x42, 0xd6, 0x07, 0x21,
        0x7f, 0xce, 0x6e, 0x04, 0x2f, 0x68, 0x15, 0xea,
        0x00, 0x0c, 0x85, 0x74, 0x1d, 0xe5, 0xc8, 0xda,
        0x11, 0x44, 0xa6, 0xa1, 0xab, 0xa7, 0xf9, 0x6d,
        0xe4, 0x25, 0x05, 0xd7, 0xa7, 0x29, 0x85, 0x24,
        0xfd, 0xa5, 0x38, 0xfc, 0xcb, 0xbb, 0x75, 0x4f,
        0x57, 0x8c, 0x1c, 0xad, 0x10, 0xd5, 0x4d, 0x0d,
        0x54, 0x28, 0x40, 0x7e, 0x85, 0xdc, 0xbc, 0x98,
        0xa4, 0x91, 0x55, 0xc1, 0x37, 0x64, 0xe6, 0x6c,
        0x3c, 0x00};

    test_assert("ed448 signing", std::vector<uint8_t>(sig.begin(), sig.end()) == expected_sig);
    test_assert("ed448 verify", c.eddsa_verify(pk, msg.data(), msg.size(), sig, ctx.data(), ctx.size()));
    measure_very_slow(std::string("ed448 eddsa verify"),
                      [&]() {
                          c.eddsa_verify(pk, msg.data(), msg.size(), sig, ctx.data(), ctx.size());
                      });
}

void test_ecc_basic()
{
    std::vector<uint8_t> a(512 / 8);
    (*rng)(a.data(), a.size());

    secp256r1 c(rng);

    mpi<512> ma(a.data(), a.size());
    mpi<256> r1;
    measure(std::string("benchmark raw reduction"), [&]() { fp_reduce(ma, c.dp_p(), r1); });

    mpi<256> r2;
    measure(std::string("p-256 reduction"), [&]() { c.reduce(ma, r2); });
    test_assert("p-256 reduction", r1.compare(r2) == 0);
}

int main(int argc, char* argv[])
{
    test_perf();

    test_mul();

    test_gcd();

    test_div();

    test_mod();

    //test_montg();

    test_mod_exp();

    test_inverse();

    test_naf();

    test_ecc_basic();

    test_short_weierstrass();

    test_curve25519();

    test_curve448();

    test_ecdsa();

    test_ecdh();

    test_sha512();

    test_ed25519();

    test_ed448();
}
