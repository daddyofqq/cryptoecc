// SPDX-License-Identifier: GPL-3.0-only
/*
 *  ecc.h
 *
 *  Copyright (C) 2019 Daniel Hu <daddy.of.qq@gmail.com>
 *
 *  ECC crypto in mordern C++ (need c++17 or later)
 *
 */

/* The code is heavily templated to allow compile-time optimization
 *
 * Three types of curve supported
 * Short Weierstrass : y^2 = x^3 + ax + b mod p
 * Montgomery curve : y^2 = x^3 + A x^2 + x mod p
 * Twisted Edward Curve : a * x^2 + y^2 = 1 + d x^2 y^2
 *
 * algorithm for ECDSA, EDDSA, and ECDH are implemented
 *
 * The code support curve parameters for NIST p256, curve25519/curve448, and ed25519/ed448 (for above curve type respectively)
 * However, it should be quite easy to support other parameters as long as they
 * belong to to one of the above three curves
 *
 *
 * Reference:
 *  for overall ECC crypto : Guide to Elliptic Curve Cryptography
 *  for Edward curve: Fast Adddition and doubling on elliptic curves
 *               (Daniel J.Bernstein & Tanja Lange)
 *               and https://en.wikipedia.org/wiki/Twisted_Edwards_curve
 *               and RFC8032
 *
 */
#ifndef ECC_H
#define ECC_H

#include "hash.h"
#include "mpi.h"

#define ECC_STAT

namespace ecc_ns
{

using namespace mpi_ns;

enum class CurveType {
    NONE = 0,
    SHORT_WEIERSTRASS,
    MONTGOMERY,
    EDWARD
};

template <CurveType CT, typename DP, unsigned N>
struct ec_curve;

enum class CurveId {
    NONE = 0,
    SECP256R1,
    CURVE25519,
    CURVE448,
    ED25519,
    ED448
};

// an utility class to provide 32bit view into mpi limb,
// which may be 32- or 64-bit.
// it is used by the fast modulation of NIST curve
template <typename MPI>
struct mpi32bit_viewer {

    mpi32bit_viewer(MPI& m) : m(m){};

    inline uint32_t operator[](int index) const
    {
        static_assert(sizeof(uint32_t) == sizeof(limb_t) ||
                          sizeof(uint64_t) == sizeof(limb_t),
                      "your architecture is not supported");

        if constexpr (sizeof(uint32_t) == sizeof(limb_t)) {
            return m[index];
        } else {
            limb_t p = m[index >> 1];
            if (index & 1) {
                return static_cast<uint32_t>(p >> (bits_per_limb / 2));
            } else {
                return static_cast<uint32_t>(p);
            }
        }
    }

    template <unsigned IDX>
    inline uint32_t get() const
    {
        static_assert(sizeof(uint32_t) == sizeof(limb_t) ||
                          sizeof(uint64_t) == sizeof(limb_t),
                      "your architecture is not supported");

        if constexpr (sizeof(uint32_t) == sizeof(limb_t)) {
            return m[IDX];
        } else {
            limb_t p = m[IDX >> 1];
            if (IDX & 1) {
                return static_cast<uint32_t>(p >> (bits_per_limb / 2));
            } else {
                return static_cast<uint32_t>(p);
            }
        }
    }

    MPI& m;
};

template <typename MPI>
auto get_mpi32_viewer(MPI& m)
{
    return mpi32bit_viewer(m);
};

inline void do_assemble_mpi32(limb_t* p, uint32_t v)
{
    *p = static_cast<limb_t>(v);
}

template <typename... ARGS>
inline void do_assemble_mpi32(limb_t* p, uint32_t v1, uint32_t v2, ARGS... args)
{
    if constexpr (sizeof(limb_t) == sizeof(uint32_t)) {
        *p++ = v1;
        *p++ = v2;
    } else {
        *p++ = static_cast<limb_t>(v2) << (bits_per_limb / 2) | static_cast<limb_t>(v1);
    }

    if constexpr (sizeof...(args) > 0) {
        do_assemble_mpi32(p, args...);
    }
}

template <typename... ARGS>
auto assemble_mpi32(ARGS... args)
{
    mpi<sizeof...(args) * 32 + 1> m;

    do_assemble_mpi32(&m[0], args...);

    return m;
}

// Affine coordinates
template <unsigned N>
struct ecp_aff {
    ecp_aff() : X(1), Y(0){}; // ZERO by default
    ecp_aff(mpi<N> const& x, mpi<N> const& y) : X(x), Y(y){};
    ecp_aff(ecp_aff<N> const& point) : X(point.X), Y(point.Y){};
    ecp_aff(const uint8_t* x, size_t xlen, const uint8_t* y, size_t ylen) : X(x, xlen), Y(y, ylen){};
    inline bool is_zero() const
    {
        return Y.equal(0);
    }

    void set_zero()
    {
        X = mpi<N>(1);
        Y = mpi<N>(0);
    }

    ecp_aff& operator=(ecp_aff<N> const& other)
    {
        X = other.X;
        Y = other.Y;
        return *this;
    }

    mpi<N> X;
    mpi<N> Y;
};

// projective coordinates where (X, Y, Z) -> (X/Z, Y/Z)
template <unsigned N>
struct ecp_proj {
    ecp_proj() : X(0), Y(1), Z(1){};
    ecp_proj(mpi<N> const& x, mpi<N> const& y, mpi<N> const& z) : X(x), Y(y), Z(z){};
    ecp_proj(ecp_proj<N> const& point) : X(point.X), Y(point.Y), Z(point.Z){};
    ecp_proj(ecp_aff<N> const& point) : X(point.X), Y(point.Y), Z(1){};

    inline bool is_zero() const
    {
        return X.equal(0);
    }

    void set_zero()
    {
        X = mpi<N>(0);
        Y = mpi<N>(1);
        Z = mpi<N>(1);
    }

    ecp_proj& operator=(ecp_proj<N> const& other)
    {
        X = other.X;
        Y = other.Y;
        Z = other.Z;
        return *this;
    }

    mpi<N> X;
    mpi<N> Y;
    mpi<N> Z;
};

// extended homogenious coordinates here (X, Y, Z, T)
// where x = X/Z, y = Y/Z, x * y = T/Z
template <unsigned N>
struct ecp_ext {
    ecp_ext() : X(0), Y(1), Z(1), T(0){};

    inline bool is_zero() const
    {
        return X.equal(0);
    }

    void set_zero()
    {
        X = mpi<N>(0);
        Y = mpi<N>(1);
        Z = mpi<N>(1);
        T = mpi<N>(0);
    }

    ecp_ext& operator=(ecp_ext<N> const& other)
    {
        X = other.X;
        Y = other.Y;
        Z = other.Z;
        T = other.T;
        return *this;
    }

    mpi<N> X;
    mpi<N> Y;
    mpi<N> Z;
    mpi<N> T;
};

// XZ coordinates
template <unsigned N>
struct ecp_xz {
    ecp_xz() : X(1), Z(0){};
    ecp_xz(mpi<N> const& x, mpi<N> const& z) : X(x), Z(z){};
    ecp_xz(mpi<N> const& x) : X(x), Z(mpi<N>(1)){};
    ecp_xz(ecp_xz<N> const& point) : X(point.X), Z(point.Z){};

    mpi<N> X;
    mpi<N> Z;
};

template <unsigned N>
inline bool operator==(ecp_aff<N> const& p1, ecp_aff<N> const& p2)
{
    return p1.X.compare(p2.X) == 0 && p1.Y.compare(p2.Y) == 0;
}

// Jacobian coordinates
template <unsigned N>
struct ecp_jac {
    ecp_jac() : X(1), Y(1), Z(0){}; // ZERO by default
    ecp_jac(ecp_aff<N> const& point) : X(point.X), Y(point.Y), Z(1){};
    ecp_jac(mpi<N> const& x, mpi<N> const& y, mpi<N> const& z) : X(x), Y(y), Z(z){};
    ecp_jac(ecp_jac<N> const& point) : X(point.X), Y(point.Y), Z(point.Z){};

    ecp_jac<N>& operator=(ecp_jac<N> const& other)
    {
        X = other.X;
        Y = other.Y;
        Z = other.Z;
        return *this;
    }

    ecp_jac<N>& operator=(ecp_aff<N> const& other)
    {
        if (other.is_zero()) {
            *this = ecp_jac<N>{};
        } else {
            X = other.X;
            Y = other.Y;
            Z = mpi<N>(1);
        }
        return *this;
    }

    inline bool is_zero() const { return Z.equal(0); }

    mpi<N> X;
    mpi<N> Y;
    mpi<N> Z;
};

#ifdef ECC_STAT
#define INC_I()   \
    do {          \
        stat.I++; \
    } while (0)

#define INC_M()   \
    do {          \
        stat.M++; \
    } while (0)

#define INC_A()   \
    do {          \
        stat.A++; \
    } while (0)

#define INC_FA()   \
    do {           \
        stat.FA++; \
    } while (0)

#define INC_R()   \
    do {          \
        stat.R++; \
    } while (0)

#define INC_F()   \
    do {          \
        stat.F++; \
    } while (0)
#else
#define INC_I()
#define INC_M()
#define INC_A()
#define INC_FA()
#define INC_R()
#define INC_F()
#endif

struct ops_stat {
    int I{};
    int M{};
    int A{};
    int R{};
    int FA{};
    int F{};

    void clear()
    {
        I = M = A = FA = R = F = 0;
    }
};

struct curve_info {
    const char* name;
};

template <unsigned N, typename DP>
class curve_base
{
protected:
    mutable ops_stat stat;

public:
    ops_stat get_statis() const
    {
        return stat;
    }

protected:
    randomizer* rng;
    mpi<N> zero;

    curve_base() {}

    void normalize(ecp_proj<N> const& p1, ecp_aff<N>& p) const
    {
        mpi<N> z1;
        ffp_inverse(p1.Z, z1);
        ffp_mul(z1, p1.X, p.X);
        ffp_mul(z1, p1.Y, p.Y);
    }

    void normalize(ecp_ext<N> const& p1, ecp_aff<N>& p) const
    {
        mpi<N> z1;
        ffp_inverse(p1.Z, z1);
        ffp_mul(z1, p1.X, p.X);
        ffp_mul(z1, p1.Y, p.Y);
    }

    int randomize(ecp_jac<N>& P) const
    {
        mpi<N> r;
        mpi<N> one(1);
        int count = 0;

        do {
            r.randomize(rng);
            while (r.compare(DP::dp_p) >= 0) {
                r.rshift(1);
            }
            if (++count >= 20)
                return -1;
        } while (r.compare(one) <= 0);

        ffp_mul(P.Z, r, P.Z);

        mpi<N> rr;
        ffp_mul(r, r, rr);
        ffp_mul(P.X, rr, P.X);

        ffp_mul(rr, r, rr);
        ffp_mul(P.Y, rr, P.Y);
        return 0;
    }

    int randomize(ecp_xz<N>& P) const
    {
        mpi<N> r;
        mpi<N> one(1);
        int count = 0;

        do {
            r.randomize(rng);
            while (r.compare(DP::dp_p) >= 0) {
                r.rshift(1);
            }
            if (++count >= 20)
                return -1;
        } while (r.compare(one) <= 0);

        // X = r X, Z = r Z
        ffp_mul(P.X, r, P.X);
        ffp_mul(P.Z, r, P.Z);
        return 0;
    }

    int randomize(ecp_proj<N>& P) const
    {
        mpi<N> r;
        mpi<N> one(1);
        int count = 0;

        do {
            r.randomize(rng);
            while (r.compare(DP::dp_p) >= 0) {
                r.rshift(1);
            }
            if (++count >= 20)
                return -1;
        } while (r.compare(one) <= 0);

        // X = r X, Y = r Y, Z = r Z
        ffp_mul(P.X, r, P.X);
        ffp_mul(P.Y, r, P.Y);
        ffp_mul(P.Z, r, P.Z);
        return 0;
    }

    int randomize(ecp_ext<N>& P) const
    {
        mpi<N> r;
        mpi<N> one(1);
        int count = 0;

        do {
            r.randomize(rng);
            while (r.compare(DP::dp_p) >= 0) {
                r.rshift(1);
            }
            if (++count >= 20)
                return -1;
        } while (r.compare(one) <= 0);

        ffp_mul(P.Z, r, P.Z);
        ffp_mul(P.X, r, P.X);
        ffp_mul(P.Y, r, P.Y);
        ffp_mul(P.T, r, P.T);
        return 0;
    }

    // once single addition for reduction
    // to counter side channel attack
    template <unsigned M>
    inline void fast_reduce_once(mpi<M>& A) const
    {
        if (A.compare(DP::dp_p) >= 0) {
            A.sub(DP::dp_p);
        } else if (A.sign < 0 && !A.equal(0)) {
            A.add(DP::dp_p);
        } else {
            A.add(zero);
        }
    }

    // fast_reduce_once seemly constant execution for fast reduce
    template <unsigned M>
    inline void fast_reduce(mpi<M>& A, mpi<N>& r) const
    {
        fast_reduce_once(A);

        resize(A, r);
    }

    template <unsigned N1, unsigned N2>
    inline void ffp_mul(mpi<N1> const& A, mpi<N2> const& B, mpi<N>& r) const
    {
        mpi<N * 2> C;
        mul(A, B, C);
        ffp_reduce(C, r);
        INC_M();
    }

    template <unsigned K>
    inline void ffp_lshift(mpi<N> const& A, mpi<N>& r) const
    {
        mpi<N * 2> tmp(A);
        tmp.lshift(K);
        ffp_reduce(tmp, r);
        INC_F();
    }

    // ffp_reduce replies on the quality of fast modulus for specific curve
    // this can cause trouble for side-channel attack when we have to do
    // a number of input-dependent sub/add to further reduce result
    // as a consequence we try to compensate for sca by adding dummy addition
    inline void ffp_reduce(mpi<N * 2> const& A, mpi<N>& r) const
    {
        auto x = DP::reduce(A);
        INC_R();

        if constexpr (DP::sca) {
            while (x.compare(DP::dp_p) >= 0) {
                INC_FA();
                x.sub(DP::dp_p);
            }

            while (x.sign < 0 && !x.equal(0)) {
                INC_FA();
                x.add(DP::dp_p);
            }

#ifdef ECC_STAT
            if (stat.FA < stat.R) {
                INC_FA();
                x.add(zero);
            }
#endif
        }

        resize(x, r);
    }

    // both A/B are smaller than p, fast reduction is possible
    inline void ffp_add(mpi<N> const& A, mpi<N> const& B, mpi<N>& r) const
    {
        mpi<N + 1> tmp(A);
        tmp.add(B);
        fast_reduce(tmp, r);

        INC_A();
    }

    // both A/B are smaller than p, fast reduction is possible
    inline void ffp_sub(mpi<N> const& A, mpi<N> const& B, mpi<N>& r) const
    {
        mpi<N + 1> tmp(A);
        tmp.sub(B);
        fast_reduce(tmp, r);

        INC_A();
    }

    inline void ffp_inverse(mpi<N> const& A, mpi<N>& r) const
    {
        DP::iop(A, r);
        INC_I();
    }
};

// primary template for common edward curve
// It is not used in current code setup
// but it can be a good start point to support other parameters
template <CurveId CID, typename DP, unsigned N>
struct ed_point_op : public curve_base<N, DP> {
    using ecbase = curve_base<N, DP>;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_mul;
    using ecbase::ffp_reduce;
    using ecbase::ffp_sub;
    using ecp_target = ecp_proj<N>;
    static inline constexpr unsigned int encoded_size = mpi<N>::num_of_bytes;
    using prvkey_t = mpi<N>;

    ed_point_op() : ecbase(){};

protected:
    bool genprivkey(prvkey_t& prv) const
    {
        mpi<N> d;
        int count = 0;

        do {
            d.randomize(ecbase::rng);
            d.rshift(N - DP::dp_n.bitlen());
            if (++count > 20) {
                return false;
            }
        } while (d.equal(0) || d.compare(DP::dp_n) >= 0);

        prv = d;
        return true;
    }

    // p3 = p1 + p2
    void point_addition(ecp_proj<N> const& p1, ecp_proj<N> const& p2, ecp_proj<N>& p3) const
    {
        mpi<N> A;
        ffp_mul(p1.Z, p2.Z, A);

        mpi<N> B;
        ffp_mul(A, A, B);

        mpi<N> C;
        ffp_mul(p1.X, p2.X, C);

        mpi<N> D;
        ffp_mul(p1.Y, p2.Y, D);

        mpi<N> E;
        ffp_mul(C, D, E);
        ffp_mul(DP::dp_d, E, E);

        mpi<N> F;
        ffp_sub(B, E, F);

        mpi<N> G;
        ffp_add(B, E, G);

        mpi<N> t1, t2;
        ffp_add(p1.X, p1.Y, t1);
        ffp_add(p2.X, p2.Y, t2);
        ffp_mul(t1, t2, t1);
        ffp_sub(t1, C, t1);
        ffp_sub(t1, D, t1);
        ffp_mul(t1, F, t1);
        ffp_mul(t1, A, p3.X);

        mpi<N> Y3;
        if constexpr (DP::special_a == -1) {
            ffp_sub(DP::dp_p, C, C);
        } else if constexpr (DP::special_a != 1) {
            ffp_mul(C, DP::a, C);
        }
        ffp_sub(D, C, Y3);
        ffp_mul(G, Y3, Y3);
        ffp_mul(A, Y3, p3.Y);

        ffp_mul(F, G, p3.Z);
    }

    // P3 = 2P1
    void point_double(ecp_proj<N> const& p1, ecp_proj<N>& p3) const
    {
        mpi<N> B;
        ffp_add(p1.X, p1.Y, B);
        ffp_mul(B, B, B);

        mpi<N> C;
        ffp_mul(p1.X, p1.X, C);

        mpi<N> D;
        ffp_mul(p1.Y, p1.Y, D);

        mpi<N> E;
        // E = a * C
        if constexpr (DP::special_a == -1) {
            ffp_sub(DP::dp_p, C, E);
        } else if constexpr (DP::special_a == 1) {
            E = C;
        } else {
            ffp_mul(DP::dp_a, C, E);
        }

        mpi<N> F;
        ffp_add(E, D, F);

        mpi<N> H;
        ffp_mul(p1.Z, p1.Z, H);

        mpi<N> J;
        ffp_add(H, H, J);
        ffp_sub(F, J, J);

        mpi<N> X3;
        ffp_sub(B, C, X3);
        ffp_sub(X3, D, X3);
        ffp_mul(X3, J, p3.X);

        mpi<N> Y3;
        ffp_sub(E, D, Y3);
        ffp_mul(F, Y3, p3.Y);

        ffp_mul(F, J, p3.Z);
    }

    void point_mul(mpi<N> const& k, ecp_aff<N> const& P, ecp_aff<N>& r) const
    {
        ecp_proj<N> R;
        point_mul(k, P, R);
        normalize(R, r);
    }

    void point_mul(mpi<N> const& k, ecp_aff<N> const& P, ecp_target& r) const
    {
        ecbase::stat.clear();
        auto Qr = make_target(P);
        ecbase::randomize(Qr);

        // dummy operation for zero bits
        ecp_target Q = Qr;
        bit_iterator it(k);
        while (!it.end() && !it) {
            point_double(Q, Q);
            point_addition(Q, Qr, Q);
            it.next();
        }
        if (!it.end())
            it.next();

        Q = Qr;
        ecp_target Q1;

        while (!it.end()) {
            point_double(Q, Q);
            if (it) {
                point_addition(Q, Qr, Q);
            } else {
                point_addition(Q, Qr, Q1);
            }

            it.next();
        }

        r = Q;
    }

    inline ecp_target make_target(ecp_aff<N> const& aff) const
    {
        return ecp_target(aff);
    }
};

// specialization for ed448 operation
template <typename DP, unsigned N>
class ed_point_op<CurveId::ED448, DP, N> : public curve_base<N, DP>
{
public:
    using ecbase = curve_base<N, DP>;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_mul;
    using ecbase::ffp_reduce;
    using ecbase::ffp_sub;
    using ecp_target = ecp_proj<N>;
    using ecbase::normalize;

    // E448 key is a weirdo, longer than prime p
    static inline constexpr unsigned int encoded_size = 57;
    using prvkey_t = mpi<encoded_size * 8>;

    ed_point_op() : ecbase(){};

    bool genprivkey(prvkey_t& prv) const
    {
        prvkey_t d;
        int count = 0;

        do {
            d.randomize(ecbase::rng);
            if (++count > 20) {
                return false;
            }
        } while (d.equal(0));

        prv = d;
        return true;
    }

    // refer to RFC8032
    bool recover_x(mpi<N> const& Y, int x0, mpi<N>& X) const
    {
        mpi<N> u, v;
        mpi<N> one(1);
        ffp_mul(Y, Y, u);
        ffp_mul(u, DP::dp_d, v);
        ffp_sub(u, one, u);
        ffp_sub(v, DP::dp_a, v);

        // x is the square root of u/v
        // calcuate the candidate root x = (u/v)^((p + 1)/4)
        // which is : u^3 v (u^5 v^3)^((p - 3)/4)
        mpi<N> v3;
        ffp_mul(v, v, v3);
        ffp_mul(v3, v, v3);
        mpi<N> u3, u5;
        {
            mpi<N> u2;
            ffp_mul(u, u, u2);
            ffp_mul(u2, u, u3);
            ffp_mul(u2, u3, u5);
        }

        mpi<N> g;
        ffp_mul(u5, v3, g);
        mpi<N> e = DP::dp_p;
        e.sub(mpi<1>(3));
        e.rshift(2);
        mpi<N> x;
        fp_exp_mod(g, e, DP::dp_p, x);
        ffp_mul(x, v, x);
        ffp_mul(x, u3, x);

        mpi<N> vx2;
        ffp_mul(x, x, vx2);
        ffp_mul(v, vx2, vx2);
        if (vx2.compare(u) != 0) {
            return false;
        }

        if (x0 && x.equal(0)) {
            return false;
        }

        if (x0 != (x[0] & 1)) {
            ffp_sub(DP::dp_p, x, x);
        }

        X = x;
        return true;
    }

protected:
    struct internal_hasher {
        static inline constexpr int output_size = 114;

        default_shake256 sha;

        internal_hasher() : sha(output_size)
        {
            init();
        }

        // ED25519 context should always be empty
        void init(const uint8_t* ctx = NULL, size_t ctxlen = 0)
        {
            sha.init();
            if (ctxlen != 0) {
                static const char prefix[] = "SigEd448";
                sha.update(reinterpret_cast<const uint8_t*>(prefix), 8);
                uint8_t hflag = DP::phflag;
                sha.update(&hflag, 1);
                uint8_t ctxsize = static_cast<uint8_t>(ctxlen);
                sha.update(&ctxsize, 1);
                sha.update(ctx, ctxlen);
            }
        }

        void update(const uint8_t* data, size_t len)
        {
            sha.update(data, len);
        }

        auto finish()
        {
            FixedSizeBuf<output_size> buf;
            finish(buf);
            return buf;
        }

        void finish(FixedSizeBuf<output_size>& buf)
        {
            sha.finish(&buf[0]);
        }
    };

    ecp_target make_target(ecp_aff<N> const& aff) const
    {
        return ecp_proj<N>(aff);
    }

    void point_double(ecp_proj<N> const& p1, ecp_proj<N>& p3) const
    {
        mpi<N> B;
        ffp_add(p1.X, p1.Y, B);
        ffp_mul(B, B, B);

        mpi<N> C;
        ffp_mul(p1.X, p1.X, C);

        mpi<N> D;
        ffp_mul(p1.Y, p1.Y, D);

        mpi<N> E;
        ffp_add(C, D, E);

        mpi<N> H;
        ffp_mul(p1.Z, p1.Z, H);

        mpi<N> J;
        ffp_add(H, H, J);
        ffp_sub(E, J, J);

        mpi<N> X3;
        ffp_sub(B, E, X3);
        ffp_mul(X3, J, p3.X);

        mpi<N> Y3;
        ffp_sub(C, D, Y3);
        ffp_mul(E, Y3, p3.Y);

        ffp_mul(E, J, p3.Z);
    }

    void point_addition(ecp_proj<N> const& p1, ecp_proj<N> const& p2, ecp_proj<N>& p3) const
    {
        mpi<N> A;
        ffp_mul(p1.Z, p2.Z, A);

        mpi<N> B;
        ffp_mul(A, A, B);

        mpi<N> C;
        ffp_mul(p1.X, p2.X, C);

        mpi<N> D;
        ffp_mul(p1.Y, p2.Y, D);

        mpi<N> E;
        ffp_mul(DP::dp_d, C, E);
        ffp_mul(E, D, E);

        mpi<N> F;
        ffp_sub(B, E, F);

        mpi<N> G;
        ffp_add(B, E, G);

        mpi<N> t;
        ffp_add(p1.X, p1.Y, t);
        mpi<N> H;
        ffp_add(p2.X, p2.Y, H);
        ffp_mul(t, H, H);

        mpi<N> X3;
        ffp_sub(H, C, X3);
        ffp_sub(X3, D, X3);
        ffp_mul(F, X3, X3);
        ffp_mul(A, X3, p3.X);

        mpi<N> Y3;
        ffp_sub(D, C, Y3);
        ffp_mul(G, Y3, Y3);
        ffp_mul(A, Y3, p3.Y);

        ffp_mul(F, G, p3.Z);
    }
};

// specialization for ed25519 operation
template <typename DP, unsigned N>
class ed_point_op<CurveId::ED25519, DP, N> : public curve_base<N, DP>
{
public:
    using ecbase = curve_base<N, DP>;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_mul;
    using ecbase::ffp_reduce;
    using ecbase::ffp_sub;
    using ecp_target = ecp_ext<N>;
    using ecbase::normalize;
    static inline constexpr unsigned int encoded_size = mpi<N>::num_of_bytes;
    using prvkey_t = mpi<N>;

    ed_point_op() : ecbase(){};

    bool genprivkey(prvkey_t& prv) const
    {
        mpi<N> d;
        int count = 0;

        do {
            d.randomize(ecbase::rng);
            d.rshift(N - DP::dp_n.bitlen());
            if (++count > 20) {
                return false;
            }
        } while (d.equal(0) || d.compare(DP::dp_n) >= 0);

        prv = d;
        return true;
    }

    // refer to RFC8032
    bool recover_x(mpi<N> const& Y, int x0, mpi<N>& X) const
    {
        mpi<N> u, v;
        mpi<N> one(1);
        ffp_mul(Y, Y, u);
        ffp_mul(u, DP::dp_d, v);
        ffp_sub(u, one, u);
        ffp_sub(v, DP::dp_a, v);

        // x is the square root of u/v, here we go
        // calcuate the candidate root x = (u/v)^((p + 3)/8)
        // which is : u v^3 (u v^7)^((p - 5)/8)
        mpi<N> v3;
        ffp_mul(v, v, v3);
        ffp_mul(v3, v, v3);
        mpi<N> v7;
        ffp_mul(v3, v3, v7);
        ffp_mul(v7, v, v7);

        mpi<N> g;
        ffp_mul(u, v7, g);
        mpi<N> e = DP::dp_p;
        e.sub(mpi<1>(5));
        e.rshift(3);
        mpi<N> x;
        fp_exp_mod(g, e, DP::dp_p, x);
        ffp_mul(x, v3, x);
        ffp_mul(x, u, x);

        mpi<N> vx2;
        ffp_mul(x, x, vx2);
        ffp_mul(v, vx2, vx2);
        if (vx2.compare(u) == 0) {
        } else if (vx2.compare(DP::dp_p - u) == 0) {
            e = DP::dp_p;
            e.sub(mpi<1>(1));
            e.rshift(2);
            mpi<N> r;
            fp_exp_mod(mpi<N>(2), e, DP::dp_p, r);
            ffp_mul(x, r, x);
        } else {
            return false;
        }

        if (x0 && x.equal(0)) {
            return false;
        }

        if (x0 != (x[0] & 1)) {
            ffp_sub(DP::dp_p, x, x);
        }

        X = x;
        return true;
    }

protected:
    struct internal_hasher {
        default_sha512 sha;

        internal_hasher()
        {
            init();
        }

        // ED25519 context should always be empty
        void init(const uint8_t* ctx = NULL, size_t ctxlen = 0)
        {
            sha.init();
        }

        void update(const uint8_t* data, size_t len)
        {
            sha.update(data, len);
        }

        auto finish()
        {
            FixedSizeBuf<64> buf;
            sha.finish(&buf[0]);
            return buf;
        }

        void finish(FixedSizeBuf<64>& buf)
        {
            sha.finish(&buf[0]);
        }
    };

    ecp_target make_target(ecp_aff<N> const& aff) const
    {
        ecp_ext<N> ext;
        ext.X = aff.X;
        ext.Y = aff.Y;
        ext.Z = mpi<N>(1);
        ffp_mul(aff.X, aff.Y, ext.T);
        return ext;
    }

    void point_double(ecp_ext<N> const& p1, ecp_ext<N>& p3) const
    {
        mpi<N> A;
        ffp_mul(p1.X, p1.X, A);

        mpi<N> B;
        ffp_mul(p1.Y, p1.Y, B);

        mpi<N> C;
        ffp_mul(p1.Z, p1.Z, C);
        ffp_add(C, C, C);

        mpi<N> H;
        ffp_add(A, B, H);

        mpi<N> E;
        ffp_add(p1.X, p1.Y, E);
        ffp_mul(E, E, E);
        ffp_sub(H, E, E);

        mpi<N> G;
        ffp_sub(A, B, G);

        mpi<N> F;
        ffp_add(C, G, F);

        ffp_mul(E, F, p3.X);
        ffp_mul(G, H, p3.Y);
        ffp_mul(E, H, p3.T);
        ffp_mul(F, G, p3.Z);
    }

    void point_addition(ecp_ext<N> const& p1, ecp_ext<N> const& p2, ecp_ext<N>& p3) const
    {
        mpi<N> A;
        mpi<N> T;
        ffp_sub(p1.Y, p1.X, T);
        ffp_sub(p2.Y, p2.X, A);
        ffp_mul(T, A, A);

        mpi<N> B;
        ffp_add(p1.Y, p1.X, T);
        ffp_add(p2.Y, p2.X, B);
        ffp_mul(T, B, B);

        mpi<N> C;
        ffp_add(p1.T, p1.T, C);
        ffp_mul(C, DP::dp_d, C);
        ffp_mul(C, p2.T, C);

        mpi<N> D;
        ffp_add(p1.Z, p1.Z, D);
        ffp_mul(D, p2.Z, D);

        mpi<N> E;
        ffp_sub(B, A, E);

        mpi<N> F;
        ffp_sub(D, C, F);

        mpi<N> G;
        ffp_add(D, C, G);

        mpi<N> H;
        ffp_add(B, A, H);

        ffp_mul(E, F, p3.X);
        ffp_mul(G, H, p3.Y);
        ffp_mul(E, H, p3.T);
        ffp_mul(F, G, p3.Z);
    }
};

// specialization for Edward curve
template <typename DP, unsigned N>
class ec_curve<CurveType::EDWARD, DP, N> : public ed_point_op<DP::curve_id, DP, N>
{
public:
    using edop = ed_point_op<DP::curve_id, DP, N>;
    using ecp_target = typename edop::ecp_target;
    using ecbase = typename edop::ecbase;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_mul;
    using ecbase::ffp_reduce;
    using ecbase::ffp_sub;
    using ecbase::normalize;
    using edop::genprivkey;
    using edop::point_addition;
    using edop::point_double;
    using edop::recover_x;

    using encoded_point = FixedSizeBuf<edop::encoded_size>;
    using eddsa_sig = FixedSizeBuf<edop::encoded_size * 2>;
    using signature_t = eddsa_sig;

    using prvkey_t = typename edop::prvkey_t;
    using pubkey_t = encoded_point;

    struct keypair {
        prvkey_t d;
        pubkey_t Q;
    };

    ec_curve() : edop(){};

    ecp_aff<N> getG()
    {
        return ecp_aff<N>(DP::dp_gx, DP::dp_gy);
    }

    template <unsigned M>
    bool scalar_mul(mpi<M> const& d, pubkey_t const& pk, ecp_aff<N>& r) const
    {
        pubkey_t pub;
        prvkey_t s;
        genpk(d, s, pub);

        ecp_aff<N> A;
        if (!decoding(pk, A)) {
            return false;
        }

        point_mul(s, A, r);
        return true;
    }

    void genpk(prvkey_t const& d, prvkey_t& s, pubkey_t& pk) const
    {
        typename edop::internal_hasher sha{};
        mpi_hash(d, sha);
        auto hash = sha.finish();

        // hash treated used as little endian format
        // while our mpi always import in big-endian
        shuffle(&hash[0], s.num_of_bytes);

        if constexpr (DP::curve_id == CurveId::ED448) {
            hash[0] = 0;
            hash[1] |= 0x80;
            hash[s.num_of_bytes - 1] &= 0xFC;
        } else if constexpr (DP::curve_id == CurveId::ED25519) {
            hash[s.num_of_bytes - 1] &= 0xF8;
            hash[0] &= ~0x80;
            hash[0] |= 0x40;
        }

        s = prvkey_t(&hash[0], s.num_of_bytes);

        ecp_aff<N> Q;
        point_mul(s, Q);
        pk = encoding(Q);
    }

    auto genpk(prvkey_t const& d) const
    {
        prvkey_t s;
        pubkey_t pk;
        genpk(d, s, pk);
        return pk;
    }

    bool genkeypair(keypair& kp) const
    {
        prvkey_t d;
        auto ret = genprivkey(d);

        if (ret) {
            prvkey_t s;
            genpk(d, s, kp.Q);
            kp.d = d;
        }

        return ret;
    }

    bool eddsa_sign(prvkey_t const& d, const uint8_t* m, size_t mlen, signature_t& sig, const uint8_t* ctx = nullptr, size_t ctxlen = 0)
    {
        typename edop::internal_hasher sha{};
        mpi_hash(d, sha);
        auto hash = sha.finish();

        // hash treated used as little endian format
        // while our mpi always import in big-endian
        shuffle(&hash[0], d.num_of_bytes);

        if constexpr (DP::curve_id == CurveId::ED448) {
            hash[0] = 0;
            hash[1] |= 0x80;
            hash[prvkey_t::num_of_bytes - 1] &= 0xFC;
        } else if constexpr (DP::curve_id == CurveId::ED25519) {
            hash[prvkey_t::num_of_bytes - 1] &= 0xF8;
            hash[0] &= ~0x80;
            hash[0] |= 0x40;
        }
        prvkey_t s(&hash[0], prvkey_t::num_of_bytes);

        ecp_aff<N> Q;
        point_mul(s, Q);
        auto A = encoding(Q); // pk

        // HASH(dom(F, C) || prefix || PH(M)
        sha.init(ctx, ctxlen);
        sha.update(&hash[hash.buf_size / 2], hash.buf_size / 2);
        sha.update(m, mlen);
        sha.finish(hash);

        // R = r B
        mpi<N> r;
        standard_modulus_op mop(DP::dp_n);
        shuffle(&hash[0], hash.buf_size);
        {
            mpi<hash.buf_size * 8> tmp(&hash[0], hash.buf_size);
            mop(tmp, r);
        }

        encoded_point R;
        {
            ecp_aff<N> _R;
            point_mul(r, _R);
            R = encoding(_R);
        }

        // HASH(dom(F, C) || R || A || PH(M))
        sha.init(ctx, ctxlen);
        buf_hash(R, sha);
        buf_hash(A, sha);
        sha.update(m, mlen);
        sha.finish(hash);

        mpi<N> k;
        shuffle(&hash[0], hash.buf_size);
        {
            mpi<hash.buf_size * 8> tmp(&hash[0], hash.buf_size);
            mop(tmp, k);
        }

        // S = (r + k * s) mod L
        FixedSizeBuf<edop::encoded_size> S;
        {
            mpi<N> _S;
            mop(k * s, _S);
            mpi<N + 1> tmp(_S);
            tmp.add(r);
            mop(tmp, _S);

            // S is one byte longer than _S for ed448
            S.copy(_S.output(false));
            if constexpr (DP::curve_id == CurveId::ED448) {
                S[S.buf_size - 2] &= 0x3F;
            } else if constexpr (DP::curve_id == CurveId::ED25519) {
                S[S.buf_size - 1] &= 0x1F;
            }
        }

        // concatenate R & S
        sig.copy(R);
        sig.copy(S, R.size());
        return true;
    }

    bool eddsa_verify(encoded_point const& pk, const uint8_t* m, size_t mlen, signature_t const& sig, const uint8_t* ctx = nullptr, size_t ctxlen = 0)
    {
        // first half R, second half S
        encoded_point r;

        // NOTE for ed448 S is actually mpi<N>::num_of_bytes + 1
        // however, the last bytes is always zero
        FixedSizeBuf<mpi<N>::num_of_bytes> s;
        s.copy(&sig[r.size()], s.size());
        shuffle(&s[0], s.size());
        mpi<N> S(&s[0], s.size());
        if (S.compare(DP::dp_n) >= 0) {
            return false;
        }

        ecp_aff<N> A;
        if (!decoding(pk, A)) {
            return false;
        }

        ecp_aff<N> R;
        r.copy(&sig[0], r.size());
        if (!decoding(r, R)) {
            return false;
        }

        // compuate SHA512(dom2(F, C) || R || A || PH(M)
        typename edop::internal_hasher sha{};
        sha.init(ctx, ctxlen);
        buf_hash(r, sha);
        buf_hash(pk, sha);
        sha.update(m, mlen);
        auto hash = sha.finish();

        shuffle(&hash[0], hash.buf_size);
        mpi<N> k;
        {
            mpi<hash.buf_size * 8> k2(&hash[0], hash.buf_size);
            fp_reduce(k2, DP::dp_n, k);
        }

        // check [S]B = R + [k]A
        ecp_target R1, R2;
        point_mul(S, R1);
        point_mul(k, A, R2);
        edop::point_addition(edop::make_target(R), R2, R2);
        return equal(R1, R2);
    }

    // a * x^2 + y^2 = 1 + d * x^2 y^2
    bool verify_point(ecp_aff<N> const& Q) const
    {
        mpi<N> XX;
        ffp_mul(Q.X, Q.X, XX);
        mpi<N> YY;
        ffp_mul(Q.Y, Q.Y, YY);
        mpi<N> lhs;

        // lhs = a * XX + YY
        if constexpr (DP::special_a == -1) {
            ffp_sub(DP::dp_p, XX, lhs);
        } else if constexpr (DP::special_a == 1) {
            lhs = XX;
        } else {
            ffp_mul(DP::dp_a, XX, lhs);
        }
        ffp_add(lhs, YY, lhs);

        mpi<N> rhs;
        ffp_mul(DP::dp_d, XX, rhs);
        ffp_mul(rhs, YY, rhs);
        ffp_add(mpi<N>(1), rhs, rhs);

        return lhs.compare(rhs) == 0;
    }

    bool verify_pk(encoded_point const& p) const
    {
        ecp_aff<N> R;
        return decoding(p, R);
    }

    // conversion:
    // (u, v) = ((1+y)/(1-y), sqrt(-486664)*u/x)
    ecp_aff<N> to_montgomery(ecp_aff<N> const& ed)
    {
        mpi<N> one(1);
        mpi<N> t;
        ffp_sub(one, ed.Y, t);

        ffp_inverse(t, t);
        mpi<N> u;
        ffp_add(one, ed.Y, u);
        ffp_mul(t, u, u);

        mpi<N> v;
        ffp_inverse(ed.X, v);
        ffp_mul(u, v, v);
        ffp_mul(v, DP::sqrt486664, v);

        return ecp_aff<N>(u, v);
    }

    mpi<N> get_sqrt(mpi<N> const& u, mpi<N> const& v, mpi<N> const& x)
    {
        mpi<N> t;
        ffp_inverse(u, t);
        ffp_mul(t, x, t);
        ffp_mul(t, v, t);
        return t;
    }

    bool decoding(encoded_point const& p, ecp_aff<N>& R) const
    {
        encoded_point pp = p;
        shuffle(&pp[0], pp.buf_size);
        int x0 = !!(pp[0] & 0x80);
        pp[0] &= 0x7F;

        // ED448 encodded size is longer than Y coordinates
        // so the start offset to copy Y might not be zero
        mpi<N> Y(&pp[pp.buf_size - mpi<N>::num_of_bytes], mpi<N>::num_of_bytes);
        mpi<N> X;

        if (Y.compare(DP::dp_p) >= 0) {
            return false;
        }

        auto ret = recover_x(Y, x0, X);
        if (ret) {
            R = ecp_aff<N>(X, Y);
        }

        return ret;
    }

private:
    template <unsigned M>
    inline void point_mul(mpi<M> const& k, ecp_target& r) const
    {
        point_mul(k, ecp_aff<N>(DP::dp_gx, DP::dp_gy), r);
    }

    template <unsigned M>
    inline void point_mul(mpi<M> const& k, ecp_aff<N>& r) const
    {
        point_mul(k, ecp_aff<N>(DP::dp_gx, DP::dp_gy), r);
    }

    template <unsigned M>
    void point_mul(mpi<M> const& k, ecp_aff<N> const& P, ecp_aff<N>& r) const
    {
        ecp_target R;
        point_mul(k, P, R);
        normalize(R, r);
    }

    template <unsigned M>
    void point_mul(mpi<M> const& k, ecp_aff<N> const& P, ecp_target& r) const
    {
        ecbase::stat.clear();
        auto Qr = edop::make_target(P);
        ecbase::randomize(Qr);

        // dummy operation for zero bits
        ecp_target Q = Qr;
        bit_iterator it(k);
        while (!it.end() && !it) {
            point_double(Q, Q);
            point_addition(Q, Qr, Q);
            it.next();
        }
        if (!it.end())
            it.next();

        Q = Qr;
        ecp_target Q1;

        while (!it.end()) {
            point_double(Q, Q);
            if (it) {
                point_addition(Q, Qr, Q);
            } else {
                point_addition(Q, Qr, Q1);
            }

            it.next();
        }

        r = Q;
    }

    void shuffle(uint8_t* p, size_t len) const
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

    auto encoding(ecp_aff<N> const& Q) const
    {
        mpi<edop::encoded_size * 8> Y(Q.Y);
        Y.set_bit(Y.num_of_bits - 1, Q.X.get_bit(0));
        return Y.output(false);
    }

    template <unsigned M, typename H>
    inline void mpi_hash(mpi<M> const& m, H& op) const
    {
        auto o = m.output(false);
        op.update(&o[0], o.size());
    }

    template <unsigned SIZE, typename H>
    inline void buf_hash(FixedSizeBuf<SIZE> const& buf, H& op)
    {
        op.update(&buf[0], buf.size());
    }

    // this works for projective coordinates and
    // extended homogenious coordiates
    bool equal(ecp_target const& p1, ecp_target const& p2)
    {
        mpi<N> t1, t2;
        ffp_mul(p1.X, p2.Z, t1);
        ffp_mul(p2.X, p1.Z, t2);
        if (t1.compare(t2) != 0)
            return false;

        ffp_mul(p1.Y, p2.Z, t1);
        ffp_mul(p2.Y, p1.Z, t2);
        if (t1.compare(t2) != 0)
            return false;

        return true;
    }
};

// specialization for montgomery curve
template <typename DP, unsigned N>
class ec_curve<CurveType::MONTGOMERY, DP, N> : public curve_base<N, DP>
{
public:
    using ecbase = curve_base<N, DP>;
    using ecbase::fast_reduce;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_mul;
    using ecbase::ffp_sub;
    using prvkey_t = mpi<N>;
    using pubkey_t = mpi<N>;
    struct keypair {
        prvkey_t d;
        pubkey_t Q;
    };
    using signature_t = void;

    ec_curve() : ecbase(){};

    const mpi<N>& getG()
    {
        return DP::dp_gx;
    }

    template <unsigned M>
    bool scalar_mul(mpi<M> const& d, mpi<N> const& X, ecp_aff<N>& r) const
    {
        ecp_xz<N> P(X, mpi<N>(1));
        ecp_xz<N> xz;
        point_mul(d, P, xz);
        normalize(xz, r.X);
        return true;
    }

    void clampkey(prvkey_t& d) const
    {
        static constexpr int nbits = DP::nbits;
        // private key clamping
        d.set_bit(nbits, 1);
        d.set_bit(0, 0);
        d.set_bit(1, 0);

        if (nbits == 254)
            d.set_bit(2, 0);
    }

    bool genprivkey(prvkey_t& prv) const
    {
        static constexpr int nbits = DP::nbits;

        mpi<N> d;
        int count = 0;

        do {
            d.randomize(ecbase::rng);
            if (++count >= 20)
                return false;
        } while (d.bitlen() == 0);

        // private key clamping
        auto b = d.bitlen() - 1;
        if (b > nbits) {
            d.rshift(b - nbits);
        } else {
            d.set_bit(nbits, 1);
        }

        // set last two bits for both 25519/448
        d.set_bit(0, 0);
        d.set_bit(1, 0);

        // do this only for 25519
        if (nbits == 254)
            d.set_bit(2, 0);

        prv = d;

        return true;
    }

    bool genkeypair(keypair& key) const
    {
        mpi<N> d;
        auto ret = genprivkey(d);

        if (ret) {
            point_mul(d, key.Q);
            key.d = d;
        }

        return ret;
    }

    auto genpk(prvkey_t const& d) const
    {
        pubkey_t pub;
        point_mul(d, pub);
        return pub;
    }

    bool verify_point(ecp_aff<N>& Q) const
    {
        if (!Q.X.equal(0) && !Q.Y.equal(0)) {
            mpi<N> lhs;
            ffp_mul(Q.Y, Q.Y, lhs);

            mpi<N> XX;
            ffp_mul(Q.X, Q.X, XX);
            mpi<N> XXX;
            ffp_mul(XX, Q.X, XXX);
            ffp_mul(XX, mpi<256>(DP::_A), XX);
            mpi<N> X(Q.X);
            fast_reduce(X, X);
            mpi<N> rhs;
            ffp_add(XXX, XX, rhs);
            ffp_add(rhs, X, rhs);
            return lhs.compare(rhs) == 0;
        }

        return true;
    }

    bool verify_pk(mpi<N> const& pk) const
    {
        return true;
    }

private:
    void point_double(ecp_xz<N> const& P, ecp_xz<N>& Double) const
    {
        // addition
        mpi<N> A;
        ffp_add(P.X, P.Z, A);
        mpi<N> B;
        ffp_sub(P.X, P.Z, B);

        mpi<N> AA, BB;
        ffp_mul(A, A, AA);
        ffp_mul(B, B, BB);
        mpi<N> T;
        ffp_sub(AA, BB, T);
        ffp_mul(AA, BB, Double.X);
        ffp_mul(T, DP::dp_a, Double.Z);
        ffp_add(BB, Double.Z, Double.Z);
        ffp_mul(T, Double.Z, Double.Z);
    }

    // Guide to Elliptic Curve Cryptography 13.2.3
    // to compute
    // Double = 2P
    // Sum = P + Q
    // where delta = X(P - Q)
    void point_double_addition(ecp_xz<N> const& P,
                               ecp_xz<N> const& Q,
                               ecp_xz<N> const& delta,
                               ecp_xz<N>& Double,
                               ecp_xz<N>& Sum) const
    {
        // addition
        mpi<N> A;
        ffp_add(P.X, P.Z, A);
        mpi<N> B;
        ffp_sub(P.X, P.Z, B);

        mpi<N> C;
        ffp_add(Q.X, Q.Z, C);
        mpi<N> D;
        ffp_sub(Q.X, Q.Z, D);

        mpi<N> AD;
        ffp_mul(A, D, AD);
        mpi<N> BC;
        ffp_mul(B, C, BC);

        mpi<N> X;
        ffp_add(AD, BC, X);
        ffp_mul(X, X, Sum.X);
        if (!delta.Z.equal(1)) {
            ffp_mul(delta.Z, Sum.X, Sum.X);
        }

        mpi<N> Z;
        ffp_sub(AD, BC, Z);
        ffp_mul(Z, Z, Z);
        ffp_mul(delta.X, Z, Sum.Z);

        // double
        mpi<N> AA, BB;
        ffp_mul(A, A, AA);
        ffp_mul(B, B, BB);
        mpi<N> T;
        ffp_sub(AA, BB, T);
        ffp_mul(AA, BB, Double.X);
        ffp_mul(T, DP::dp_a, Double.Z);
        ffp_add(BB, Double.Z, Double.Z);
        ffp_mul(T, Double.Z, Double.Z);
    }

    // montgomery ladder
    inline void point_mul(mpi<N> const& k, ecp_xz<N>& r) const
    {
        point_mul(k, ecp_xz<N>(DP::dp_gx, mpi<N>(1)), r);
    }

    inline void point_mul(mpi<N> const& k, mpi<N>& X) const
    {
        ecp_xz<N> r;
        point_mul(k, r);
        normalize(r, X);
    }

    void point_mul(mpi<N> const& k, ecp_xz<N> const& P, ecp_xz<N>& r) const
    {
        ecbase::stat.clear();
        ecp_xz<N> Pr(P);
        ecbase::randomize(Pr);

        ecp_xz<N> P1(Pr);
        ecp_xz<N> P2;

        point_double(P, P2);

        bit_iterator it(k);
        // dummy loop for zero bits
        while (!it.end() && !it) {
            point_double_addition(Pr, Pr, P, Pr, Pr);
            it.next();
        }

        if (!it.end())
            it.next();

        while (!it.end()) {
            if (!it) {
                point_double_addition(P1, P2, P, P1, P2);
            } else {
                point_double_addition(P2, P1, P, P2, P1);
            }
            it.next();
        }

        r = P1;
    }

    void normalize(ecp_xz<N> const& p1, mpi<N>& X) const
    {
        mpi<N> z1;
        ffp_inverse(p1.Z, z1);
        ffp_mul(z1, p1.X, X);
    }
};

// specialization for short weierstrass curve
template <typename DP, unsigned N>
class ec_curve<CurveType::SHORT_WEIERSTRASS, DP, N> : public curve_base<N, DP>
{
    static constexpr int curve_bits = N;
    static constexpr int WSIZE = 4;
    ecp_aff<N> presets[1 << WSIZE];

    // using a second table is more cost-effective than increasing the W size
    // double window size from 4 to 8 increase buffering from 1K to 16K!!!! (2X performance)
    // second table halve the double operation at cost of 2X memory (rougly +30% performance)
    ecp_aff<N> alt_presets[1 << WSIZE];

    static constexpr int comb_columns = (N + WSIZE - 1) / WSIZE;

public:
    using ecbase = curve_base<N, DP>;
    using ecbase::ffp_add;
    using ecbase::ffp_inverse;
    using ecbase::ffp_lshift;
    using ecbase::ffp_mul;
    using ecbase::ffp_sub;

    struct ecdsa_sig {
        ecdsa_sig() : r(0), s(0){};
        ecdsa_sig(mpi<N> const& r, mpi<N> const& s) : r(r), s(s){};
        ecdsa_sig(const uint8_t* rp, size_t rlen, const uint8_t* sp, size_t slen) : r(rp, rlen), s(sp, slen) {}
        mpi<N> r;
        mpi<N> s;
    };
    using signature_t = ecdsa_sig;
    using prvkey_t = mpi<N>;
    using pubkey_t = ecp_aff<N>;

    struct keypair {
        prvkey_t d;
        pubkey_t Q;
    };

    ec_curve() : ecbase()
    {
        precompute();
    }

    ecp_aff<N> getG()
    {
        return ecp_aff(DP::dp_gx, DP::dp_gy);
    }

    // scalar multiplication for unknown base
    bool scalar_mul(mpi<N> const& k, ecp_aff<N> const& P, ecp_aff<N>& r) const
    {
        if (!verify_pk(P))
            return false;

        ecp_jac<N> Q;
        ecp_jac<N> Qr(P);
        ecbase::randomize(Qr);

        bit_iterator it(k);

        // dummy operation for zero bits
        Q = Qr;
        while (!it.end() && !it) {
            point_double(Q, Q);
            point_addition(Q, P, Q);
            it.next();
        }

        if (!it.end())
            it.next();

        Q = Qr;
        ecp_jac<N> Q1;
        while (!it.end()) {
            point_double(Q, Q);
            if (it) {
                point_addition(Q, P, Q);
            } else {
                point_addition(Q, P, Q1);
            }
            it.next();
        }

        normalize(Q, r);
        return true;
    }

    void scalar_mul(mpi<N> const& k, ecp_aff<N>& r) const
    {
        point_mul(k, r);
    }

    template <unsigned M>
    bool ecdsa_sign(mpi<N> const& d, mpi<M> const& m, signature_t& sig)
    {
        mpi<N> r;
        mpi<N> s;

        // TODO perhaps we need only one helper to save startup time
        barrett_helper bh(DP::dp_n);
        montg_inverse_helper mih{DP::dp_n};
        keypair ephem; // ephemeral key

        do {
            do {
                if (!genkeypair(ephem))
                    return false;
                modn(ephem.Q.X, r);
            } while (r.equal(0));

            mpi<N> e;
            transform_msg(m, e);

            mpi<N> Rm;
            if (!genprivkey(Rm))
                return false;

            // s = (e + r * d) / k = Rm (e + r * d) / (k * Rm) mod n
            auto x = r * d;
            x.add(e);
            bh.reduce(x, s);
            bh.reduce(s * Rm, s);
            bh.reduce(ephem.d * Rm, Rm);
            mih(Rm, Rm);
            bh.reduce(Rm * s, s);
        } while (s.equal(0));

        sig.s = s;
        sig.r = r;
        return true;
    }

    template <unsigned M>
    bool ecdsa_sign(keypair const& key, mpi<M> const& m, signature_t& sig)
    {
        return ecdsa_sign(key.d, m, sig);
    }

    inline void modn(mpi<N> const& m, mpi<N>& r)
    {
        r = m;
        while (r.compare(DP::dp_n) >= 0) {
            r.sub(DP::dp_n);
        }
    }

    template <unsigned M>
    bool ecdsa_verify(ecp_aff<N> const& Q, mpi<M> const& m,
                      signature_t const& sig)
    {
        if (sig.r.compare(DP::dp_n) >= 0 ||
            sig.r.compare(mpi<N>(1)) < 0 ||
            sig.s.compare(DP::dp_n) >= 0 ||
            sig.s.compare(mpi<N>(1)) < 0) {
            return false;
        }

        mpi<N> e;
        transform_msg(m, e);

        barrett_helper bh(DP::dp_n);
        montg_inverse_helper mih{DP::dp_n};

        mpi<N> w;
        mih(sig.s, w);

        mpi<N> u1, u2;
        bh.reduce(e * w, u1);
        bh.reduce(sig.r * w, u2);
        ecp_jac<N> R;
        add_scalar_mul(u1, u2, Q, R);

        if (R.is_zero()) {
            return false;
        }

        mpi<N> x1;
        modn(R.X, x1);
        return x1.compare(sig.r) == 0;
    }

    bool genprivkey(mpi<N>& prv) const
    {
        mpi<N> d;
        int count = 0;

        do {
            d.randomize(ecbase::rng);
            d.rshift(N - DP::dp_n.bitlen());
            if (++count > 20) {
                return false;
            }
        } while (d.equal(0) || d.compare(DP::dp_n) >= 0);

        prv = d;
        return true;
    }

    int genkeypair(keypair& key) const
    {
        mpi<N> d;
        auto ret = genprivkey(d);

        if (ret) {
            point_mul(d, key.Q);
            key.d = d;
        }

        return ret;
    }

    auto genpk(prvkey_t const& d) const
    {
        pubkey_t pk;
        point_mul(d, pk);
        return pk;
    }

    bool verify_pk(ecp_aff<N> const& Q) const
    {
        auto& X = Q.X;
        auto& Y = Q.Y;

        if (X.compare(DP::dp_p) >= 0 || Y.compare(DP::dp_p) >= 0)
            return false;

        // verify Y^2 = X(X^2 + A) + B
        mpi<N> rhs;
        ffp_mul(X, X, rhs);

        if constexpr (DP::special_a == -3) {
            ffp_sub(rhs, mpi<N>(3), rhs);
        } else {
            ffp_add(rhs, DP::dp_a, rhs);
        }

        ffp_mul(rhs, X, rhs);
        ffp_add(rhs, DP::dp_b, rhs);

        mpi<N> YY;
        ffp_mul(Y, Y, YY);
        return (rhs.compare(YY) == 0);
    }

private:
    template <unsigned M>
    void transform_msg(mpi<M> const& m, mpi<N>& e)
    {
        resize<N>(m, e);
        auto bits1 = e.bitlen();
        auto bits2 = DP::dp_n.bitlen();
        if (bits1 > bits2) {
            e.rshift(bits1 - bits2);
        }

        while (e.compare(DP::dp_n) >= 0) {
            e.sub(DP::dp_n);
        }
    }

    // compute R = l P + k Q
    void add_scalar_mul(mpi<N> const& l, mpi<N> const& k,
                        ecp_aff<N> const& Q,
                        ecp_jac<N>& R) const
    {
        ecp_aff<N> P1, P2;
        point_mul(l, P1);

        point_mul(k, Q, P2);
        point_addition(ecp_jac<N>(P1), P2, R);
        if (!R.is_zero()) {
            normalize(R);
        }
    }

    void normalize(ecp_jac<N>& p) const
    {
        ecp_aff<N> p2;
        normalize(p, p2);
        p = ecp_jac<N>(p2);
    }

    void normalize(ecp_jac<N> const& p1, ecp_aff<N>& p2) const
    {
        if (p1.is_zero()) {
            p2.set_zero();
            return;
        }

        mpi<N> z1, z2;

        ffp_inverse(p1.Z, z1);

        ffp_mul(z1, z1, z2);

        // X = X / Z^2
        ffp_mul(p1.X, z2, p2.X);

        // Y = Y / Z^3
        ffp_mul(p1.Y, z2, p2.Y);
        ffp_mul(p2.Y, z1, p2.Y);
    }

    ecp_aff<N> normalize(ecp_jac<N> const& p1) const
    {
        ecp_aff<N> p2;
        normalize(p1, p2);
        return p2;
    }

    // save a bunch of inverse by normalizing together
    template <unsigned K>
    void normalize(std::array<ecp_jac<N>, K> const& p1, std::array<ecp_aff<N>, K>& p2) const
    {
        if constexpr (K == 1) {
            normalize(p1[0], p2[0]);
            return;
        }

        std::array<mpi<N>, K> c;
        c[0] = p1[0].Z;

        // Ci = Z0 * Z1 * ... * Zi
        for (unsigned i = 1; i < p1.size(); i++) {
            ffp_mul(c[i - 1], p1[i].Z, c[i]);
        }

        // u = 1 / (Z0 * ... * Zn) mod P
        mpi<N> u;
        DP::iop(c[K - 1], u);

        mpi<N> z1, z2;

        for (unsigned i = K; i;) {
            if (--i == 0) {
                z1 = u;
            } else {
                ffp_mul(c[i - 1], u, z1);
                ffp_mul(u, c[i], u);
            }

            ffp_mul(z1, z1, z2);

            // X = X / Z^2
            ffp_mul(p1[i].X, z2, p2[i].X);

            // Y = Y / Z^3
            ffp_mul(p1[i].Y, z2, p2[i].Y);
            ffp_mul(p2[i].Y, z1, p2[i].Y);
        }
    }

    // double a point
    // p3 = 2 * p1
    // Guide to Elliptic Curve Cryptography (3.13)
    // x3 = (3X^2 + a * Z^4)^2 - 8XY^2
    // y2 = (3X^2 + a * Z^4)(4XY^2 - x3) - 8Y^4
    // z3 = 2YZ
    //
    // in the special case a == -3
    // (3X^2 + a * Z^4) == 3(X + Z^2)(X - Z^2)
    void point_double(ecp_jac<N> const& p1, ecp_jac<N>& p3) const
    {
        // A = Y^2
        mpi<N> A;
        ffp_mul(p1.Y, p1.Y, A);

        // B = 4 X A
        mpi<N> B;
        ffp_mul(p1.X, A, B);
        //ffp_lshift<2>(B, B);
        ffp_mul(B, mpi<1>(4), B);

        // C = 8 A^2
        mpi<N> C;
        ffp_mul(A, A, C);
        //ffp_lshift<3>(C, C);
        ffp_mul(C, mpi<1>(8), C);

        mpi<N> D;
        if constexpr (DP::special_a == -3) {
            // D = 3(X + Z^2)(X - Z^2)
            mpi<N> z2;
            ffp_mul(p1.Z, p1.Z, z2);
            mpi<N> r1, r2;
            ffp_add(p1.X, z2, r1);
            ffp_sub(p1.X, z2, r2);
            ffp_mul(r1, r2, D);
            ffp_mul(D, mpi<1>(3), D);
        } else {
            // D = 3X^2 + a * Z^4
            ffp_mul(p1.X, p1.X, D);
            ffp_mul(D, mpi<1>(3), D);

            if constexpr (DP::special_a != 0) {
                mpi<N> z4;
                ffp_mul(p1.Z, p1.Z, z4);
                ffp_mul(z4, z4, z4);
                ffp_mul(DP::dp_a, z4, z4);
                ffp_add(D, z4, D);
            }
        }

        // x3 = D^2 - 2B
        mpi<N> X3;
        ffp_mul(D, D, X3);
        mpi<N> T;
        //ffp_lshift<1>(B, T);
        ffp_mul(B, mpi<1>(2), T);
        ffp_sub(X3, T, X3);

        // y3 = D (B - x3) - C
        mpi<N> Y3;
        ffp_sub(B, X3, T);
        ffp_mul(D, T, T);
        ffp_sub(T, C, Y3);

        // z3 = 2 Y Z
        mpi<N> Z3;
        ffp_mul(p1.Y, p1.Z, T);
        //ffp_lshift<1>(T, Z3);
        ffp_mul(T, mpi<1>(2), Z3);

        p3.X = X3;
        p3.Y = Y3;
        p3.Z = Z3;
    }

    // Guide to Elliptic Curve Cryptography (3.14)
    // p3 = p1 + p2 (result in jacobian coordinates)
    //
    // x3 = (y2 z1^3 - y1)^2 - (x2 z1^2 - x1)^2 (x1 + x2 z1^2)
    // y3 = (y2 z1^3 - y1)(x1(x2 z1^2 - x1)^2 - x3) - y1(x2 z1^2 - x1)^3
    // z3 = (x2 z1^2 - x1) z1
    void point_addition(ecp_jac<N> const& p1, ecp_aff<N> const& p2, ecp_jac<N>& p3) const
    {
        // TODO we could use less temp variables
        mpi<N> A;
        ffp_mul(p1.Z, p1.Z, A);

        mpi<N> B;
        ffp_mul(p1.Z, A, B);

        mpi<N> C;
        ffp_mul(p2.X, A, C);

        mpi<N> D;
        ffp_mul(p2.Y, B, D);

        mpi<N> E;
        ffp_sub(C, p1.X, E);

        mpi<N> F;
        ffp_sub(D, p1.Y, F);

        mpi<N> G;
        ffp_mul(E, E, G);

        mpi<N> H;
        ffp_mul(G, E, H);

        mpi<N> I;
        ffp_mul(p1.X, G, I);

        // x3 = F^2 - (H + 2 I)
        mpi<N> X3, T;
        ffp_mul(F, F, X3);
        //ffp_lshift<1>(I, T);
        ffp_mul(I, mpi<1>(2), T);
        ffp_add(H, T, T);
        ffp_sub(X3, T, X3);

        // y3 = F (I - x3) - y1 H
        mpi<N> Y3;
        ffp_sub(I, X3, T);
        ffp_mul(F, T, Y3);
        ffp_mul(p1.Y, H, T);
        ffp_sub(Y3, T, Y3);

        // z3 = z1 E
        mpi<N> Z3;
        ffp_mul(p1.Z, E, Z3);

        // this routine should always execute with constant time
        // regardless of whether adding zero
        if (p2.is_zero()) {
            p3 = p1;
        } else if (p1.is_zero()) {
            p3 = p2;
        } else {
            p3.X = X3;
            p3.Y = Y3;
            p3.Z = Z3;
        }
    }

    // for unknown base
    // r = kP
    // this routine should only be used for pre-computing and anything not private-key related
    // since it is subject to timing attack
    template <unsigned M>
    void point_mul(mpi<M> const& k, ecp_aff<N> const& P, ecp_aff<N>& r) const
    {
        ecbase::stat.clear();
        bit_iterator it(k);
        while (!it.end() && !it) {
            it.next();
        }

        if (it.end()) {
            r.set_zero();
            return;
        }

        ecp_jac<N> Q(P);
        it.next();

        while (!it.end()) {
            point_double(Q, Q);
            if (it) {
                point_addition(Q, P, Q);
            }

            it.next();
        }

        normalize(Q, r);
    }

    // r = kG
    // Guide to Ecliptic Crypgrography Alg. 3.45
    void point_mul(mpi<N> const& k, ecp_aff<N>& r) const
    {
        using comb_type = uint8_t;
        comb_type combs[comb_columns];

        static_assert(sizeof(comb_type) * 8 >= WSIZE, "combs has wrong unit size");

        ecbase::stat.clear();

        for (auto& x : combs) {
            x = 0;
        }

        int bit = 0;
        bit_iterator it(k, false); // iterate from low bit

        while (!it.end()) {
            if (it) {
                combs[bit % comb_columns] |= (static_cast<comb_type>(1) << (bit / comb_columns));
            } else {
                // technically no need to reset the zero bit
                // only intended to minimize the timing attack
                combs[bit % comb_columns] &= ~(static_cast<comb_type>(1) << (bit / comb_columns));
            }
            it.next();
            bit++;
        }

        auto e = (comb_columns >> 1);

        ecp_jac<N> Q{};
        int c = e - 1;

        // loop until we find the first non-zero comb
        // and randomize it
        while (c >= 0) {
            if (combs[c] != 0) {
                Q = ecp_jac<N>(presets[combs[c]]);
                ecbase::randomize(Q);
                point_addition(Q, alt_presets[combs[c + e]], Q);
                break;
            } else if (combs[c + e] != 0) {
                Q = ecp_jac<N>(alt_presets[combs[c + e]]);
                ecbase::randomize(Q);
                point_addition(Q, alt_presets[combs[c]], Q);
                break;
            }
            c--;
        }

        while (c > 0) {
            c--;
            point_double(Q, Q);
            point_addition(Q, presets[combs[c]], Q);
            point_addition(Q, alt_presets[combs[c + e]], Q);
        }

        normalize(Q, r);
    }

    void precompute()
    {
        ecp_aff<N> powers[WSIZE];
        powers[0] = ecp_aff(DP::dp_gx, DP::dp_gy);

        mpi<N> x(1);
        mpi<N> y(1);
        x.lshift(comb_columns);
        y.lshift(comb_columns / 2);

        for (unsigned i = 1; i != WSIZE; i++) {
            point_mul(x, powers[i - 1], powers[i]);
        }

        for (unsigned i = 0; i < (1 << WSIZE); i++) {
            int bit = WSIZE - 1;
            ecp_jac<N> Q;
            for (int bit = WSIZE - 1; bit >= 0; bit--) {
                if ((1 << bit) & i) {
                    point_addition(Q, powers[bit], Q);
                }
            }

            // TODO use batch normalising to save time
            // at the cost of space
            normalize(Q, presets[i]);
            point_mul(y, presets[i], alt_presets[i]);
        }
    }

    void negate(ecp_jac<N>& p) const
    {
        if (!p.Y.equal(0)) {
            p.Y = DP::dp_p - p.Y;
        }
    }
};

struct dp_secp256r1 {
    static inline constexpr char curve_name[] = "secp256r1";

    static inline constexpr CurveId curve_id = CurveId::SECP256R1;

    static inline constexpr CurveType curve_type = CurveType::SHORT_WEIERSTRASS;

    // all domain parameters
    static inline const mpi<256> dp_p{std::array<uint8_t, 32>{0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x01,
                                                              0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                                              0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF,
                                                              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}};

    static inline constexpr int special_a = -3;

    static inline const mpi<256> dp_b{std::array<uint8_t, 32>{0x5A, 0xC6, 0x35, 0xD8, 0xAA, 0x3A, 0x93, 0xE7,
                                                              0xB3, 0xEB, 0xBD, 0x55, 0x76, 0x98, 0x86, 0xBC,
                                                              0x65, 0x1D, 0x06, 0xB0, 0xCC, 0x53, 0xB0, 0xF6,
                                                              0x3B, 0xCE, 0x3C, 0x3E, 0x27, 0xD2, 0x60, 0x4B}};

    static inline const mpi<256> dp_gx{std::array<uint8_t, 32>{0x6B, 0x17, 0xD1, 0xF2, 0xE1, 0x2C, 0x42, 0x47,
                                                               0xF8, 0xBC, 0xE6, 0xE5, 0x63, 0xA4, 0x40, 0xF2,
                                                               0x77, 0x03, 0x7D, 0x81, 0x2D, 0xEB, 0x33, 0xA0,
                                                               0xF4, 0xA1, 0x39, 0x45, 0xD8, 0x98, 0xC2, 0x96}};

    static inline const mpi<256> dp_gy{std::array<uint8_t, 32>{0x4F, 0xE3, 0x42, 0xE2, 0xFE, 0x1A, 0x7F, 0x9B,
                                                               0x8E, 0xE7, 0xEB, 0x4A, 0x7C, 0x0F, 0x9E, 0x16,
                                                               0x2B, 0xCE, 0x33, 0x57, 0x6B, 0x31, 0x5E, 0xCE,
                                                               0xCB, 0xB6, 0x40, 0x68, 0x37, 0xBF, 0x51, 0xF5}};

    static inline const mpi<256> dp_n{std::array<uint8_t, 32>{0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x00,
                                                              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
                                                              0xBC, 0xE6, 0xFA, 0xAD, 0xA7, 0x17, 0x9E, 0x84,
                                                              0xF3, 0xB9, 0xCA, 0xC2, 0xFC, 0x63, 0x25, 0x51}};

    static inline limb_t dp_h = 1;

    static inline constexpr bool sca = true;

    struct uint32calc {
        uint32calc(uint32_t* p) : ptr(p){};

        template <typename... P>
        inline uint32calc& add(uint32_t t, P... p)
        {
            *ptr += t;
            *(ptr + 1) += (*ptr < t);
            if constexpr (sizeof...(p) > 0) {
                add(p...);
            }
            return *this;
        }

        template <typename... P>
        inline uint32calc& sub(uint32_t t, P... p)
        {
            borrow += (*ptr < t);
            *ptr -= t;
            if constexpr (sizeof...(p) > 0) {
                sub(p...);
            }
            return *this;
        }

        inline uint32calc& next()
        {
            ptr++;
            uint32_t tmp = borrow;
            borrow = *ptr < tmp;
            *ptr -= tmp;
            return *this;
        }

        void finalize()
        {
            ptr++;

            if (*ptr >= borrow) {
                *ptr -= borrow;
                borrow = 0;
            }
        }

        uint32_t borrow = 0;
        uint32_t* ptr;
    };

// two implementations below are effectively the same alg
// but first one appears to be faster
#if 1
    // p-256 modulus
    // Guide to Eclliptic Curve Cryptography 2.29
    static auto reduce(mpi<512> const& a)
    {
        auto c = get_mpi32_viewer(a);
        uint32_t A[9] = {c.get<0>(), c.get<1>(), c.get<2>(), c.get<3>(),
                         c.get<4>(), c.get<5>(), c.get<6>(), c.get<7>(), 0};

        uint32calc calc(&A[0]);

        calc.add(c.get<8>(), c.get<9>()) // A[0]
            .sub(c.get<11>(), c.get<12>(), c.get<13>(), c.get<14>())
            .next() // A[1]
            .add(c.get<9>(), c.get<10>())
            .sub(c.get<12>(), c.get<13>(), c.get<14>(), c.get<15>())
            .next() // A[2]
            .add(c.get<10>(), c.get<11>())
            .sub(c.get<13>(), c.get<14>(), c.get<15>())
            .next() // A[3]
            .add(c.get<11>(), c.get<11>(), c.get<12>(), c.get<12>(), c.get<13>())
            .sub(c.get<15>(), c.get<8>(), c.get<9>())
            .next() // A[4]
            .add(c.get<12>(), c.get<12>(), c.get<13>(), c.get<13>(), c.get<14>())
            .sub(c.get<9>(), c.get<10>())
            .next() // A[5]
            .add(c.get<13>(), c.get<13>(), c.get<14>(), c.get<14>(), c.get<15>())
            .sub(c.get<10>(), c.get<11>())
            .next() // A[6]
            .add(c.get<14>(), c.get<14>(), c.get<15>(), c.get<15>(),
                 c.get<14>(), c.get<13>())
            .sub(c.get<8>(), c.get<9>())
            .next() // A[7]
            .add(c.get<15>(), c.get<15>(), c.get<15>(), c.get<8>())
            .sub(c.get<10>(), c.get<11>(), c.get<12>(), c.get<13>())
            .finalize();

        auto s1 = assemble_mpi32(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]);
        if (calc.borrow) {
            auto s2 = assemble_mpi32(0, 0, 0, 0, 0, 0, 0, 0, calc.borrow);
            s2.sub(s1);
            if (a.sign > 0)
                s2.sign = -1;
            return s2;
        } else {
            if (a.sign < 0)
                s1.sign = -1;
            return s1;
        }
    }

#else
    static decltype(auto) reduce(mpi<512> const& a)
    {
        auto c = get_mpi32_viewer(a);

        auto s1 = assemble_mpi32(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);
        auto s2 = assemble_mpi32(0, 0, 0, c[11], c[12], c[13], c[14], c[15]);
        s2.lshift<1>();
        auto s3 = assemble_mpi32(0, 0, 0, c[12], c[13], c[14], c[15], 0);
        s3.lshift<1>();
        auto s4 = assemble_mpi32(c[8], c[9], c[10], 0, 0, 0, c[14], c[15]);
        auto s5 = assemble_mpi32(c[9], c[10], c[11], c[13], c[14], c[15], c[13], c[8]);
        auto s6 = assemble_mpi32(c[11], c[12], c[13], 0, 0, 0, c[8], c[10]);
        auto s7 = assemble_mpi32(c[12], c[13], c[14], c[15], 0, 0, c[9], c[11]);
        auto s8 = assemble_mpi32(c[13], c[14], c[15], c[8], c[9], c[10], 0, c[12]);
        auto s9 = assemble_mpi32(c[14], c[15], 0, c[9], c[10], c[11], 0, c[13]);

        s1.add(s2);
        s1.add(s3);
        s1.add(s4);
        s1.add(s5);

        s6.add(s7);
        s6.add(s8);
        s6.add(s9);

        s1.sub(s6);

        //s1.sign = s1.sign * a.sign;
        if (s1.sign != a.sign) {
            s1.sign = -1;
        } else {
            s1.sign = 1;
        }

        return s1;
    }
#endif

    static inline montg_inverse_helper iop{dp_p};
};

// birational map between this edward curve and montgomery is
// (u, v) = ((1+y)/(1-y), sqrt(486664)*u/x)
// (x, y) = (sqrt(486664)*u/v, (u-1)/(u+1))
//
struct dp_ed25519 {
    static inline constexpr char curve_name[] = "ed25519";

    static inline constexpr CurveId curve_id = CurveId::ED25519;

    static inline constexpr CurveType curve_type = CurveType::EDWARD;

    // all domain parameters
    // P = 2^255 - 19
    static inline const mpi<256> dp_p = (mpi<256>(1) << 255) - mpi<256>(19);

    static inline const mpi<256> sqrt486664{std::array<uint8_t, 32>{
        0x70, 0xD9, 0x12, 0x0B, 0x9F, 0x5F, 0xF9, 0x44,
        0x2D, 0x84, 0xF7, 0x23, 0xFC, 0x03, 0xB0, 0x81,
        0x3A, 0x5E, 0x2C, 0x2E, 0xB4, 0x82, 0xE5, 0x7D,
        0x33, 0x91, 0xFB, 0x55, 0x00, 0xBA, 0x81, 0xE7}};

    static inline constexpr unsigned int p255_width = (255 / bits_per_limb + 1);

    static inline const mpi<256> dp_a = dp_p - mpi<256>(1);

    static inline constexpr int special_a = -1;

    static inline const mpi<256> dp_d{std::array<uint8_t, 32>{
        0x52, 0x03, 0x6C, 0xEE, 0x2B, 0x6F, 0xFE, 0x73,
        0x8C, 0xC7, 0x40, 0x79, 0x77, 0x79, 0xE8, 0x98,
        0x00, 0x70, 0x0A, 0x4D, 0x41, 0x41, 0xD8, 0xAB,
        0x75, 0xEB, 0x4D, 0xCA, 0x13, 0x59, 0x78, 0xA3}};

    static inline const mpi<256> dp_gx{std::array<uint8_t, 32>{
        0x21, 0x69, 0x36, 0xD3, 0xCD, 0x6E, 0x53, 0xFE,
        0xC0, 0xA4, 0xE2, 0x31, 0xFD, 0xD6, 0xDC, 0x5C,
        0x69, 0x2C, 0xC7, 0x60, 0x95, 0x25, 0xA7, 0xB2,
        0xC9, 0x56, 0x2D, 0x60, 0x8F, 0x25, 0xD5, 0x1A}};

    static inline const mpi<256> dp_gy{std::array<uint8_t, 32>{
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x58}};

    // N = 2^252 + 27742317777372353535851937790883648493
    static inline const mpi<256> dp_n = (mpi<256>(1) << 252) +
                                        mpi<128>(std::array<uint8_t, 16>{0x14, 0xDE, 0xF9, 0xDE, 0xA2, 0xF7, 0x9C, 0xD6,
                                                                         0x58, 0x12, 0x63, 0x1A, 0x5C, 0xF5, 0xD3, 0xED});

    static inline constexpr limb_t dp_h = 1;

    static inline constexpr bool sca = false;

    // A = a0 + 2^255 a1
    // return a0 + 19 * a1
    template <unsigned M>
    inline static decltype(auto) internal_reduce(mpi<M> const& A)
    {
        static constexpr size_t hi_width = mpi<M>::num_of_limbs - p255_width + 1;
        mpi<256> a0(&A[0], p255_width);
        a0.set_bit(255, 0);

        mpi<hi_width * bits_per_limb> a1(&A[p255_width - 1], hi_width);
        a1.rshift(255 & bit_pos_mask);
        auto x = a1 * mpi<8>(19);
        return x + a0;
    }

    static decltype(auto) reduce(mpi<512> const& a)
    {
        auto r1 = internal_reduce(a);
        auto r2 = internal_reduce(r1); // further reduction to 256 bit

        r2.sign = a.sign;

        return r2;
    }

    static inline montg_inverse_helper iop{dp_p};
};

struct dp_curve25519 {

    static inline constexpr char curve_name[] = "curve25519";

    static inline constexpr CurveId curve_id = CurveId::CURVE25519;

    static inline constexpr CurveType curve_type = CurveType::MONTGOMERY;

    // all domain parameters
    // P = 2^255 - 19
    static inline const mpi<256> dp_p = (mpi<256>(1) << 255) - mpi<256>(19);

    static inline constexpr unsigned int p255_width = (255 / bits_per_limb + 1);

    // (486662 + 2) / 4
    static inline constexpr limb_t A = (486662 + 2) / 4;

    static inline constexpr int _A = 486662;

    static inline constexpr int nbits = 254;

    static inline const mpi<256> dp_a{std::array<uint8_t, 3>{0x01, 0xDB, 0x42}};

    static inline const mpi<256> dp_gx{9};
    static inline const mpi<256> dp_gy{std::array<uint8_t, 32>{
        0x20, 0xAE, 0x19, 0xA1, 0xB8, 0xA0, 0x86, 0xB4,
        0xE0, 0x1E, 0xDD, 0x2C, 0x77, 0x48, 0xD1, 0x4C,
        0x92, 0x3D, 0x4D, 0x7E, 0x6D, 0x7C, 0x61, 0xB2,
        0x29, 0xE9, 0xC5, 0xA2, 0x7E, 0xCE, 0xD3, 0xD9}};

    // N = 2^252 + 27742317777372353535851937790883648493
    static inline const mpi<256> dp_n = (mpi<256>(1) << 252) +
                                        mpi<128>(std::array<uint8_t, 16>{0x14, 0xDE, 0xF9, 0xDE, 0xA2, 0xF7, 0x9C, 0xD6,
                                                                         0x58, 0x12, 0x63, 0x1A, 0x5C, 0xF5, 0xD3, 0xED});

    static inline constexpr limb_t dp_h = 1;

    static inline constexpr bool sca = false;

    // A = a0 + 2^255 a1
    // return a0 + 19 * a1
    template <unsigned M>
    inline static decltype(auto) internal_reduce(mpi<M> const& A)
    {
        static constexpr size_t hi_width = mpi<M>::num_of_limbs - p255_width + 1;
        mpi<256> a0(&A[0], p255_width);
        a0.set_bit(255, 0);

        mpi<hi_width * bits_per_limb> a1(&A[p255_width - 1], hi_width);
        a1.rshift(255 & bit_pos_mask);
        auto x = a1 * mpi<5>(19);
        return x + a0;
    }

    static decltype(auto) reduce(mpi<512> const& a)
    {
        auto r1 = internal_reduce(a);
        auto r2 = internal_reduce(r1); // further reduction to 256 bit

        r2.sign = a.sign;

        return r2;
    }

    static inline montg_inverse_helper iop{dp_p};
};

struct dp_curve448 {

    static inline constexpr char curve_name[] = "curve448";

    static inline constexpr CurveId curve_id = CurveId::CURVE448;

    static inline constexpr CurveType curve_type = CurveType::MONTGOMERY;

    // all domain parameters
    // P = 2^448 - 2^224 - 1
    static inline const mpi<448> dp_p = resize<448>((mpi<449>(1) << 448) - (mpi<448>(1) << 224) - mpi<1>(1));

    // (156326 + 2) / 4
    static inline constexpr limb_t A = (156326 + 2) / 4;

    static inline constexpr int _A = 156326;

    static inline constexpr int nbits = 447;

    static inline const mpi<448> dp_a{std::array<uint8_t, 2>{0x98, 0xAA}};

    static inline const mpi<448> dp_gx{5};
    static inline const mpi<448> dp_gy{std::array<uint8_t, 56>{
        0x7D, 0x23, 0x5D, 0x12, 0x95, 0xF5, 0xB1, 0xF6,
        0x6C, 0x98, 0xAB, 0x6E, 0x58, 0x32, 0x6F, 0xCE,
        0xCB, 0xAE, 0x5D, 0x34, 0xF5, 0x55, 0x45, 0xD0,
        0x60, 0xF7, 0x5D, 0xC2, 0x8D, 0xF3, 0xF6, 0xED,
        0xB8, 0x02, 0x7E, 0x23, 0x46, 0x43, 0x0D, 0x21,
        0x13, 0x12, 0xC4, 0xB1, 0x50, 0x67, 0x7A, 0xF7,
        0x6F, 0xD7, 0x22, 0x3D, 0x45, 0x7B, 0x5B, 0x1A}};

    // N = 2^446 - 0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
    static inline const mpi<448> dp_n = (mpi<448>(1) << 446) -
                                        mpi<224>(std::array<uint8_t, 28>{0x83, 0x35, 0xdc, 0x16,
                                                                         0x3b, 0xb1, 0x24, 0xb6,
                                                                         0x51, 0x29, 0xc9, 0x6f,
                                                                         0xde, 0x93, 0x3d, 0x8d,
                                                                         0x72, 0x3a, 0x70, 0xaa,
                                                                         0xdc, 0x87, 0x3d, 0x6d,
                                                                         0x54, 0xa7, 0xbb, 0x0d});

    static inline constexpr limb_t dp_h = 4;

    static inline constexpr bool sca = false;

    // A = A0 + 2^448 A1, A1 = B0 + 2^224 B1
    // return A0 + A1 + B1 + (B0 + B1) * 2^224
    template <unsigned M>
    inline static decltype(auto) internal_reduce(mpi<M> const& A)
    {
        if constexpr (M <= 448) {
            return A;
        } else {
            mpi<448 + 1> a0(&A[0], 448 / bits_per_limb);
            mpi<448> a1(&A[448 / bits_per_limb], A.num_of_limbs - 448 / bits_per_limb);
            a0.add(a1);

            if constexpr (M - 448 > 224) {
                mpi<a1.num_of_bits> b1 = (a1 >> 224);
                a0.add(b1);
                b1.lshift(224);
                a0.add(b1);
            }

            a1.lshift(224);
            a0.add(a1);
            return a0;
        }
    }

    static decltype(auto) reduce(mpi<448 * 2> const& a)
    {
        auto r1 = internal_reduce(a);
        auto r2 = internal_reduce(r1); // further reduction to 256 bit

        r2.sign = a.sign;
        return r2;
    }

    static inline montg_inverse_helper iop{dp_p};
};

struct dp_ed448 {
    static inline constexpr char curve_name[] = "ed448";

    static inline constexpr CurveId curve_id = CurveId::ED448;

    static inline constexpr CurveType curve_type = CurveType::EDWARD;

    // all domain parameters
    // P = 2^448 - 2^224 - 1
    static inline const mpi<448> dp_p = resize<448>((mpi<449>(1) << 448) - (mpi<448>(1) << 224) - mpi<1>(1));

    static inline const mpi<256> sqrt486664{std::array<uint8_t, 32>{
        0x70, 0xD9, 0x12, 0x0B, 0x9F, 0x5F, 0xF9, 0x44,
        0x2D, 0x84, 0xF7, 0x23, 0xFC, 0x03, 0xB0, 0x81,
        0x3A, 0x5E, 0x2C, 0x2E, 0xB4, 0x82, 0xE5, 0x7D,
        0x33, 0x91, 0xFB, 0x55, 0x00, 0xBA, 0x81, 0xE7}};

    static inline constexpr unsigned int p448_width = (448 / bits_per_limb);

    static inline const mpi<448> dp_a = mpi<448>(1);

    static inline constexpr int special_a = 1;

    static inline constexpr uint8_t phflag = 0;

    static inline const mpi<448> dp_d = dp_p - mpi<32>(39081);

    // 224580040295924300187604334099896036246789641632564134246125461686950415467406032909029192869357953282578032075146446173674602635247710
    static inline const mpi<448> dp_gx{std::array<uint8_t, 56>{
        0x4F, 0x19, 0x70, 0xC6, 0x6B, 0xED, 0x0D, 0xED,
        0x22, 0x1D, 0x15, 0xA6, 0x22, 0xBF, 0x36, 0xDA,
        0x9E, 0x14, 0x65, 0x70, 0x47, 0x0F, 0x17, 0x67,
        0xEA, 0x6D, 0xE3, 0x24, 0xA3, 0xD3, 0xA4, 0x64,
        0x12, 0xAE, 0x1A, 0xF7, 0x2A, 0xB6, 0x65, 0x11,
        0x43, 0x3B, 0x80, 0xE1, 0x8B, 0x00, 0x93, 0x8E,
        0x26, 0x26, 0xA8, 0x2B, 0xC7, 0x0C, 0xC0, 0x5E}};

    // 298819210078481492676017930443930673437544040154080242095928241372331506189835876003536878655418784733982303233503462500531545062832660
    static inline const mpi<448> dp_gy{std::array<uint8_t, 56>{
        0x69, 0x3F, 0x46, 0x71, 0x6E, 0xB6, 0xBC, 0x24,
        0x88, 0x76, 0x20, 0x37, 0x56, 0xC9, 0xC7, 0x62,
        0x4B, 0xEA, 0x73, 0x73, 0x6C, 0xA3, 0x98, 0x40,
        0x87, 0x78, 0x9C, 0x1E, 0x05, 0xA0, 0xC2, 0xD7,
        0x3A, 0xD3, 0xFF, 0x1C, 0xE6, 0x7C, 0x39, 0xC4,
        0xFD, 0xBD, 0x13, 0x2C, 0x4E, 0xD7, 0xC8, 0xAD,
        0x98, 0x08, 0x79, 0x5B, 0xF2, 0x30, 0xFA, 0x14}};

    // N = 2^446 - 13818066809895115352007386748515426880336692474882178609894547503885
    static inline const mpi<448> dp_n{std::array<uint8_t, 56>{
        0x3f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0x7c, 0xca, 0x23, 0xe9,
        0xc4, 0x4e, 0xdb, 0x49, 0xae, 0xd6, 0x36, 0x90,
        0x21, 0x6c, 0xc2, 0x72, 0x8d, 0xc5, 0x8f, 0x55,
        0x23, 0x78, 0xc2, 0x92, 0xab, 0x58, 0x44, 0xf3}};

    static inline constexpr bool sca = false;

    // N = A0 + 2^448 A1
    // A1 = B0 + 2^224 B1
    // return A0 + A1 + B1 + (B0 + B1) * 2^224
    template <unsigned M>
    inline static decltype(auto) internal_reduce(mpi<M> const& A)
    {
        if constexpr (M <= 448) {
            return A;
        } else {
            mpi<448 + 1> a0(&A[0], 448 / bits_per_limb);
            mpi<448> a1(&A[448 / bits_per_limb], A.num_of_limbs - 448 / bits_per_limb);
            a0.add(a1);

            if constexpr (M - 448 > 224) {
                mpi<a1.num_of_bits> b1 = (a1 >> 224);
                a0.add(b1);
                b1.lshift(224);
                a0.add(b1);
            }

            a1.lshift(224);
            a0.add(a1);
            return a0;
        }
    }

    static decltype(auto) reduce(mpi<448 * 2> const& a)
    {
        auto r1 = internal_reduce(a);
        auto r2 = internal_reduce(r1);

        r2.sign = a.sign;

        return resize<448>(r2);
    }

    static inline montg_inverse_helper iop{dp_p};
};

template <CurveType CT, typename DP, unsigned N>
class ec_alg : public ec_curve<CT, DP, N>
{
    using eccurve = ec_curve<CT, DP, N>;
    using eccurve::ecbase::ffp_reduce;

public:
    using dp_template = DP;
    using keypair = typename eccurve::keypair;
    using pubkey_t = typename eccurve::pubkey_t;
    using prvkey_t = typename eccurve::prvkey_t;
    using signature_t = typename eccurve::signature_t;
    using eccurve::scalar_mul;
    using shared_secret_t = mpi<N>;

    ec_alg(randomizer* rng) : eccurve()
    {
        eccurve::ecbase::rng = rng;
    }

    inline const mpi<N>& dp_p() const
    {
        return DP::dp_p;
    }

    inline void reduce(mpi<N * 2> const& a, mpi<N>& r) const
    {
        ffp_reduce(a, r);
    }

    bool ecdh_compute_shared(prvkey_t const& d, pubkey_t const& Q, shared_secret_t& secret) const
    {
        ecp_aff<N> R;
        bool ret = scalar_mul(d, Q, R);
        if (ret) {
            secret = R.X;
        }
        return ret;
    }

    inline bool ecdh_compute_shared(keypair const& key, pubkey_t const& Q, shared_secret_t& secret) const
    {
        return ecdh_compute_shared(key.d, Q, secret);
    }
};

using curve25519 = ec_alg<CurveType::MONTGOMERY, dp_curve25519, 256>;
using curve448 = ec_alg<CurveType::MONTGOMERY, dp_curve448, 448>;
using secp256r1 = ec_alg<CurveType::SHORT_WEIERSTRASS, dp_secp256r1, 256>;
using ed25519 = ec_alg<CurveType::EDWARD, dp_ed25519, 256>;
using ed448 = ec_alg<CurveType::EDWARD, dp_ed448, 448>;

}; // namespace ecc_crypto

#endif
