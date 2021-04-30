// SPDX-License-Identifier: GPL-3.0-only
/*
 *  mpi.h
 *
 *  Copyright (C) 2019 Daniel Hu <daddy.of.qq@gmail.com>
 *
 *  multi-precision math in mordern C++ (need c++17 or later)
 */

#ifndef MPI_H
#define MPI_H

#include <algorithm>
#include <array>
#include <cstdint>

#ifdef DEBUG_CAST
#include <exception>
#endif

#ifndef CACHE_LINE_ALIGN
#define CACHE_LINE_ALIGN 32
#endif

template <typename Dest, typename Src>
inline Dest narrow_cast(Src src)
{
    Dest d = static_cast<Dest>(src);
#ifdef DEBUG_CAST
    if (src != static_cast<Src>(d)) {
        throw std::bad_cast();
    }
#endif
    return d;
};

namespace mpi_ns
{

using limb_t = size_t;
static inline constexpr int bytes_per_limb = sizeof(limb_t);
static inline constexpr unsigned bits_per_limb = (bytes_per_limb << 3);
static inline constexpr unsigned bit_pos_mask = (bits_per_limb - 1);

constexpr unsigned compute_limb_pos_shift()
{
    unsigned shift = 0;
    auto mask = bit_pos_mask;
    while (mask) {
        shift++;
        mask >>= 1;
    }
    return shift;
};
static inline constexpr unsigned limb_pos_shift = compute_limb_pos_shift();

template <typename T, typename... P>
inline void longadd(limb_t& carry, limb_t& sum, T t, P... p)
{
    sum += t;
    carry += (sum < t);
    if constexpr (sizeof...(p) > 0) {
        longadd(carry, sum, p...);
    }
}

template <unsigned S>
struct double_precision;

template <>
struct double_precision<2> {
    using type = uint32_t;
};

template <>
struct double_precision<4> {
    using type = uint64_t;
};

// for a 32-bit CPU (like ARM v7, 128 bit is neither defined, but not needed for the implementation)
// for 64-bit cpu (like x86-64 or ARM v8, we need 128 bit to perform efficient multiplication)
template <>
struct double_precision<8> {
#if defined(__SIZEOF_INT128__)
    using type = __uint128_t;
#endif
};

using double_precision_t = typename double_precision<sizeof(limb_t)>::type;

// multiplication with optional shift
// we are assuming the platform has double-precision multiplication
// this is not the best code for best performance
// but it serve our need without writing assembly
template <unsigned SHIFT = 0>
static inline limb_t longmul(limb_t a, limb_t b, limb_t& hi, limb_t& lo)
{
    limb_t ret{0};
    double_precision_t result = static_cast<double_precision_t>(a) * b;
    if constexpr (SHIFT > 0) {
        ret = static_cast<limb_t>(result >> (sizeof(result) * 8 - SHIFT));
        result <<= SHIFT;
    };
    lo = static_cast<limb_t>(result);
    hi = static_cast<limb_t>(result >> bits_per_limb);
    return ret;
}

static inline limb_t longdiv(limb_t& hi, limb_t& lo, limb_t d)
{
    limb_t ret{0};
    double_precision_t result = (static_cast<double_precision_t>(hi) << (bits_per_limb)) | static_cast<double_precision_t>(lo);
    result /= d;
    return static_cast<limb_t>(result);
}

// TODO do we need a side-channel safe comparison??
static inline int compare_helper(const limb_t* p1, size_t len1,
                                 const limb_t* p2, size_t len2)
{
    p1 += len1;
    p2 += len2;

    while (len1 > len2) {
        if (*--p1 != 0)
            return 1;
        len1--;
    }

    while (len2 > len1) {
        if (*--p2 != 0)
            return -1;
        len2--;
    }

    while (len1 && *--p1 == *--p2) {
        len1--;
    }

    if (len1 == 0) {
        return 0;
    } else if (*p1 > *p2) {
        return 1;
    } else {
        return -1;
    }
}

static inline limb_t add_helper(limb_t* d, const limb_t* a, const limb_t* b, size_t size, size_t extra)
{
    unsigned c = 0;

    while (size--) {
        limb_t x = *a + c;
        c = (x < c);
        x += *b;
        c += (x < *b);
        a++;
        b++;
        *d++ = x;
    }

    while (c && extra--) {
        limb_t x = *a + c;
        c = (x < c);
        a++;
        *d++ = x;
    };

    return c;
}

static inline limb_t sub_helper(limb_t* d, const limb_t* a, const limb_t* b, size_t size, size_t extra)
{
    unsigned c = 0, z;

    while (size--) {
        z = (*a < c);
        limb_t x = *a++ - c;

        c = (x < *b) + z;
        *d++ = x - *b++;
    }

    while (c && extra--) {
        z = (*a < c);
        *d++ = *a++ - c;
        c = z;
    }

    return c;
}

static limb_t inc_mul_helper(limb_t x, const limb_t* src, size_t size, limb_t* dest, size_t extra, limb_t carry = 0)
{
    while (size--) {
        limb_t high, low;
        longmul(x, *src, high, low);
        limb_t sum{0};
        longadd(high, sum, low, *dest, carry);
        *dest++ = sum;
        carry = high;
        src++;
    };

    while (carry && extra--) {
        *dest += carry;
        carry = (*dest < carry);
        dest++;
    };

    return carry;
};

// the only difference with above is no addition to destination
static limb_t mul_helper(limb_t x, const limb_t* src, size_t size, limb_t* dest)
{
    limb_t carry{0};

    while (size--) {
        limb_t high, low;
        longmul(x, *src, high, low);
        limb_t sum{0};
        longadd(high, sum, low, carry);
        *dest++ = sum;
        carry = high;
        src++;
    };

    return carry;
};

template <unsigned N, unsigned M>
struct Copyable {
    const static bool value = (M <= N);
};

struct randomizer {
    virtual void operator()(uint8_t*, size_t) const = 0;
};

template <unsigned SIZE>
class FixedSizeBuf
{
public:
    static inline constexpr size_t buf_size = SIZE;
    using iterator = uint8_t*;
    using const_iterator = const uint8_t*;

    FixedSizeBuf()
    {
        for (auto& x : buf)
            x = 0;
    }

    // assume the ptr will always point to same amount of data
    FixedSizeBuf(const uint8_t* ptr)
    {
        unsigned i = 0;
        while (i != buf_size) {
            buf[i++] = *ptr++;
        }
    }

    void set(const uint8_t* ptr)
    {
        unsigned i = 0;
        while (i != buf_size) {
            buf[i++] = *ptr++;
        }
    }

    FixedSizeBuf(FixedSizeBuf<SIZE> const& other)
    {
        unsigned i = 0;
        while (i != buf_size) {
            buf[i] = other[i];
            i++;
        }
    }

    FixedSizeBuf& operator=(FixedSizeBuf<SIZE> const& other)
    {
        unsigned i = 0;
        while (i != buf_size) {
            buf[i] = other[i];
            i++;
        }
        return *this;
    }

    inline uint8_t& operator[](size_t index)
    {
        return buf[index];
    }

    inline const uint8_t& operator[](size_t index) const
    {
        return buf[index];
    }

    static inline constexpr size_t size()
    {
        return SIZE;
    }

    inline void copy(const uint8_t* ptr, size_t len, size_t index = 0)
    {
        while (index < SIZE && len--) {
            buf[index++] = *ptr++;
        }
    }

    template <unsigned S>
    inline void copy(FixedSizeBuf<S> const& other, size_t index = 0)
    {
        copy(&other[0], other.size(), index);
    }

    auto data()
    {
        return &(*this)[0];
    }

    const uint8_t* begin() const
    {
        return &buf[0];
    }

    const uint8_t* end() const
    {
        return &buf[SIZE];
    }

    uint8_t* begin()
    {
        return &buf[0];
    }

    uint8_t* end()
    {
        return &buf[SIZE];
    }

    const_iterator cbegin() const
    {
        return &(*this)[0];
    }

    const_iterator cend() const
    {
        return &(*this)[SIZE];
    }

    ~FixedSizeBuf()
    {
        for (auto& x : buf)
            x = 0;
    }

    friend inline bool operator==(FixedSizeBuf const& buf1, FixedSizeBuf const& buf2)
    {
        for (unsigned i = 0; i != SIZE; i++) {
            if (buf1.buf[i] != buf2.buf[i]) {
                return false;
            }
        }

        return true;
    }

private:
    uint8_t buf[SIZE];
};

template <unsigned N> // N is the number of bits in bignum
class mpi
{
public:
    static constexpr unsigned num_of_bits = N;
    static constexpr unsigned num_of_bytes = (num_of_bits + 7) / 8;
    static constexpr unsigned num_of_limbs = (num_of_bytes + sizeof(limb_t) - 1) / sizeof(limb_t);

    explicit mpi(limb_t v = 0) : sign(1)
    {
        zeroize();
        buff[0] = v;
    };

    explicit mpi(int v) : sign(v >= 0 ? 1 : -1)
    {
        zeroize();
        buff[0] = v >= 0 ? v : -v;
    };

    inline void zeroize()
    {
        for (auto& x : buff) {
            x = 0;
        }
    }

    mpi(const uint8_t* ptr, size_t size) : sign(1)
    {
        zeroize();
        import(ptr, size > num_of_bytes ? num_of_bytes : size);
    };

    template <unsigned long int ARRSIZE>
    explicit mpi(std::array<uint8_t, ARRSIZE> const& a) : sign(1)
    {
        static_assert(ARRSIZE * 8 <= N, "array length must be less than mpi");
        import(a.data(), a.size());
    };

    template <unsigned ARRSIZE>
    explicit mpi(std::array<uint8_t, ARRSIZE> const& a) : sign(1)
    {
        static_assert(ARRSIZE * 8 <= N, "array length must be less than mpi");
        import(a.data(), a.size());
    };

    mpi(const limb_t* ptr, size_t limbs, int offset = 0) : sign(1)
    {
        for (unsigned i = 0; i < num_of_limbs; i++) {
            if (i >= offset && limbs) {
                buff[i] = *ptr++;
                limbs--;
            } else {
                buff[i] = 0;
            }
        }
    }

    template <unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
    explicit mpi(mpi<M> const& other) : sign(other.sign)
    {
        auto p = &buff[0];
        auto q = &other[0];

        for (unsigned i = mpi<M>::num_of_limbs; i; i--) {
            *p++ = *q++;
        }

        for (unsigned i = num_of_limbs - mpi<M>::num_of_limbs; i; i--) {
            *p++ = 0;
        }

        return;
    };

    template <unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
    mpi<N>& operator=(mpi<M> const& other)
    {
        auto p = &buff[0];
        auto q = &other.buff[0];

        sign = other.sign;

        for (auto i = mpi<M>::num_of_limbs; i; i--) {
            *p++ = *q++;
        }

        for (auto i = num_of_limbs - mpi<M>::num_of_limbs; i; i--) {
            *p++ = 0;
        }

        return *this;
    };

    inline void negate()
    {
        sign = -sign;
    }

    inline void randomize(randomizer* rng)
    {
        (*rng)(reinterpret_cast<uint8_t*>(&buff[0]), num_of_bytes);
    }

    auto output(bool be = true) const
    {
        FixedSizeBuf<num_of_bytes> buf;
        output(&buf[0], buf.size(), be);
        return buf;
    }

    void output(uint8_t* ptr, size_t size, bool be = true) const
    {
        if (be)
            ptr += size - 1;
        auto left = size > num_of_bytes ? num_of_bytes : size;
        auto p = &buff[0];

        while (left >= sizeof(limb_t)) {
            auto x = *p;
            for (unsigned i = sizeof(limb_t); i; i--) {
                *ptr = static_cast<uint8_t>(x & 0xFF);
                if (be)
                    ptr--;
                else
                    ptr++;
                x >>= 8;
            }

            left -= sizeof(limb_t);
            p++;
        };

        auto remain = *p;
        while (left--) {
            *ptr = remain & 0xFF;
            if (be)
                ptr--;
            else
                ptr++;
            remain >>= 8;
        }
    }

    inline limb_t& operator[](int index)
    {
        return buff[index];
    }

    inline const limb_t& operator[](int index) const
    {
        return buff[index];
    }

    inline bool is_odd() const
    {
        return ((buff[0] & 1) != 0);
    }

    inline bool is_even() const
    {
        return !is_odd();
    }

    bool equal(limb_t v) const
    {
        if (buff[0] != v)
            return false;

        for (unsigned i = 1; i < num_of_limbs; i++) {
            if (buff[i])
                return false;
        }

        return true;
    };

    // great than or equal
    template <unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
    int compare(const mpi<M>& other) const
    {
        auto diff = sign - other.sign;
        if (diff != 0)
            return diff;

        int ret = compare_helper(&buff[0], num_of_limbs,
                                 &other[0], mpi<M>::num_of_limbs);

        //return ret * sign;
        if (sign > 0) {
            return ret;
        }

        return -ret;
    }

    template <unsigned K>
    void rshift()
    {
        constexpr unsigned limbs = K >> limb_pos_shift;
        constexpr unsigned shift1 = K & bit_pos_mask;

        if (limbs) {
            for (unsigned i = 0; i != num_of_limbs; i++) {
                if (i + limbs < num_of_limbs) {
                    buff[i] = buff[i + limbs];
                } else {
                    buff[i] = 0;
                }
            }
        }

        if (shift1) {
            constexpr unsigned shift2 = bits_per_limb - shift1;
            limb_t high{0};
            unsigned i = num_of_limbs;
            auto p = &buff[i];
            while (i--) {
                auto low = *--p;
                *p = (high << shift2) | (low >> shift1);
                high = low;
            }
        }
    };

    void rshift(unsigned count)
    {
        unsigned limbs = count >> limb_pos_shift;
        unsigned shift1 = count & bit_pos_mask;

        if (limbs) {
            for (unsigned i = 0; i != num_of_limbs; i++) {
                if (i + limbs < num_of_limbs) {
                    buff[i] = buff[i + limbs];
                } else {
                    buff[i] = 0;
                }
            }
        }

        if (shift1) {
            unsigned shift2 = bits_per_limb - shift1;
            limb_t high{0};
            unsigned i = num_of_limbs;
            auto p = &buff[i];
            while (i--) {
                auto low = *--p;
                *p = (high << shift2) | (low >> shift1);
                high = low;
            }
        }
    };

    template <unsigned K>
    void lshift()
    {
        constexpr unsigned limbs = K >> limb_pos_shift;
        constexpr unsigned shift1 = K & bit_pos_mask;

        if (limbs) {
            for (unsigned i = num_of_limbs; i;) {
                if (--i >= limbs) {
                    buff[i] = buff[i - limbs];
                } else {
                    buff[i] = 0;
                }
            }
        }

        if (shift1) {
            constexpr unsigned shift2 = bits_per_limb - shift1;
            limb_t low{0};
            unsigned i = num_of_limbs;
            auto p = &buff[0];

            while (i--) {
                auto high = *p;
                *p++ = (low >> shift2) | (high << shift1);
                low = high;
            }
        }
    };

    void lshift(unsigned count)
    {
        unsigned limbs = count >> limb_pos_shift;
        unsigned shift1 = count & bit_pos_mask;

        if (limbs) {
            for (unsigned i = num_of_limbs; i;) {
                if (--i >= limbs) {
                    buff[i] = buff[i - limbs];
                } else {
                    buff[i] = 0;
                }
            }
        }

        if (shift1) {
            unsigned shift2 = bits_per_limb - shift1;
            limb_t low{0};
            unsigned i = num_of_limbs;
            auto p = &buff[0];

            while (i--) {
                auto high = *p;
                *p++ = (low >> shift2) | (high << shift1);
                low = high;
            }
        }
    };

    template <unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
    limb_t add(const mpi<M>& other)
    {
        if (sign == other.sign)
            return add_helper(&buff[0], &buff[0], &other[0], mpi<M>::num_of_limbs, num_of_limbs - mpi<M>::num_of_limbs);

        if (compare_helper(&buff[0], num_of_limbs,
                           &other[0], mpi<M>::num_of_limbs) >= 0) {
            return sub_helper(&buff[0], &buff[0], &other[0], mpi<M>::num_of_limbs, num_of_limbs - mpi<M>::num_of_limbs);
        }

        sign = other.sign;
        return sub_helper(&buff[0], &other[0], &buff[0], mpi<M>::num_of_limbs, 0);
    }

    limb_t add(limb_t c)
    {
        auto p = &buff[0];
        for (auto i = num_of_limbs; c && i; i--, p++) {
            *p += c;
            c = (*p < c);
        };

        return c;
    };

    template <unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
    limb_t sub(const mpi<M>& other)
    {
        if (sign != other.sign)
            return add_helper(&buff[0], &buff[0], &other[0], mpi<M>::num_of_limbs, num_of_limbs - mpi<M>::num_of_limbs);

        if (compare_helper(&buff[0], num_of_limbs,
                           &other[0], mpi<M>::num_of_limbs) >= 0) {
            return sub_helper(&buff[0], &buff[0], &other[0], mpi<M>::num_of_limbs, num_of_limbs - mpi<M>::num_of_limbs);
        } else {
            sign = -sign;
            return sub_helper(&buff[0], &other[0], &buff[0], mpi<M>::num_of_limbs, 0);
        }
    }

    int bitlen() const
    {
        int pos = N;
        while (--pos > 0) {
            limb_t mask = static_cast<limb_t>(1) << (pos & bit_pos_mask);

            if (buff[pos >> limb_pos_shift] & mask)
                break;
        }
        return pos + 1;
    }

    void set_bit(int pos, int v)
    {
        limb_t mask = static_cast<limb_t>(1) << (pos & bit_pos_mask);

        if (v) {
            buff[pos >> limb_pos_shift] |= mask;
        } else {
            buff[pos >> limb_pos_shift] &= ~mask;
        }
    }

    int get_bit(int pos) const
    {
        limb_t mask = static_cast<limb_t>(1) << (pos & bit_pos_mask);

        return !!(buff[pos >> limb_pos_shift] & mask);
    }

    void import(const uint8_t* ptr, size_t size)
    {
        auto p = &buff[0];
        for (; size >= sizeof(limb_t); p++) {
            auto _ptr = ptr + size - sizeof(limb_t);
            for (unsigned i = sizeof(limb_t); i; i--) {
                *p <<= 8;
                *p |= *_ptr++;
            };
            size -= sizeof(limb_t);
        }

        while (size--) {
            *p <<= 8;
            *p |= *ptr++;
        };
    };

public:
    alignas(CACHE_LINE_ALIGN) limb_t buff[num_of_limbs];
    int sign;
};

// resize can make copy of mpi to create arbitrarily new size (cut or extend)
// this is not the same like copy constructor of mpi, where you can only copy smaller
// because we want programmer to know what they are doing
template <unsigned M, unsigned N>
static inline void resize(mpi<N> const& src, mpi<M>& dest)
{
    constexpr int limbs = ((mpi<N>::num_of_limbs > mpi<M>::num_of_limbs) ? mpi<M>::num_of_limbs : mpi<N>::num_of_limbs);
    dest.sign = src.sign;

    unsigned i;
    for (i = 0; i < limbs; i++) {
        dest[i] = src[i];
    }

    while (i < dest.num_of_limbs) {
        dest[i++] = 0;
    }
};

template <unsigned M, unsigned N>
static inline mpi<M> resize(mpi<N> const& m)
{
    mpi<M> r{};
    resize(m, r);
    return r;
};

template <typename MPI>
struct bit_iterator {
    bit_iterator(MPI& n, bool left_to_right = true) : n(n)
    {
        if (left_to_right) {
            pos = bits_per_limb * n.num_of_limbs - 1;
            step = -1;
        } else {
            pos = 0;
            step = 1;
        }
    };

    bit_iterator& operator=(bit_iterator const& other) = default;

    void next() { pos += step; };
    void prev() { pos -= step; };
    bool end() const
    {
        return pos < 0 || pos > (bits_per_limb * n.num_of_limbs - 1);
    }

    inline int get_pos() const { return pos; };

    // for efficiency, we don't check boundary here
    // so it relies on caller to check whether they go beyond end
    operator bool() const
    {
        return n[pos >> limb_pos_shift] & (static_cast<limb_t>(1) << (pos & bit_pos_mask));
    };

    void set(bool v)
    {
        limb_t mask = static_cast<limb_t>(1) << (pos & bit_pos_mask);

        if (v) {
            n[pos >> limb_pos_shift] |= mask;
        } else {
            n[pos >> limb_pos_shift] &= ~mask;
        }
    }

    MPI& n;
    int pos;
    int step;
};

template <unsigned N>
struct secure_mpi : public mpi<N> {
    using mpi<N>::mpi;
    ~secure_mpi()
    {
        for (auto& x : mpi<N>::buff) {
            x = 0;
        }
    }
};

template <unsigned N, unsigned M, unsigned K, typename V = std::enable_if_t<Copyable<K, N + M>::value>>
void mul(mpi<N> const& A, mpi<M> const& B, mpi<K>& C)
{
    C.zeroize();
    C.sign = (A.sign == B.sign) ? 1 : -1;

    for (unsigned i = 0; i < A.num_of_limbs; i++) {
        limb_t u{0};
        for (unsigned j = 0; j < B.num_of_limbs; j++) {
            // (UV) <- C[i + j] + A[i]*B[j] + U
            limb_t high, low;
            longmul(A[i], B[j], high, low);
            limb_t sum{0};
            longadd(high, sum, low, C[i + j], u);
            u = high;
            C[i + j] = sum;
        };
        C[i + B.num_of_limbs] = u;
    }
}

// unsafe version of lshift without carry
template <unsigned N>
static inline mpi<N> operator<<(mpi<N> const& A, unsigned int count)
{
    mpi<N> C(A);
    C.lshift(count);
    return C;
}

template <unsigned N>
static inline mpi<N> operator>>(mpi<N> const& A, unsigned int count)
{
    mpi<N> C(A);
    C.rshift(count);
    return C;
}

template <unsigned N, unsigned M>
static inline mpi<N + M> operator*(mpi<N> const& A, mpi<M> const& B)
{
    mpi<N + M> C;
    mul(A, B, C);
    return C;
}

// X = q*Y + r
template <unsigned N, unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
void div(mpi<N> const& X, mpi<M> const& Y, mpi<N>& q, mpi<M>& r)
{
    q.zeroize();
    mpi<N> x(X);

    bit_iterator yi(Y);
    bit_iterator xi(x);
    int delta = (X.num_of_limbs - Y.num_of_limbs) * bits_per_limb;

    while (!yi.end() && !yi)
        yi.next();

    while (1) {
        while (!xi.end() && !xi)
            xi.next();

        if (xi.end())
            break;

        int diff = xi.get_pos() - yi.get_pos();
        // continue normalization of X until
        // we have X within delta bits of Y
        if (diff > delta) {
            diff--;
        } else if (diff == delta) {
            mpi<N> xy(Y);
            xy.lshift(delta);
            if (x.compare(xy) < 0)
                break;
        } else {
            break;
        }

        auto bits = (diff & bit_pos_mask);
        int pos = (diff >> limb_pos_shift);
        limb_t addition = static_cast<limb_t>(1) << bits;
        // q += (1 << diff)
        add_helper(&q[pos], &q[pos], &addition, 1, q.num_of_limbs - pos - 1);

        mpi<M + bits_per_limb> y(Y);
        y.lshift(bits);
        // x -= y * (1 << diff)
        sub_helper(&x[pos], &x[pos], &y[0], y.num_of_limbs, x.num_of_limbs - y.num_of_limbs);
    }

    for (int i = X.num_of_limbs - 1; i >= Y.num_of_limbs; i--) {
        if (x[i] == Y[Y.num_of_limbs - 1]) {
            q[i - Y.num_of_limbs] = ~((limb_t)0);
        } else {
            q[i - Y.num_of_limbs] = longdiv(x[i], x[i - 1], Y[Y.num_of_limbs - 1]);
        }

        do {
            mpi<3 * bits_per_limb> t1;
            t1[0] = Y.num_of_limbs > 1 ? Y[Y.num_of_limbs - 2] : 0;
            t1[1] = Y[Y.num_of_limbs - 1];
            t1[2] = 0;
            mul_helper(q[i - Y.num_of_limbs], &t1[0], 3, &t1[0]);
            mpi<3 * bits_per_limb> t2;
            t2[0] = i > 1 ? x[i - 2] : 0;
            t2[1] = x[i - 1];
            t2[2] = x[i];
            if (t1.compare(t2) > 0) {
                q[i - Y.num_of_limbs]--;
            } else {
                break;
            }
        } while (1);

        mpi<N> yy;

        for (auto k = 0; k != Y.num_of_limbs; k++) {
            yy[k + i - Y.num_of_limbs] = Y[k];
        }

        mul_helper(q[i - Y.num_of_limbs], &yy[0], yy.num_of_limbs, &yy[0]);

        if (x.compare(yy) >= 0) {
            x.sub(yy);
        } else {
            add_helper(&x[i - Y.num_of_limbs], &x[i - Y.num_of_limbs], &Y[0], Y.num_of_limbs, x.num_of_limbs - i - 1);
            x.sub(yy);
            q[i - Y.num_of_limbs]--;
        }
    }

    for (auto i = 0; i != r.num_of_limbs; i++) {
        r[i] = x[i];
    }
}

template <unsigned N, unsigned M, typename V = std::enable_if_t<Copyable<N, M>::value>>
static inline mpi<N> operator/(mpi<N> const& X, mpi<M> const& Y)
{
    mpi<N> q;
    mpi<M> r;
    div(X, Y, q, r);
    return q;
}

// r = A mod m
template <unsigned N, unsigned M>
static inline void fp_reduce_by_div(mpi<M> const& A, mpi<N> const& m, mpi<N>& r)
{
    mpi<M> q;
    div(A, m, q, r);
}

template <unsigned N>
struct barrett_helper {
    static constexpr unsigned K = mpi<N>::num_of_limbs;

    barrett_helper(mpi<N> const& m) : m(m)
    {
        mpi<K * 2 * bits_per_limb + 1> b2k(1);
        b2k.lshift(K * 2 * bits_per_limb);
        mpi<N> r;
        div(b2k, m, u, r);
    }

    void reduce(mpi<N * 2> const& x, mpi<N>& r) const
    {
        // q1 = x/b^(K - 1)
        mpi<(K + 1) * bits_per_limb> q1(&x[K - 1], K + 1);

        auto q2 = q1 * u;

        // q3 = q2/b^(K + 1)
        decltype(u) q3(&q2[K + 1], u.num_of_limbs);

        // quick modulus, r1 = x mod b^(k + 1)
        mpi<(K + 1) * bits_per_limb> r1(&x[0], K + 1);

        // quick modulus, r2 = q3 * m mod b^(k + 1)
        auto r2 = resize<(K + 1) * bits_per_limb>(q3 * m);

        if (r1.compare(r2) >= 0) {
            r1.sub(r2);
        } else {
            // r = r1 + b^(K + 1) - r2
            mpi<(K + 2) * bits_per_limb> tmp(r1);
            tmp[K + 1] = 1;
            tmp.sub(r2);
            resize(tmp, r1);
        }

        if (r1.compare(m) >= 0) {
            r1.sub(m);
        }

        resize(r1, r);
    }

    const mpi<N>& m;
    mpi<K * 2 * bits_per_limb + 1> u; // b^2k/m
};

template <unsigned N>
struct standard_modulus_op {
    standard_modulus_op(mpi<N> const& m) : m(m), bh(m){};

    template <unsigned M>
    inline void operator()(mpi<M> const& A, mpi<N>& r) const
    {
        if constexpr (M > N * 2) {
            fp_reduce_by_div(A, m, r);
        } else if constexpr (M == N * 2) {
            bh.reduce(A, r);
        } else if constexpr (M == N) {
            // a couple of addition/subtraction will do for this simple case
            if (m[m.num_of_limbs - 1] & (static_cast<limb_t>(1) << (bits_per_limb - 1))) {
                r = A;
                while (r.compare(m) >= 0) {
                    r.sub(m);
                }
                while (r.sign < 0 && !r.equal(0)) {
                    r.add(m);
                }
                return;
            }

            bh.reduce(resize<N * 2>(A), r);
        } else {
            bh.reduce(resize<N * 2>(A), r);
        }
    }

    mpi<N> const& m;
    barrett_helper<N> bh;
};

// unsafe version of addition without carry
template <unsigned N, unsigned M>
static inline decltype(auto) operator+(mpi<N> const& A, mpi<M> const& B)
{
    if constexpr (N >= M) {
        mpi<N> tmp(A);
        tmp.add(B);
        return tmp;
    } else {
        mpi<M> tmp(B);
        tmp.add(A);
        return tmp;
    }
}

// unsafe version of substract without carry
template <unsigned N, unsigned M>
static inline decltype(auto) operator-(mpi<N> const& A, mpi<M> const& B)
{
    if constexpr (N >= M) {
        mpi<N> tmp(A);
        tmp.sub(B);
        return tmp;
    } else {
        mpi<N> tmp(B);
        tmp.sub(A);
        tmp.negate();
        return tmp;
    }
}

// A must >= B
template <unsigned N>
void __gcd(mpi<N> const& A, mpi<N> const& B, mpi<N>& g)
{
    mpi<N> u{A};
    mpi<N> v{B};
    mpi<N> e(1);

    while (u.is_even() && v.is_even()) {
        u.rshift(1);
        v.rshift(1);
        e.lshift(1);
    }

    while (!u.equal(0)) {
        while (u.is_even()) {
            u.rshift(1);
        }

        while (v.is_even()) {
            v.rshift(1);
        }

        if (u.compare(v) >= 0) {
            u.sub(v);
        } else {
            v.sub(u);
        }
    }

    g = resize<N>(e * v);
}

template <unsigned N>
static inline decltype(auto) gcd(mpi<N> const& A, mpi<N> const& B, mpi<N>& g)
{
    if (A.compare(B) >= 0) {
        return __gcd(A, B, g);
    } else {
        return __gcd(B, A, g);
    }
}

template <unsigned N>
static inline void square(mpi<N> const& A, mpi<N * 2>& C)
{
    mul(A, A, C);
}

template <unsigned N>
class montg_helper
{
public:
    montg_helper(mpi<N> const& m) : m(m)
    {
        limb_t x, m0 = m[0];
        unsigned int i;

        x = m0;
        x += ((m0 + 2) & 4) << 1;

        for (i = bits_per_limb; i >= 8; i /= 2)
            x *= (2 - (m0 * x));

        mm = ~x + 1;

        // R^2 mod m
        mpi<N * 2 + 1> rr(1);
        rr.lshift(N * 2);
        fp_reduce_by_div(rr, m, R2m);
    }

    // r = X * Y * R^-1 mod m
    void montred(mpi<N> const& X, mpi<N> const& Y, mpi<N>& r) const
    {
        mpi<(mpi<N>::num_of_limbs * 2 + 1) * bits_per_limb> T;
        limb_t* A = &T[0];

        for (unsigned i = 0; i < m.num_of_limbs; i++) {
            // compute A = (A + Xi*Y + u*m) / b
            // where u = (a0 + XiY0) * mm mod b
            limb_t u = (A[0] + X[i] * Y[0]) * mm;

            limb_t c = inc_mul_helper(X[i], &Y[0], Y.num_of_limbs, A, T.num_of_limbs - Y.num_of_limbs - i);

            inc_mul_helper(u, &m[0], m.num_of_limbs, A, T.num_of_limbs - m.num_of_limbs - i, c);

            A++; // divided by radix b
            A[m.num_of_limbs + 1] = 0;
        }

        if (compare_helper(A, m.num_of_limbs + 1, &m[0], m.num_of_limbs) >= 0) {
            sub_helper(A, A, &m[0], m.num_of_limbs, 1);
        }

        for (int i = 0; i != r.num_of_limbs; i++) {
            r[i] = *A++;
        }
    }

    inline void montred(mpi<N> const& T, mpi<N>& r) const
    {
        montred(T, mpi<N>(1), r);
    }

protected:
    limb_t mm;
    mpi<N> R2m; // R^2 mod m
    const mpi<N>& m;
};

template <unsigned N>
class montg_inverse_helper : public montg_helper<N>
{
    using base = montg_helper<N>;

public:
    montg_inverse_helper(mpi<N> const& m) : base(m){};

    // compute R = A^-1 * R^2 mod m
    void montinv(mpi<N> const& a, mpi<N>& R)
    {
        mpi<N> r(a);
        mpi<N + 1> s(1);
        mpi<N> t(base::m);
        mpi<N + 1> v(0);
        unsigned k = 0;

        while (!r.equal(0)) {
            if (t.is_even()) {
                t.rshift(1);
                s.lshift(1);
            } else if (r.is_even()) {
                r.rshift(1);
                v.lshift(1);
            } else if (t.compare(r) > 0) {
                t.sub(r);
                t.rshift(1);
                v.add(s);
                s.lshift(1);
            } else {
                r.sub(t);
                r.rshift(1);
                s.add(v);
                v.lshift(1);
            }

            k++;
        }

        if (v.compare(base::m) >= 0) {
            v.sub(base::m);
        }

        // v = p - v
        mpi<N> x = resize<N>(v);
        x.sub(base::m);
        x.negate();

        if (k < N) {
            base::montred(x, base::R2m, x);
            k += N;
        }

        base::montred(x, base::R2m, x);
        mpi<N> tmp(1);
        tmp.lshift(2 * N - k);
        base::montred(x, tmp, R);
    }

    // compute R = A^-1 mod m
    void operator()(mpi<N> const& a, mpi<N>& R)
    {
        mpi<N> r(a);
        mpi<N + 1> s(1);
        mpi<N> t(base::m);
        mpi<N + 1> v(0);
        unsigned k = 0;

        while (!r.equal(0)) {
            if (t.is_even()) {
                t.rshift(1);
                s.lshift(1);
            } else if (r.is_even()) {
                r.rshift(1);
                v.lshift(1);
            } else if (t.compare(r) > 0) {
                t.sub(r);
                t.rshift(1);
                v.add(s);
                s.lshift(1);
            } else {
                r.sub(t);
                r.rshift(1);
                s.add(v);
                v.lshift(1);
            }
            k++;
        }

        if (v.compare(base::m) >= 0) {
            v.sub(base::m);
        }

        // v = p - v
        mpi<N> x = resize<N>(v);
        x.sub(base::m);
        x.negate();

        if (k > N) {
            base::montred(x, x);
            k -= N;
        }

        mpi<N> tmp(1);
        tmp.lshift(N - k);
        base::montred(x, tmp, R);
    }
};

template <unsigned N>
class montg_exp_helper : public montg_helper<N>
{
    using base = montg_helper<N>;

public:
    montg_exp_helper(mpi<N> const& m) : base(m)
    {
        // R^2 * R^-1 mod N = R mod N
        base::montred(base::R2m, Rm);
    }

    template <unsigned M>
    void exp_mod_simple(mpi<N> const& X, mpi<M> const& E, mpi<M>& r) const
    {
        mpi<N> xx;
        base::montred(X, base::R2m, xx);

        mpi<N> A(Rm);
        bit_iterator p(E);

        while (!p.end() && !p) {
            p.next();
        }

        while (!p.end()) {
            base::montred(A, A, A);
            if (p) {
                base::montred(A, xx, A);
            }
            p.next();
        }

        base::montred(A, A);
        r = A;
    }

    // sliding window algorithm
    // R = g^e mod m
    template <unsigned M, unsigned K = 3>
    void operator()(mpi<N> const& g, mpi<M> const& e, mpi<N>& R) const
    {
        // precomputation
        // for a k bit window, the least significant bit is always 1
        // so we save a space of 2^(k - 1) only (with implict 1 at least significant bit)
        // TODO, if we use fixed g for a number of time, it's worthing buffering G somewhere
        mpi<N> G[1 << (K - 1)];
        base::montred(g, base::R2m, G[0]); // G[0] = A * R mod N

        mpi<N> G2;
        base::montred(G[0], G[0], G2);

        for (int i = 1; i < (1 << (K - 1)); i++) {
            base::montred(G[i - 1], G2, G[i]);
        };

        mpi<N> A(Rm);
        bit_iterator p(e);

        while (!p.end() && !p) {
            p.next(); // find the first bit 1
        }

        while (!p.end()) {
            if (!p) {
                base::montred(A, A, A); // squaring
                p.next();
                continue;
            }

            // find the longest bitstring <= K (that start with '1' and end with '1')
            unsigned bitstring = 1;
            unsigned bits = 1;
            do {
                p.next();
                if (p.end()) {
                    p.prev();
                    break;
                }
                bits++;
                bitstring <<= 1;
                if (p)
                    bitstring++;
                if (bits == K)
                    break;
            } while (1);

            while (!p) {
                p.prev();
                bitstring >>= 1;
                bits--;
            };

            while (bits--) {
                base::montred(A, A, A); // squaring
            }

            base::montred(A, G[bitstring >> 1], A); // multiplication

            p.next();
        };

        base::montred(A, A);
        R = A;
    }

private:
    mpi<N> Rm;  // R mod m
};

template <unsigned M, unsigned N>
static inline void fp_reduce(mpi<M> const& A, mpi<N> const& m, mpi<N>& r)
{
    standard_modulus_op mop(m);
    mop(A, r);
}

template <unsigned N, unsigned M>
static inline void fp_exp_mod(mpi<N> const& g, mpi<M> const& e, mpi<N> const& m, mpi<N>& R)
{
    montg_exp_helper meh(m);
    meh(g, e, R);
}

template <unsigned N>
void fp_inverse(mpi<N> const& A, mpi<N> const& p, mpi<N>& r)
{
    montg_inverse_helper iop(p);
    iop(A, r);
}

} // namespace mpi_ns

#endif
