// SPDX-License-Identifier: GPL-3.0-only
/*
 *  hash.h
 *
 *  Copyright (C) 2019 Daniel Hu <daddy.of.qq@gmail.com>
 */
#ifndef HASH_H
#define HASH_H

// This header file provides the necessary interface for sha512 / shake256
// to be used by ECC crypto(to be more specific, EDDSA)
// a reference implementation based on open source project mbedtls / tiny_sha3
// is provided in test folder
// what users want to choose for their own implementation is out of scope
//
#include <cstddef>

struct hash_op {
    virtual size_t hashlen() const = 0;
    virtual size_t blocksize() const = 0;
    virtual void init() = 0;
    virtual void update(const unsigned char* input, size_t ilen) = 0;
    virtual void finish(unsigned char* output) = 0;
    virtual ~hash_op(){};
};

class default_shake256 : public hash_op
{
private:
    void* pimp;

public:
    static inline constexpr unsigned int block_size = 1088 / 8;
    default_shake256(size_t hashlen);
    virtual size_t hashlen() const override
    {
        return hash_size;
    }
    virtual size_t blocksize() const override
    {
        return block_size;
    }
    virtual void init() override;
    virtual void update(const unsigned char* input, size_t ilen) override;
    virtual void finish(unsigned char* output) override;
    ~default_shake256();
    const size_t hash_size;
};

class default_sha512 : public hash_op
{
private:
    void* pimp;

public:
    static inline constexpr unsigned int hash_size = 64;
    static inline constexpr unsigned int block_size = 1024 / 8;
    default_sha512();
    ~default_sha512();
    virtual void init() override;
    virtual size_t hashlen() const override
    {
        return hash_size;
    }
    virtual size_t blocksize() const override
    {
        return block_size;
    }
    virtual void update(const unsigned char* input, size_t ilen) override;
    virtual void finish(unsigned char* output) override;
};

class default_sha384 : public hash_op
{
private:
    void* pimp;

public:
    static inline constexpr unsigned int hash_size = 48;
    static inline constexpr unsigned int block_size = 1024 / 8;

    default_sha384();
    ~default_sha384();
    virtual void init() override;
    virtual size_t hashlen() const override
    {
        return hash_size;
    }
    virtual size_t blocksize() const override
    {
        return block_size;
    }
    virtual void update(const unsigned char* input, size_t ilen) override;
    virtual void finish(unsigned char* output) override;
};

class default_sha256 : public hash_op
{
private:
    void* pimp;

public:
    static inline constexpr unsigned int hash_size = 32;
    static inline constexpr unsigned int block_size = 512 / 8;

    default_sha256();
    ~default_sha256();
    virtual void init() override;
    virtual size_t hashlen() const override
    {
        return hash_size;
    }
    virtual size_t blocksize() const override
    {
        return block_size;
    }
    virtual void update(const unsigned char* input, size_t ilen) override;
    virtual void finish(unsigned char* output) override;
};

#endif
