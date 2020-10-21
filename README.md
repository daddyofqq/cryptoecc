# Introduction

The implementation is in two header files only, mpi.h and ecc.h, which you can easily include in your C++ project. No linking of library code is needed. (but currently it does have some dependency on external source code, see Porting Guide)

Three types of ECC curve supported
* Short Weierstrass : y^2 = x^3 + ax + b mod p
* Montgomery curve : y^2 = x^3 + A x^2 + x mod p
* Twisted Edward Curve : a * x^2 + y^2 = 1 + d x^2 y^2

Algorithm for ECDSA, EDDSA, and ECDH are implemented.

Support for Well-defined curve parameters:

* NIST p256
* curve25519
* curve448
* ed25519
* ed448

Other parameters can be easily supported if you bother to look at the domain parameter template in the ecc.h.

You can figure out how to use these stuff mostly by looking at the demo.cpp of project.


# Compile code and how to use it

When you pull code, make sure you also pull git submodules (git submodule update --init --recursive)

run 'make' in root folder and test folder will compile demo code and test code respectively.

demo code shows basics of signature signing and ECDH for all supported curves. test code has more intense stuff, including some performance benchmark. It shows some interesting data regarding mbedtls and this implementation.

# Porting Guide

The implementation does not have much dependency on target system except that you will have to create a decent random number generator (used by implementation to generate secret keys and for counter-measure to side-channel attack). For a solution that needs to get PCI level security, the default rng implementation in default_rng.h is not good enough. You are on your own to go extra miles.

On the other hand, EDDSA needs sha512 and shake256 for signature scheme. It is not the goal of this project to implement such hashing alg from scratch. So you are welcome to use the reference implementation wrapped on mbedtls and tiny_sha3 project (in test/hash.cpp), or you can DIY.

