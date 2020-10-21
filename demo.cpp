// SPDX-License-Identifier: GPL-3.0-only
/*
 *  demo.cpp
 *
 *  Copyright (C) 2020 Daniel Hu <daddy.of.qq@gmail.com>
 *
 */

#include "default_rng.h"
#include "ecc.h"
#include <iostream>

using namespace ecc_ns;

template <typename C>
void print_statis(C const& curve)
{
    auto statis = curve.get_statis();
    std::cout << std::string(C::dp_template::curve_name) << " scalar mul oprs ";
    std::cout << std::string(" [INV ") << statis.I;
    std::cout << std::string(", MUL ") << statis.M;
    std::cout << std::string(", ADD ") << statis.A;
    std::cout << std::string(", FADD ") << statis.FA;
    std::cout << std::string(", RED ") << statis.R;
    std::cout << std::string(", SHIFT ") << statis.F;
    std::cout << std::string("]") << std::endl;
};

default_rng rng;
ed448 ed448_inst(&rng);
ed25519 ed25519_inst(&rng);
curve25519 curve25519_inst(&rng);
secp256r1 secp256r1_inst(&rng);

static bool gen_ed448_key(std::vector<uint8_t>& prv, std::vector<uint8_t>& pub)
{
    ed448::keypair kp;
    auto ret = ed448_inst.genkeypair(kp);
    if (ret) {
        // ed448 output keypair includes the private key as mpi<>
        // and public key as encoded point, that is, FixedSizeBuf<>
        // you can transform them into a format easy to save and transfer
        // like here in std::vector
        auto d = kp.d.output(); // mpi<> into FixedSizeBuf<>
        prv = std::vector<uint8_t>(d.begin(), d.end());
        pub = std::vector<uint8_t>(kp.Q.begin(), kp.Q.end());
    }
    print_statis(ed448_inst);
    return ret;
}

static void demo_ed448()
{
    std::vector<uint8_t> prv1, pub1;
    std::string msg = "see you on friday";
    std::string ctx = "on the beach";

    auto ret = gen_ed448_key(prv1, pub1);
    if (!ret) {
        std::cout << std::string("fail to generate ed448 key\n");
        exit(1);
    }

    ed448::signature_t sig;
    ed448::prvkey_t d(prv1.data(), prv1.size()); // restore kp.d

    // ed448 context is optional. if empty, pass nullptr or omit from argument list
    ret = ed448_inst.eddsa_sign(d,
                                reinterpret_cast<const uint8_t*>(msg.data()), msg.size(),
                                sig,
                                reinterpret_cast<const uint8_t*>(ctx.data()), ctx.size());
    if (!ret) {
        std::cout << std::string("fail to sign ed448 msg\n");
        exit(1);
    }

    ed448::pubkey_t pk(pub1.data()); // restore kp.Q
    ret = ed448_inst.eddsa_verify(pk,
                                  reinterpret_cast<const uint8_t*>(msg.data()), msg.size(),
                                  sig,
                                  reinterpret_cast<const uint8_t*>(ctx.data()), ctx.size());

    if (!ret) {
        std::cout << std::string("fail to verify ed448 signature\n");
        exit(1);
    }

    std::cout << std::string("ed448 signature test completed successfully\n");

    // generate second keypair
    std::vector<uint8_t> prv2, pub2;
    ret = gen_ed448_key(prv2, pub2);
    if (!ret) {
        std::cout << std::string("fail to generate ed448 key\n");
        exit(1);
    }

    // generate secret from private key 2 and public key 1
    ed448::shared_secret_t secret1;
    ret = ed448_inst.ecdh_compute_shared(ed448::prvkey_t(prv2.data(), prv2.size()),
                                         ed448::pubkey_t(pub1.data()),
                                         secret1);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // generate secret from private key 1 and public key 2
    ed448::shared_secret_t secret2;
    ret = ed448_inst.ecdh_compute_shared(ed448::prvkey_t(prv1.data(), prv1.size()),
                                         ed448::pubkey_t(pub2.data()),
                                         secret2);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // of course the two secrets should be exactly the same
    // that is the whole point of EDCH
    if (secret1.compare(secret2) != 0) {
        std::cout << std::string("fail to derive the same secret key\n");
        exit(1);
    }

    std::cout << std::string("ed448 edch tested successfully\n");
}

static bool gen_ed25519_key(std::vector<uint8_t>& prv, std::vector<uint8_t>& pub)
{
    ed25519::keypair kp;
    auto ret = ed25519_inst.genkeypair(kp);
    if (ret) {
        // ed25519 output keypair includes the private key as mpi<>
        // and public key as encoded point, that is, FixedSizeBuf<>
        // you can transform them into a format easy to save and transfer
        // like here in std::vector
        auto d = kp.d.output();
        prv = std::vector<uint8_t>(d.begin(), d.end());
        pub = std::vector<uint8_t>(kp.Q.begin(), kp.Q.end());
    }
    print_statis(ed25519_inst);
    return ret;
}

static void demo_ed25519()
{
    std::vector<uint8_t> prv1, pub1;
    std::string msg = "see you on friday";

    auto ret = gen_ed25519_key(prv1, pub1);
    if (!ret) {
        std::cout << std::string("fail to generate ed25519 key\n");
        exit(1);
    }

    ed25519::signature_t sig;
    ed25519::prvkey_t d(prv1.data(), prv1.size()); // restore kp.d

    // Note: ed25519 do not use ctx for signature at all. It will always be ignored
    ret = ed25519_inst.eddsa_sign(d,
                                  reinterpret_cast<const uint8_t*>(msg.data()), msg.size(),
                                  sig);
    if (!ret) {
        std::cout << std::string("fail to sign ed25519 msg\n");
        exit(1);
    }

    ed25519::pubkey_t pk(pub1.data());
    ret = ed25519_inst.eddsa_verify(pk,
                                    reinterpret_cast<const uint8_t*>(msg.data()), msg.size(),
                                    sig);

    if (!ret) {
        std::cout << std::string("fail to verify ed25519 signature\n");
        exit(1);
    }

    std::cout << std::string("ed25519 signature test completed successfully\n");

    // generate second keypair
    std::vector<uint8_t> prv2, pub2;
    ret = gen_ed25519_key(prv2, pub2);
    if (!ret) {
        std::cout << std::string("fail to generate ed25519 key\n");
        exit(1);
    }

    // generate secret from private key 2 and public key 1
    ed25519::shared_secret_t secret1;
    ret = ed25519_inst.ecdh_compute_shared(ed25519::prvkey_t(prv2.data(), prv2.size()),
                                           ed25519::pubkey_t(pub1.data()),
                                           secret1);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // generate secret from private key 1 and public key 2
    ed25519::shared_secret_t secret2;
    ret = ed25519_inst.ecdh_compute_shared(ed25519::prvkey_t(prv1.data(), prv1.size()),
                                           ed25519::pubkey_t(pub2.data()),
                                           secret2);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // of course the two secrets should be exactly the same
    // that is the whole point of EDCH
    if (secret1.compare(secret2) != 0) {
        std::cout << std::string("fail to derive the same secret key\n");
        exit(1);
    }

    std::cout << std::string("ed25519 edch tested successfully\n");

    // this is bonus feature for someone interested
    // you are probably already aware that ed25519 and curve25519 are equivalent
    // there is a birational mapping between them
    // following code demonstrate this
    ecp_aff<256> point_on_ed25519;
    ret = ed25519_inst.decoding(ed25519::pubkey_t(pub2.data()), point_on_ed25519);
    if (!ret) {
        std::cout << std::string("fail to decode public point\n");
        exit(1);
    }
    auto point_on_curve25519 = ed25519_inst.to_montgomery(point_on_ed25519);
    if (curve25519_inst.verify_point(point_on_curve25519)) {
        std::cout << std::string("Yes, ed25519/curve25519 are equivalent!\n");
    } else {
        std::cout << std::string("point validation failed\n");
        exit(1);
    }
}

static bool gen_curve25519_key(std::vector<uint8_t>& prv, std::vector<uint8_t>& pub)
{
    curve25519::keypair kp;
    auto ret = curve25519_inst.genkeypair(kp);
    if (ret) {
        // curve25519 output keypair includes the private key as mpi<>
        // and public key as mpi<>, which is simply the X coordinates
        // you can transform them into a format easy to save and transfer
        // like here in std::vector
        auto d = kp.d.output();
        auto Q = kp.Q.output();
        prv = std::vector<uint8_t>(d.begin(), d.end());
        pub = std::vector<uint8_t>(Q.begin(), Q.end());
    }
    print_statis(curve25519_inst);
    return ret;
}

static void demo_curve25519()
{
    std::vector<uint8_t> prv1, pub1;

    auto ret = gen_curve25519_key(prv1, pub1);
    if (!ret) {
        std::cout << std::string("fail to generate curve25519 key\n");
        exit(1);
    }

    // generate second keypair
    std::vector<uint8_t> prv2, pub2;
    ret = gen_curve25519_key(prv2, pub2);
    if (!ret) {
        std::cout << std::string("fail to generate curve25519 key\n");
        exit(1);
    }

    // generate secret from private key 2 and public key 1
    curve25519::shared_secret_t secret1;
    ret = curve25519_inst.ecdh_compute_shared(curve25519::prvkey_t(prv2.data(), prv2.size()),
                                              curve25519::pubkey_t(pub1.data(), pub1.size()),
                                              secret1);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // generate secret from private key 1 and public key 2
    curve25519::shared_secret_t secret2;
    ret = curve25519_inst.ecdh_compute_shared(curve25519::prvkey_t(prv1.data(), prv1.size()),
                                              curve25519::pubkey_t(pub2.data(), pub2.size()),
                                              secret2);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // of course the two secrets should be exactly the same
    // that is the whole point of EDCH
    if (secret1.compare(secret2) != 0) {
        std::cout << std::string("fail to derive the same secret key\n");
        exit(1);
    }

    std::cout << std::string("curve25519 edch tested successfully\n");
}

static bool gen_secp256r1_key(std::vector<uint8_t>& prv,
                              std::vector<uint8_t>& pub_x,
                              std::vector<uint8_t>& pub_y)
{
    secp256r1::keypair kp;
    auto ret = secp256r1_inst.genkeypair(kp);
    if (ret) {
        // secp256r1 output keypair includes the private key as mpi<>
        // and public key as ecp_aff<>, which is affine coordinates of point Q
        // you can transform them into a format easy to save and transfer
        // like here in std::vector
        auto d = kp.d.output();
        auto Qx = kp.Q.X.output();
        auto Qy = kp.Q.Y.output();
        prv = std::vector<uint8_t>(d.begin(), d.end());
        pub_x = std::vector<uint8_t>(Qx.begin(), Qx.end());
        pub_y = std::vector<uint8_t>(Qy.begin(), Qy.end());
    }

    print_statis(secp256r1_inst);

    return ret;
}

static void demo_secp256r1()
{
    std::vector<uint8_t> prv1, pub1_x, pub1_y;
    auto ret = gen_secp256r1_key(prv1, pub1_x, pub1_y);
    if (!ret) {
        std::cout << std::string("fail to generate secp256r1 key\n");
        exit(1);
    }

    // typically msg contain hash
    std::vector<uint8_t> msg(32);
    rng(msg.data(), msg.size());

    secp256r1::signature_t sig;
    ret = secp256r1_inst.ecdsa_sign(secp256r1::prvkey_t(prv1.data(), prv1.size()),
                                    mpi<256>(msg.data(), msg.size()),
                                    sig);
    if (!ret) {
        std::cout << std::string("fail to sign message by secp256r1\n");
        exit(1);
    }
    // secp256r1::signature_t contains r, s components as mpi<>
    // if you want to export r, s, you can do
    // auto r = sig.r.output()
    // auto s = sig.s.output()
    // std::vector<uint8_t> rv = std::vector<uint8_t>(r.begin(), r.end())
    // std::vector<uint8_t> sv = std::vector<uint8_t>(s.begin(), send())
    //
    // secp256r1::signature_t can be rebuilt by
    // secp256r1::signature_t sig(rv.data(), rv.size(), sv.data(), sv.size())
    ret = secp256r1_inst.ecdsa_verify(secp256r1::pubkey_t(pub1_x.data(), pub1_x.size(),
                                                          pub1_y.data(), pub1_y.size()),
                                      mpi<256>(msg.data(), msg.size()),
                                      sig);
    if (!ret) {
        std::cout << std::string("fail to verify signature by secp256r1\n");
        exit(1);
    }

    std::cout << std::string("ecdsa by secp256r1 successful\n");

    // generate second keypair
    std::vector<uint8_t> prv2, pub2_x, pub2_y;
    ret = gen_secp256r1_key(prv2, pub2_x, pub2_y);
    if (!ret) {
        std::cout << std::string("fail to generate secp256r1 key\n");
        exit(1);
    }

    // generate secret from private key 2 and public key 1
    secp256r1::shared_secret_t secret1;
    ret = secp256r1_inst.ecdh_compute_shared(secp256r1::prvkey_t(prv2.data(), prv2.size()),
                                             secp256r1::pubkey_t(pub1_x.data(), pub1_x.size(),
                                                                 pub1_y.data(), pub1_y.size()),
                                             secret1);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // generate secret from private key 1 and public key 2
    secp256r1::shared_secret_t secret2;
    ret = secp256r1_inst.ecdh_compute_shared(secp256r1::prvkey_t(prv1.data(), prv1.size()),
                                             secp256r1::pubkey_t(pub2_x.data(), pub2_x.size(),
                                                                 pub2_y.data(), pub2_y.size()),
                                             secret2);
    if (!ret) {
        std::cout << std::string("fail to compute shared secret\n");
        exit(1);
    }

    // of course the two secrets should be exactly the same
    // that is the whole point of EDCH
    if (secret1.compare(secret2) != 0) {
        std::cout << std::string("fail to derive the same secret key\n");
        exit(1);
    }

    std::cout << std::string("secp256r1 edch tested successfully\n");
}

int main(int argc, char* argv[])
{
    demo_ed448();

    demo_ed25519();

    demo_curve25519();

    demo_secp256r1();
}
