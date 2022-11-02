#!/usr/bin/python3

from typing import Tuple
from consts import PRIME, Gx, Gy
from point import Point
from field import Field
from secrets import randbelow


def generate_secret_key() -> int:
    """
    Generate a random ECDSA secret key, which is ∈ [0, p) | p = PRIME
    """
    skey = randbelow(PRIME)
    return skey


def generate_public_key(skey: int) -> Point:
    """
    Given an ECDSA secret key, this routine generates corresponding ECDSA public key
    s.t. pkey = skey * G | G = secp256k1 generator point
    """
    gen = Point.fromAffine(Field.from_num(Gx), Field.from_num(Gy))
    pkey = gen.mulScalar(skey)

    return pkey


def keygen() -> Tuple[int, Point]:
    """
    Generate a random ECDSA secret, public key pair ( in order )

    Follows scheme described https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-messages#key-generation
    """
    skey = generate_secret_key()
    pkey = generate_public_key(skey)

    return skey, pkey


if __name__ == "__main__":
    print("Use `keygen` as library module")
