#!/usr/bin/python3

from typing import Tuple
from field import Gx, Gy, N
from point import Point
from field import BaseField
from secrets import randbelow


def generate_secret_key() -> int:
    """
    Generate a random ECDSA secret key, which is ∈ [0, N)
    """
    skey = randbelow(N)
    return skey


def generate_public_key(skey: int) -> Point:
    """
    Given an ECDSA secret key, this routine generates corresponding ECDSA public key
    s.t. pkey = skey * G | G = secp256k1 generator point
    """
    gen = Point.fromAffine(BaseField.from_num(Gx), BaseField.from_num(Gy))
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
