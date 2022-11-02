#!/usr/bin/python3

from typing import Tuple
from consts import PRIME, Gx, Gy
from field import Field
from point import Point
from hashlib import sha3_256


def verify(pkey: Point, msg: bytes, sig: Tuple[int, int]) -> bool:
    """
    Given ECDSA public key, message `m` and signature tuple ( i.e. (r, s) ), this routine
    attempts to verify signature.

    Returns boolean value denoting success.

    Follows scheme described https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-messages#ecdsa-verify-signature
    """
    (r, s) = sig

    h = sha3_256(msg).digest()
    h = int.from_bytes(h, byteorder="big")
    h = h % PRIME

    t0 = Field.from_num(h)
    t1 = Field.from_num(r)

    t2 = Field.from_num(s)
    t3 = t2.inv()

    t4 = t0 * t3  # = h * s1
    t5 = t1 * t3  # = r * s1

    g = Point.fromAffine(Field.from_num(Gx), Field.from_num(Gy))
    r_ = g.mulScalar(t4.to_num()) + pkey.mulScalar(t5.to_num())
    r_ = r_.toAffine()[0].to_num()

    return r == r_


if __name__ == "__main__":
    print("Use `verify` as library module")
