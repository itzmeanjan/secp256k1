#!/usr/bin/python3

from typing import Tuple
from field import Gx, Gy, N
from field import BaseField
from field import ScalarField
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
    h = h % N

    s1 = ScalarField.from_num(s).inv()

    t0 = ScalarField.from_num(h)
    t1 = ScalarField.from_num(r)

    t2 = t0 * s1
    t3 = t1 * s1

    t4 = t2.to_num()
    t5 = t3.to_num()

    g = Point.fromAffine(BaseField.from_num(Gx), BaseField.from_num(Gy))

    t6 = g.mulScalar(t4)
    t7 = pkey.mulScalar(t5)

    t8 = t6 + t7
    t9 = t8.toAffine()[0].to_num()

    return r == t9


if __name__ == "__main__":
    print("Use `verify` as library module")
