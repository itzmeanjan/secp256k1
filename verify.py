#!/usr/bin/python3

from typing import Tuple
from scalar_field_consts import Gx, Gy, N
from base_field import BaseField
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

    s1 = pow(s, -1, N)

    t0 = (h * s1) % N
    t1 = (r * s1) % N

    g = Point.fromAffine(BaseField.from_num(Gx), BaseField.from_num(Gy))

    t2 = g.mulScalar(t0)
    t3 = pkey.mulScalar(t1)

    t4 = t2 + t3
    t5 = t4.toAffine()[0].to_num()

    return r == t5


if __name__ == "__main__":
    print("Use `verify` as library module")
