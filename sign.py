#!/usr/bin/python3


from typing import Tuple
from hashlib import sha3_256
from field import Field
from point import Point
from secrets import randbelow
from consts import PRIME, Gx, Gy


def sign(skey: int, msg: bytes) -> Tuple[int, int]:
    """
    Given ECDSA secret key ( a 256 -bit integer ) and a message `m`, this routine attempts to
    perform randomized signing, while hashing the message using SHA3-256.

    Returns (r, s) two 256 -bit integers ( ∈ [0, PRIME) ), as ECDSA signature.

    Follows scheme described https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-messages#ecdsa-sign
    """
    h = sha3_256(msg).digest()
    h = int.from_bytes(h, byteorder="big")
    h = h % PRIME

    k = 1 + randbelow(PRIME - 1)

    g = Point.fromAffine(Field.from_num(Gx), Field.from_num(Gy))
    r = g.mulScalar(k)
    r = r.toAffine()[0]

    t0 = Field.from_num(skey)
    t2 = Field.from_num(h)
    t3 = Field.from_num(k)

    t4 = t3.inv() * (t2 + r * t0)
    s = t4.to_num()

    return r.to_num(), s


if __name__ == "__main__":
    print("Use `sign` as library module")
