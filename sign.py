#!/usr/bin/python3


from typing import Tuple
from hashlib import sha3_256
from field import Field
from point import Point
from secrets import randbelow
from consts import Gx, Gy, n


def sign(skey: int, msg: bytes) -> Tuple[int, int]:
    """
    Given ECDSA secret key ( a 256 -bit integer ) and a message `m`, this routine attempts to
    perform randomized signing, while hashing the message using SHA3-256.

    Returns (r, s) two 256 -bit integers ( ∈ [0, n) ), as ECDSA signature.

    Follows scheme described https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-messages#ecdsa-sign
    """
    h = sha3_256(msg).digest()
    h = int.from_bytes(h, byteorder="big")
    h = h % n

    k = 1 + randbelow(n - 1)

    g = Point.fromAffine(Field.from_num(Gx), Field.from_num(Gy))
    r = g.mulScalar(k)
    r = r.toAffine()[0].to_num()

    k_inv = pow(k, -1, n)
    s = ((h + r * skey) * k_inv) % n

    return r, s


if __name__ == "__main__":
    print("Use `sign` as library module")
