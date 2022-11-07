#!/usr/bin/python3


from typing import Tuple
from hashlib import sha3_256
from base_field import BaseField
from point import Point
from secrets import randbelow
from scalar_field_consts import Gx, Gy, N


def sign(skey: int, msg: bytes) -> Tuple[int, int]:
    """
    Given ECDSA secret key ( a 256 -bit integer ) and a message `m`, this routine attempts to
    perform randomized signing, while hashing the message using SHA3-256.

    Returns (r, s) two 256 -bit integers ( ∈ [0, n) ), as ECDSA signature.

    Follows scheme described https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-messages#ecdsa-sign
    """
    h = sha3_256(msg).digest()
    h = int.from_bytes(h, byteorder="big")
    h = h % N

    k = 1 + randbelow(N - 1)

    g = Point.fromAffine(BaseField.from_num(Gx), BaseField.from_num(Gy))
    r = g.mulScalar(k)
    r = r.toAffine()[0].to_num()

    k_inv = pow(k, -1, N)
    s = ((h + r * skey) * k_inv) % N

    return r, s


if __name__ == "__main__":
    print("Use `sign` as library module")
