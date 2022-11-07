#!/usr/bin/python3

from math import ceil
from utils import bit_count


def calculate_mu() -> int:
    """
    Calculates Montgomery magic constant for secp256k1 scalar field,
    see algorithm 3 of https://eprint.iacr.org/2017/1057.pdf
    """
    y = 1
    for i in range(2, RADIX_BIT_LEN + 1):
        if (N * y) % (1 << i) != 1:
            y = y + (1 << (i - 1))
    return RADIX - y


# Secp256k1 scalar field prime
# see section 2.4.1 of https://www.secg.org/sec2-v2.pdf#page=13
N = 2**256 - 432420386565659656852420866394968145599

# Secp256k1 curve generator point x, see section 2.4.1 of https://www.secg.org/sec2-v2.pdf#page=13
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
# Secp256k1 curve generator point y, , see section 2.4.1 of https://www.secg.org/sec2-v2.pdf#page=13
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424

RADIX_BIT_LEN: int = 32
RADIX: int = 1 << RADIX_BIT_LEN

LIMB_COUNT: int = ceil(bit_count(N) / RADIX_BIT_LEN)

# = (2 ^ 32) ^ 8 = 2 ^ 256 % n
R: int = (RADIX**LIMB_COUNT) % N
# = (2 ^ 256) ^ 2 % n
R2: int = (R * R) % N
MU: int = calculate_mu()

if __name__ == "__main__":
    print("Use `scalar_field_consts` as library module")
