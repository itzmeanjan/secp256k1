#!/usr/bin/python3

from math import ceil
from . import bit_count


def calculate_mu() -> int:
    """
    Calculates Montgomery magic constant for secp256k1 base field,
    see algorithm 3 of https://eprint.iacr.org/2017/1057.pdf
    """
    y = 1
    for i in range(2, RADIX_BIT_LEN + 1):
        if (P * y) % (1 << i) != 1:
            y = y + (1 << (i - 1))
    return RADIX - y


# Secp256k1 base field prime
# see section 2.4.1 of https://www.secg.org/sec2-v2.pdf#page=13
P: int = 0x_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFC2F

RADIX_BIT_LEN: int = 32
RADIX: int = 1 << RADIX_BIT_LEN

LIMB_COUNT: int = ceil(bit_count(P) / RADIX_BIT_LEN)

# = (2 ^ 32) ^ 8 = 2 ^ 256 % p
R: int = (RADIX**LIMB_COUNT) % P
# = (2 ^ 256) ^ 2 % p
R2: int = (R * R) % P
MU: int = calculate_mu()

if __name__ == "__main__":
    print("Use `base_field_consts` as library module")
