#!/usr/bin/python3

from math import ceil


def bit_count(num: int) -> int:
    '''
    Computes bitlength of integer, same as len(bin(num)[2:])
    '''
    cnt = 0
    num_ = num

    while(num > 0):
        cnt += 1
        num >>= 1

    assert cnt == len(bin(num_)[2:])
    return cnt


def calculate_mu() -> int:
    '''
    Calculates Montgomery magic constant for secp256k1 prime field, 
    see algorithm 3 of https://eprint.iacr.org/2017/1057.pdf
    '''
    y = 1
    for i in range(2, RADIX_BIT_LEN + 1):
        if (PRIME * y) % (1 << i) != 1:
            y = y + (1 << (i - 1))
    return RADIX - y


# = p; See https://en.bitcoin.it/wiki/Secp256k1
PRIME: int = 0x_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFC2F
PRIME_BIT_LEN: int = bit_count(PRIME)
# can be rewritten using uint32_t data type of C
RADIX_BIT_LEN: int = 32
RADIX: int = 1 << RADIX_BIT_LEN
# = 8
LIMB_COUNT: int = ceil(PRIME_BIT_LEN / RADIX_BIT_LEN)
# = (2 ^ 32) ^ 8 = 2 ^ 256 % p
R: int = (RADIX ** LIMB_COUNT) % PRIME
# = (2 ^ 256) ^ 2 % p
R2: int = (R * R) % PRIME
MU: int = calculate_mu()

TEST_CNT: int = 1 << 10
