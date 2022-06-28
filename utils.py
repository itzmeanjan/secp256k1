#!/usr/bin/python3

from consts import *
from typing import List, Tuple


def to_radix_r(num: int) -> List[int]:
    '''
    Converts large integer ( 256 -bit ) to radix-r interleaved representation | r = 2^32
    '''
    limbs = [0] * LIMB_COUNT
    idx = 0
    while num > 0:
        limbs[idx] = num % RADIX
        num //= RADIX
        idx += 1
    return limbs


def from_radix_r(limbs: List[int]) -> int:
    '''
    Converts radix-r interleaved representation, to large integer ( 256 -bit ) | r = 2^32
    '''
    cnt = len(limbs)
    num = 0
    for idx in range(cnt-1, -1, -1):
        num = num * RADIX + limbs[idx]
    return num


def adc(a: int, b: int, carry: int) -> Tuple[int, int]:
    '''
    See https://github.com/dusk-network/bls12_381/blob/ed4d87c6756c0020629edb5d8912a41e338ac85a/src/util.rs#L1-L6
    '''
    tmp = a + b + carry
    return tmp & (RADIX-1), tmp >> RADIX_BIT_LEN


def mac(a: int, b: int, c: int, carry: int) -> Tuple[int, int]:
    '''
    See https://github.com/dusk-network/bls12_381/blob/ed4d87c6756c0020629edb5d8912a41e338ac85a/src/util.rs#L15-L20
    '''
    tmp = a + (b * c) + carry
    return tmp & (RADIX-1), tmp >> RADIX_BIT_LEN


def sbb(a: int, b: int, borrow: int) -> Tuple[int, int]:
    '''
    See https://github.com/dusk-network/bls12_381/blob/2c679a284c008475b543a67ee2300ee58ffe5d11/src/util.rs#L8-L13
    '''
    wrap_at = (RADIX << RADIX_BIT_LEN) - 1
    tmp = (a - (b + (borrow >> (RADIX_BIT_LEN-1)))) % wrap_at
    return tmp & (RADIX-1), tmp >> RADIX_BIT_LEN


def bitwise_not(a: int) -> int:
    '''
    Same as `!a` in C
    '''
    return RADIX - 1 - a


def u256xu32(a: List[int], b: int, c: List[int]) -> List[int]:
    '''
    Collects inspiration from https://github.com/dusk-network/bls12_381/blob/ed4d87c6756c0020629edb5d8912a41e338ac85a/src/fp.rs#L517-L522
    '''
    assert len(a) == (LIMB_COUNT + 1) and len(c) == LIMB_COUNT

    a[0], carry = mac(a[0], b, c[0], 0)
    a[1], carry = mac(a[1], b, c[1], carry)
    a[2], carry = mac(a[2], b, c[2], carry)
    a[3], carry = mac(a[3], b, c[3], carry)
    a[4], carry = mac(a[4], b, c[4], carry)
    a[5], carry = mac(a[5], b, c[5], carry)
    a[6], carry = mac(a[6], b, c[6], carry)
    a[7], a[8] = mac(a[7], b, c[7], carry)

    return a


def montgomery_mul(a: List[int], b: List[int]) -> List[int]:
    '''
    Inspired by https://github.com/dusk-network/bls12_381/blob/ed4d87c6756c0020629edb5d8912a41e338ac85a/src/fp.rs#L437-L560
    and algorithm 2 of https://eprint.iacr.org/2017/1057.pdf
    '''
    assert len(a) == len(b)

    prime = to_radix_r(PRIME)
    cnt = len(a)
    c = [0] * (cnt << 1)

    c[0:9] = u256xu32(c[0:9], a[0], b)
    q = (MU * c[0]) % RADIX

    _, carry = mac(c[0], q, prime[0], 0)
    c[1], carry = mac(c[1], q, prime[1], carry)
    c[2], carry = mac(c[2], q, prime[2], carry)
    c[3], carry = mac(c[3], q, prime[3], carry)
    c[4], carry = mac(c[4], q, prime[4], carry)
    c[5], carry = mac(c[5], q, prime[5], carry)
    c[6], carry = mac(c[6], q, prime[6], carry)
    c[7], carry = mac(c[7], q, prime[7], carry)
    c[8], pc = adc(c[8], 0, carry)

    c[1:10] = u256xu32(c[1:10], a[1], b)
    q = (MU * c[1]) % RADIX

    _, carry = mac(c[1], q, prime[0], 0)
    c[2], carry = mac(c[2], q, prime[1], carry)
    c[3], carry = mac(c[3], q, prime[2], carry)
    c[4], carry = mac(c[4], q, prime[3], carry)
    c[5], carry = mac(c[5], q, prime[4], carry)
    c[6], carry = mac(c[6], q, prime[5], carry)
    c[7], carry = mac(c[7], q, prime[6], carry)
    c[8], carry = mac(c[8], q, prime[7], carry)
    c[9], pc = adc(c[9], pc, carry)

    c[2:11] = u256xu32(c[2:11], a[2], b)
    q = (MU * c[2]) % RADIX

    _, carry = mac(c[2], q, prime[0], 0)
    c[3], carry = mac(c[3], q, prime[1], carry)
    c[4], carry = mac(c[4], q, prime[2], carry)
    c[5], carry = mac(c[5], q, prime[3], carry)
    c[6], carry = mac(c[6], q, prime[4], carry)
    c[7], carry = mac(c[7], q, prime[5], carry)
    c[8], carry = mac(c[8], q, prime[6], carry)
    c[9], carry = mac(c[9], q, prime[7], carry)
    c[10], pc = adc(c[10], pc, carry)

    c[3:12] = u256xu32(c[3:12], a[3], b)
    q = (MU * c[3]) % RADIX

    _, carry = mac(c[3], q, prime[0], 0)
    c[4], carry = mac(c[4], q, prime[1], carry)
    c[5], carry = mac(c[5], q, prime[2], carry)
    c[6], carry = mac(c[6], q, prime[3], carry)
    c[7], carry = mac(c[7], q, prime[4], carry)
    c[8], carry = mac(c[8], q, prime[5], carry)
    c[9], carry = mac(c[9], q, prime[6], carry)
    c[10], carry = mac(c[10], q, prime[7], carry)
    c[11], pc = adc(c[11], pc, carry)

    c[4:13] = u256xu32(c[4:13], a[4], b)
    q = (MU * c[4]) % RADIX

    _, carry = mac(c[4], q, prime[0], 0)
    c[5], carry = mac(c[5], q, prime[1], carry)
    c[6], carry = mac(c[6], q, prime[2], carry)
    c[7], carry = mac(c[7], q, prime[3], carry)
    c[8], carry = mac(c[8], q, prime[4], carry)
    c[9], carry = mac(c[9], q, prime[5], carry)
    c[10], carry = mac(c[10], q, prime[6], carry)
    c[11], carry = mac(c[11], q, prime[7], carry)
    c[12], pc = adc(c[12], pc, carry)

    c[5:14] = u256xu32(c[5:14], a[5], b)
    q = (MU * c[5]) % RADIX

    _, carry = mac(c[5], q, prime[0], 0)
    c[6], carry = mac(c[6], q, prime[1], carry)
    c[7], carry = mac(c[7], q, prime[2], carry)
    c[8], carry = mac(c[8], q, prime[3], carry)
    c[9], carry = mac(c[9], q, prime[4], carry)
    c[10], carry = mac(c[10], q, prime[5], carry)
    c[11], carry = mac(c[11], q, prime[6], carry)
    c[12], carry = mac(c[12], q, prime[7], carry)
    c[13], pc = adc(c[13], pc, carry)

    c[6:15] = u256xu32(c[6:15], a[6], b)
    q = (MU * c[6]) % RADIX

    _, carry = mac(c[6], q, prime[0], 0)
    c[7], carry = mac(c[7], q, prime[1], carry)
    c[8], carry = mac(c[8], q, prime[2], carry)
    c[9], carry = mac(c[9], q, prime[3], carry)
    c[10], carry = mac(c[10], q, prime[4], carry)
    c[11], carry = mac(c[11], q, prime[5], carry)
    c[12], carry = mac(c[12], q, prime[6], carry)
    c[13], carry = mac(c[13], q, prime[7], carry)
    c[14], pc = adc(c[14], pc, carry)

    c[7:16] = u256xu32(c[7:16], a[7], b)
    q = (MU * c[7]) % RADIX

    _, carry = mac(c[7], q, prime[0], 0)
    c[8], carry = mac(c[8], q, prime[1], carry)
    c[9], carry = mac(c[9], q, prime[2], carry)
    c[10], carry = mac(c[10], q, prime[3], carry)
    c[11], carry = mac(c[11], q, prime[4], carry)
    c[12], carry = mac(c[12], q, prime[5], carry)
    c[13], carry = mac(c[13], q, prime[6], carry)
    c[14], carry = mac(c[14], q, prime[7], carry)
    c[15], pc = adc(c[15], pc, carry)

    c[8] += (pc * 977)
    c[9] += pc

    return c[8:16]


def to_montgomery(a: List[int]) -> List[int]:
    '''
    Just like https://github.com/dusk-network/bls12_381/blob/ed4d87c6756c0020629edb5d8912a41e338ac85a/src/fp.rs#L251-L253;
    for better understanding read section 2.2 of https://eprint.iacr.org/2017/1057.pdf
    '''
    return montgomery_mul(a, to_radix_r(R2))


def from_montgomery(a: List[int]) -> List[int]:
    '''
    Read section 2.2 of https://eprint.iacr.org/2017/1057.pdf
    '''
    return montgomery_mul(a, to_radix_r(1))
