#!/usr/bin/python3

from .scalar_field_consts import *
from typing import List, Tuple


def to_radix_r(num: int) -> List[int]:
    """
    Converts secp256k1 scalar field element ( represented as integer ) to radix-r
    interleaved representation | r = 2^32
    """
    limbs = [0] * LIMB_COUNT
    idx = 0
    while num > 0:
        limbs[idx] = num % RADIX
        num //= RADIX
        idx += 1
    return limbs


def from_radix_r(limbs: List[int]) -> int:
    """
    Converts radix-r interleaved representation of a secp256k1 scalar field element
    to integer | r = 2^32
    """
    cnt = len(limbs)
    num = 0
    for idx in range(cnt - 1, -1, -1):
        num = num * RADIX + limbs[idx]
    return num


def adc(a: int, b: int, carry: int) -> Tuple[int, int]:
    """
    See https://github.com/dusk-network/bls12_381/blob/ed4d87c/src/util.rs#L1-L6
    """
    tmp = a + b + carry
    return tmp & (RADIX - 1), tmp >> RADIX_BIT_LEN


def mac(a: int, b: int, c: int, carry: int) -> Tuple[int, int]:
    """
    See https://github.com/dusk-network/bls12_381/blob/ed4d87c/src/util.rs#L15-L20
    """
    tmp = a + (b * c) + carry
    return tmp & (RADIX - 1), tmp >> RADIX_BIT_LEN


def sbb(a: int, b: int, borrow: int) -> Tuple[int, int]:
    """
    See https://github.com/dusk-network/bls12_381/blob/2c679a2/src/util.rs#L8-L13
    """
    wrap_at = RADIX << RADIX_BIT_LEN
    tmp = (a - (b + (borrow >> (RADIX_BIT_LEN - 1)))) % wrap_at
    return tmp & (RADIX - 1), tmp >> RADIX_BIT_LEN


def bitwise_not(a: int) -> int:
    """
    Same as `!a` in C
    """
    return RADIX - 1 - a


def u256xu32(a: List[int], b: int, c: List[int]) -> List[int]:
    """
    Multiplies a 256 -bit number by another 32 -bit number, resulting into a 288 -bit one.

    Collects inspiration from https://github.com/dusk-network/bls12_381/blob/ed4d87c/src/fp.rs#L517-L522
    """
    assert len(a) == (LIMB_COUNT + 1) and len(c) == LIMB_COUNT

    a_ = [0] * len(c)
    carry = 0

    a_[0], carry = mac(0, c[0], b, carry)
    a_[1], carry = mac(0, c[1], b, carry)
    a_[2], carry = mac(0, c[2], b, carry)
    a_[3], carry = mac(0, c[3], b, carry)
    a_[4], carry = mac(0, c[4], b, carry)
    a_[5], carry = mac(0, c[5], b, carry)
    a_[6], carry = mac(0, c[6], b, carry)
    a_[7], a[8] = mac(0, c[7], b, carry)

    carry = 0

    a[0], carry = adc(a[0], a_[0], carry)
    a[1], carry = adc(a[1], a_[1], carry)
    a[2], carry = adc(a[2], a_[2], carry)
    a[3], carry = adc(a[3], a_[3], carry)
    a[4], carry = adc(a[4], a_[4], carry)
    a[5], carry = adc(a[5], a_[5], carry)
    a[6], carry = adc(a[6], a_[6], carry)
    a[7], carry = adc(a[7], a_[7], carry)

    a[8] += carry

    return a


def montgomery_mul(a: List[int], b: List[int]) -> List[int]:
    """
    Multiplies two 256 -bit numbers ( in Montgomery form ), resulting into a 256 -bit
    number ( in Montgomery form ), reduced by secp256k1 scalar field prime.

    Inspired by https://github.com/dusk-network/bls12_381/blob/ed4d87c/src/fp.rs#L437-L560
    and algorithm 2 of https://eprint.iacr.org/2017/1057.pdf
    """
    assert len(a) == len(b)

    prime = to_radix_r(N)
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

    c[8] += pc * 801750719
    c[9] += pc * 1076732275
    c[10] += pc * 1354194884
    c[11] += pc * 1162945305
    c[12] += pc

    return c[8:16]


def to_montgomery(a: List[int]) -> List[int]:
    """
    Converts a radix-r form secp256k1 scalar field element to Montgomery form.

    Just like https://github.com/dusk-network/bls12_381/blob/ed4d87c/src/fp.rs#L251-L253;
    for better understanding read section 2.2 of https://eprint.iacr.org/2017/1057.pdf
    """
    return montgomery_mul(a, to_radix_r(R2))


def from_montgomery(a: List[int]) -> List[int]:
    """
    Converts a Montgomery form secp256k1 scalar field element to radix-r form.

    Read section 2.2 of https://eprint.iacr.org/2017/1057.pdf
    """
    return montgomery_mul(a, to_radix_r(1))
