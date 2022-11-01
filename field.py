#!/usr/bin/python3

from functools import reduce
from typing_extensions import Self
from utils import *


class Field:
    """
    A field element of secp256k1 prime field, kept in Montgomery form
    """

    def __init__(self, limbs: List[int]):
        self._limbs = limbs

    @classmethod
    def from_num(cls, num: int) -> Self:
        """
        Given an element of secp256k1 prime field as integer, this routine returns
        it as field element in Montgomery form
        """
        return cls.from_radix_r(to_radix_r(num))

    def to_num(self) -> int:
        """
        Given field element in Montgomery form, this routine computes field element
        as integer
        """
        return from_radix_r(self.to_radix_r())

    @classmethod
    def from_radix_r(cls, limbs: List[int]) -> Self:
        """
        Given an element of secp256k1 prime field in radix-r form, this routine returns
        it as field element in Montgomery form | r = 2^32
        """
        return cls(to_montgomery(limbs))

    def to_radix_r(self) -> List[int]:
        """
        Given field element in Montgomery form, this routine computes field element
        as radix-r form | r = 2^32
        """
        return from_montgomery(self._limbs)

    def __eq__(self, rhs: Self) -> bool:
        tmp = [bool(self._limbs[i] ^ rhs._limbs[i]) for i in range(LIMB_COUNT)]
        return not reduce(lambda acc, cur: acc | cur, tmp, False)

    def __mul__(self, rhs: Self) -> Self:
        """
        Modular multiplication of two secp256k1 field elements, input/ output in Montgomery form
        """
        return Field(montgomery_mul(self._limbs, rhs._limbs))

    def __add__(self, rhs: Self) -> Self:
        """
        Modular addition of two secp256k1 field elements, input/ output in Montgomery form
        """
        c = [0] * LIMB_COUNT
        carry = 0

        c[0], carry = adc(self._limbs[0], rhs._limbs[0], carry)
        c[1], carry = adc(self._limbs[1], rhs._limbs[1], carry)
        c[2], carry = adc(self._limbs[2], rhs._limbs[2], carry)
        c[3], carry = adc(self._limbs[3], rhs._limbs[3], carry)
        c[4], carry = adc(self._limbs[4], rhs._limbs[4], carry)
        c[5], carry = adc(self._limbs[5], rhs._limbs[5], carry)
        c[6], carry = adc(self._limbs[6], rhs._limbs[6], carry)
        c[7], carry = adc(self._limbs[7], rhs._limbs[7], carry)

        c[0] += carry * 977
        c[1] += carry

        return Field(c)

    def __neg__(self) -> Self:
        """
        Negates a secp256k1 field element such that a + b = 0, if b = -a
        """
        P = to_radix_r(PRIME)

        c = [0] * LIMB_COUNT
        borrow = 0

        c[0], borrow = sbb(P[0], self._limbs[0], borrow)
        c[1], borrow = sbb(P[1], self._limbs[1], borrow)
        c[2], borrow = sbb(P[2], self._limbs[2], borrow)
        c[3], borrow = sbb(P[3], self._limbs[3], borrow)
        c[4], borrow = sbb(P[4], self._limbs[4], borrow)
        c[5], borrow = sbb(P[5], self._limbs[5], borrow)
        c[6], borrow = sbb(P[6], self._limbs[6], borrow)
        c[7], _ = sbb(P[7], self._limbs[7], borrow)

        return Field(c)

    def __sub__(self, rhs: Self) -> Self:
        """
        Modular subtraction of two secp256k1 elements, input/ output in Montgomery form
        """
        return self + (-rhs)

    def inv(self) -> Self:
        """
        Computes multiplicative inverse of secp256k1 field element
        """

        def pow(a: List[int], b: List[int]) -> List[int]:
            res = to_radix_r(R)

            for i in reversed(b):
                for j in reversed(range(RADIX_BIT_LEN)):
                    res = montgomery_mul(res, res)

                    if (i >> j) & 1:
                        res = montgomery_mul(res, a)

            return res

        return Field(pow(self._limbs, to_radix_r(PRIME - 2)))

    def __repr__(self) -> str:
        """
        Pretty print on console
        """
        return f"Fp({self.to_num()}, {PRIME})"

    def __str__(self) -> str:
        """
        Display when printed to stdout/ file
        """
        return f"{self.to_num()}"


if __name__ == "__main__":
    print("Use `field` as library module")
