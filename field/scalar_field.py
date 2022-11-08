#!/usr/bin/python3

from functools import reduce
from typing_extensions import Self
from .scalar_field_utils import *


class ScalarField:
    """
    A secp256k1 scalar field element, kept in Montgomery form
    """

    def __init__(self, limbs: List[int]):
        self._limbs = limbs

    @classmethod
    def from_num(cls, num: int) -> Self:
        """
        Given an element of secp256k1 scalar field as integer, this routine converts
        it to Montgomery form
        """
        return cls.from_radix_r(to_radix_r(num))

    def to_num(self) -> int:
        """
        Given secp256k1 scalar field element in Montgomery form, this routine computes
        it as an integer
        """
        return from_radix_r(self.to_radix_r())

    @classmethod
    def from_radix_r(cls, limbs: List[int]) -> Self:
        """
        Given an element of secp256k1 scalar field in radix-r form, this routine returns
        it in Montgomery form | r = 2^32
        """
        return cls(to_montgomery(limbs))

    def to_radix_r(self) -> List[int]:
        """
        Given a secp256k1 scalar field element in Montgomery form, this routine computes
        it in radix-r form | r = 2^32
        """
        return from_montgomery(self._limbs)

    def __eq__(self, rhs: Self) -> bool:
        """
        Checks equality of two elements of secp256k1 scalar field, when they're
        kept in their Montgomery form
        """
        tmp = [bool(self._limbs[i] ^ rhs._limbs[i]) for i in range(LIMB_COUNT)]
        return not reduce(lambda acc, cur: acc | cur, tmp, False)

    def __mul__(self, rhs: Self) -> Self:
        """
        Modular multiplication of two secp256k1 scalar field elements, input/ output
        expected to be in Montgomery form
        """
        return ScalarField(montgomery_mul(self._limbs, rhs._limbs))

    def __add__(self, rhs: Self) -> Self:
        """
        Modular addition of two secp256k1 scalar field elements, input/ output in Montgomery form
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

        one = [801750719, 1076732275, 1354194884, 1162945305, 1, 0, 0, 0]
        one = [i * carry for i in one]

        carry = 0
        c[0], carry = adc(c[0], one[0], carry)
        c[1], carry = adc(c[1], one[1], carry)
        c[2], carry = adc(c[2], one[2], carry)
        c[3], carry = adc(c[3], one[3], carry)
        c[4], carry = adc(c[4], one[4], carry)
        c[5], carry = adc(c[5], one[5], carry)
        c[6], carry = adc(c[6], one[6], carry)
        c[7], _ = adc(c[7], one[7], carry)

        return ScalarField(c)

    def __neg__(self) -> Self:
        """
        Negates a secp256k1 scalar element such that a + b = 0, if b = -a
        """
        P_ = to_radix_r(N)

        c = [0] * LIMB_COUNT
        borrow = 0

        c[0], borrow = sbb(P_[0], self._limbs[0], borrow)
        c[1], borrow = sbb(P_[1], self._limbs[1], borrow)
        c[2], borrow = sbb(P_[2], self._limbs[2], borrow)
        c[3], borrow = sbb(P_[3], self._limbs[3], borrow)
        c[4], borrow = sbb(P_[4], self._limbs[4], borrow)
        c[5], borrow = sbb(P_[5], self._limbs[5], borrow)
        c[6], borrow = sbb(P_[6], self._limbs[6], borrow)
        c[7], _ = sbb(P_[7], self._limbs[7], borrow)

        return ScalarField(c)

    def __sub__(self, rhs: Self) -> Self:
        """
        Modular subtraction of two secp256k1 scalar field elements, input/ output in Montgomery form
        """
        return self + (-rhs)

    def inv(self) -> Self:
        """
        Computes multiplicative inverse of a secp256k1 scalar field element. If operand is 0,
        returns 0, because it's not possible to compute multiplicative inverse of zero element.
        """

        def pow(a: List[int], b: List[int]) -> List[int]:
            res = to_radix_r(R)

            for i in reversed(b):
                for j in reversed(range(RADIX_BIT_LEN)):
                    res = montgomery_mul(res, res)

                    if (i >> j) & 1:
                        res = montgomery_mul(res, a)

            return res

        return ScalarField(pow(self._limbs, to_radix_r(N - 2)))

    def __repr__(self) -> str:
        """
        Pretty print on console
        """
        return f"Fp({self.to_num()}, {N})"

    def __str__(self) -> str:
        """
        Display when printed to stdout/ file
        """
        return str(self.to_num())
