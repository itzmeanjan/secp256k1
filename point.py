#!/usr/bin/python3

from typing_extensions import Self
from field import BaseField
from typing import Tuple


class Point:
    """
    A secp256k1 elliptic curve point, kept in projective coordinate system
    """

    def __init__(self, x: BaseField, y: BaseField, z: BaseField):
        self._x = x
        self._y = y
        self._z = z

    def __str__(self) -> str:
        """
        Display when printed to stdout/ file
        """
        return f"{self._x}, {self._y}, {self._z}"

    def __repr__(self) -> str:
        """
        Pretty print on console
        """
        return f"Point({self._x}, {self._y}, {self._z})"

    def __eq__(self, rhs: Self) -> bool:
        """
        First converts both of elliptic curve points to affine coordinate system & then
        checks for equality
        """
        x1, y1 = self.toAffine()
        x2, y2 = rhs.toAffine()

        return (x1 == x2) & (y1 == y2)

    @classmethod
    def zero(cls) -> Self:
        """
        Identity element of group, see https://github.com/dusk-network/bls12_381/blob/2c679a2/src/g1.rs#L587-L593
        """
        return cls(BaseField.from_num(0), BaseField.from_num(1), BaseField.from_num(0))

    @classmethod
    def fromAffine(cls, x: BaseField, y: BaseField) -> Self:
        """
        Given affine coordinate of secp256k1 elliptic curve point, this routine
        returns equivalent point in projective coordinate system
        """
        return Point(x, y, BaseField.from_num(1))

    def toAffine(self) -> Tuple[BaseField, BaseField]:
        """
        Given projective coordinate of secp256k1 elliptic curve point, this routine
        computes equivalent point in affine coordinate system
        """
        inv_z = self._z.inv()

        x = self._x * inv_z
        y = self._y * inv_z

        return x, y

    def __add__(self, rhs: Self) -> Self:
        """
        Adds two elliptic curve points in projective coordinate system, using exception-free addition
        formula provided in algorithm 7 of https://eprint.iacr.org/2015/1060.pdf
        """
        x1, y1, z1 = self._x, self._y, self._z
        x2, y2, z2 = rhs._x, rhs._y, rhs._z

        b = 7
        b3 = BaseField.from_num(3 * b)

        t0 = x1 * x2
        t1 = y1 * y2
        t2 = z1 * z2

        t3 = x1 + y1
        t4 = x2 + y2
        t3 = t3 * t4

        t4 = t0 + t1
        t3 = t3 - t4
        t4 = y1 + z1

        x3 = y2 + z2
        t4 = t4 * x3
        x3 = t1 + t2

        t4 = t4 - x3
        x3 = x1 + z1
        y3 = x2 + z2

        x3 = x3 * y3
        y3 = t0 + t2
        y3 = x3 - y3

        x3 = t0 + t0
        t0 = x3 + t0
        t2 = b3 * t2

        z3 = t1 + t2
        t1 = t1 - t2
        y3 = b3 * y3

        x3 = t4 * y3
        t2 = t3 * t1
        x3 = t2 - x3

        y3 = y3 * t0
        t1 = t1 * z3
        y3 = t1 + y3

        t0 = t0 * t3
        z3 = z3 * t4
        z3 = z3 + t0

        return Point(x3, y3, z3)

    def __neg__(self) -> Self:
        """
        Negates elliptic curve point in projective coordinate system by changing sign of Y -coordinate
        """
        return Point(self._x, -self._y, self._z)

    def __sub__(self, rhs: Self) -> Self:
        """
        Subtracts two elliptic curve points in projective coordinate system, using exception-free addition
        formula provided in algorithm 7 of https://eprint.iacr.org/2015/1060.pdf, while first negating right hand
        side operand
        """
        return self + (-rhs)

    def double(self) -> Self:
        """
        Doubles elliptic curve point `p` in projective coordinate system, using exception-free doubling
        formula provided in algorithm 9 of https://eprint.iacr.org/2015/1060.pdf | return value = p + p
        """
        x, y, z = self._x, self._y, self._z

        b = 7
        b3 = BaseField.from_num(3 * b)

        t0 = y * y
        z3 = t0 + t0
        z3 = z3 + z3

        z3 = z3 + z3
        t1 = y * z
        t2 = z * z

        t2 = b3 * t2
        x3 = t2 * z3
        y3 = t0 + t2

        z3 = t1 * z3
        t1 = t2 + t2
        t2 = t1 + t2

        t0 = t0 - t2
        y3 = t0 * y3
        y3 = x3 + y3

        t1 = x * y
        x3 = t0 * t1
        x3 = x3 + x3

        return Point(x3, y3, z3)

    def mulScalar(self, scalar: int) -> Self:
        res = Point.zero()
        tmp = self

        idx = 0
        while idx < 256:
            if (scalar >> idx) & 1:
                res += tmp

            tmp = tmp.double()
            idx += 1

        return res


if __name__ == "__main__":
    print("Use `point` as library module")
