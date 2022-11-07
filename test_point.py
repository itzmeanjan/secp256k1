#!/usr/bin/python3

from point import Point, BaseField
from scalar_field_consts import Gx, Gy
from random import randint

# execute test cases for these many rounds
TEST_CNT: int = 1 << 5


def random_point() -> Point:
    """
    Routine for generating random point on secp256k1 elliptic curve, while starting
    with curve generator points
    """
    gen = Point.fromAffine(BaseField.from_num(Gx), BaseField.from_num(Gy))

    itr = randint(0, 1 << 5)
    for _ in range(itr):
        gen = gen.double()

    return gen


def test_point_addition():
    """
    Test if elliptic curve point addition of two randomly generated secp256k1 points, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = random_point()
        b = random_point()

        c = a + b
        d = c - a

        assert b == d, f"expected {b}, found {d}"


def test_point_doubling():
    """
    Test if elliptic curve point doubling of randomly generated secp256k1 point, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = random_point()

        b = a + a
        c = a.double()

        assert b == c, f"expected {b}, found {c}"


if __name__ == "__main__":
    print("Use `pytest` to run test cases")
