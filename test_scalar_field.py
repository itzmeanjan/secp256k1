#!/usr/bin/python3

from scalar_field import ScalarField
from scalar_field_consts import N
from random import randint

# execute test cases for these many rounds
TEST_CNT: int = 1 << 10


def test_montgomery_repr():
    """
    Test with random secp256k1 scalar field elements whether convertion in between
    numeric, radix-r and Montgomery form is behaving as expected
    """
    for _ in range(TEST_CNT):
        num = randint(0, N)

        fp = ScalarField.from_num(num)
        num_ = fp.to_num()

        assert num == num_, f"expeted {num}, found {num_}"


def test_scalar_field_multiplication():
    """
    Test if modular multiplication of two randomly generated secp256k1 scalar
    field elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, N)
        b = randint(0, N)
        c = (a * b) % N

        fp_a = ScalarField.from_num(a)
        fp_b = ScalarField.from_num(b)
        fp_c = fp_a * fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_scalar_field_addition():
    """
    Test if modular addition of two randomly generated secp256k1 scalar field
    elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, N)
        b = randint(0, N)
        c = (a + b) % N

        fp_a = ScalarField.from_num(a)
        fp_b = ScalarField.from_num(b)
        fp_c = fp_a + fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_scalar_field_subtraction():
    """
    Test if modular subtraction of two randomly generated secp256k1 scalar
    field elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, N)
        b = randint(0, N)
        c = (a - b) % N

        fp_a = ScalarField.from_num(a)
        fp_b = ScalarField.from_num(b)
        fp_c = fp_a - fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_scalar_field_inversion():
    """
    Test if modular multiplicative inversion of one randomly generated secp256k1
    scalar field element, in Montgomery representation, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(1, N)

        fp_a = ScalarField.from_num(a)
        fp_a_inv = fp_a.inv()
        fp_b = fp_a * fp_a_inv
        b = fp_b.to_num()

        assert b == 1, f"expected 1, found {b}"


if __name__ == "__main__":
    print("Use `pytest` to run test cases")
