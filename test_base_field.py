#!/usr/bin/python3

from base_field import BaseField
from base_field_consts import P
from random import randint

# execute test cases for these many rounds
TEST_CNT: int = 1 << 5


def test_montgomery_repr():
    """
    Test with random secp256k1 base field elements whether convertion in between
    numeric, radix-r and Montgomery form is behaving as expected
    """
    for _ in range(TEST_CNT):
        num = randint(0, P)

        fp = BaseField.from_num(num)
        num_ = fp.to_num()

        assert num == num_, f"expeted {num}, found {num_}"


def test_base_field_multiplication():
    """
    Test if modular multiplication of two randomly generated secp256k1 base
    field elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, P)
        b = randint(0, P)
        c = (a * b) % P

        fp_a = BaseField.from_num(a)
        fp_b = BaseField.from_num(b)
        fp_c = fp_a * fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_base_field_addition():
    """
    Test if modular addition of two randomly generated secp256k1 base field
    elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, P)
        b = randint(0, P)
        c = (a + b) % P

        fp_a = BaseField.from_num(a)
        fp_b = BaseField.from_num(b)
        fp_c = fp_a + fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_base_field_subtraction():
    """
    Test if modular subtraction of two randomly generated secp256k1 base
    field elements, using Montgomery algorithm, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(0, P)
        b = randint(0, P)
        c = (a - b) % P

        fp_a = BaseField.from_num(a)
        fp_b = BaseField.from_num(b)
        fp_c = fp_a - fp_b
        c_ = fp_c.to_num()

        assert c == c_, f"expected {c}, found {c_}"


def test_base_field_inversion():
    """
    Test if modular multiplicative inversion of one randomly generated secp256k1
    base field element, in Montgomery representation, is behaving as expected
    """
    for _ in range(TEST_CNT):
        a = randint(1, P)

        fp_a = BaseField.from_num(a)
        fp_a_inv = fp_a.inv()
        fp_b = fp_a * fp_a_inv
        b = fp_b.to_num()

        assert b == 1, f"expected 1, found {b}"


if __name__ == "__main__":
    print("Use `pytest` to run test cases")
