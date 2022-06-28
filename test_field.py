#!/usr/bin/python3

from field import Field
from consts import TEST_CNT, PRIME
from random import randint


def test_to_and_from_montgomery_repr():
    '''
    Test with random secp256k1 field elements whether convertion in between numeric, radix-r and Montgomery forms
    is behaving as expected
    '''
    for _ in range(TEST_CNT):
        num = randint(0, PRIME)

        fp = Field.from_num(num)
        num_ = fp.to_num()

        assert num == num_, f'expeted {num}, found {num_}'


def test_field_multiplication():
    '''
    Test if modular multiplication of two randomly generated secp256k1 prime field elements, using Montgomery algorithm,
    is behaving as expected
    '''
    for _ in range(TEST_CNT):
        a = randint(0, PRIME)
        b = randint(0, PRIME)
        c = (a * b) % PRIME

        fp_a = Field.from_num(a)
        fp_b = Field.from_num(b)
        fp_c = fp_a * fp_b
        c_ = fp_c.to_num()

        assert c == c_, f'expected {c}, found {c_}'


def test_field_addition():
    '''
    Test if modular addition of two randomly generated secp256k1 prime field elements, using Montgomery algorithm,
    is behaving as expected
    '''
    for _ in range(TEST_CNT):
        a = randint(0, PRIME)
        b = randint(0, PRIME)
        c = (a + b) % PRIME

        fp_a = Field.from_num(a)
        fp_b = Field.from_num(b)
        fp_c = fp_a + fp_b
        c_ = fp_c.to_num()

        assert c == c_, f'expected {c}, found {c_}'


def test_field_subtraction():
    '''
    Test if modular subtraction of two randomly generated secp256k1 prime field elements, using Montgomery algorithm,
    is behaving as expected
    '''
    for _ in range(TEST_CNT):
        a = randint(0, PRIME)
        b = randint(0, PRIME)
        c = (a - b) % PRIME

        fp_a = Field.from_num(a)
        fp_b = Field.from_num(b)
        fp_c = fp_a - fp_b
        c_ = fp_c.to_num()

        assert c == c_, f'expected {c}, found {c_}'


def test_field_inversion():
    '''
    Test if modular multiplicative inversion of one randomly generated secp256k1 prime field element, in Montgomery representation,
    is behaving as expected
    '''
    for _ in range(TEST_CNT):
        a = randint(0, PRIME)

        fp_a = Field.from_num(a)
        fp_a_inv = fp_a.inv()
        fp_b = fp_a * fp_a_inv
        b = fp_b.to_num()

        assert b == 1, f'expected 1, found {b}'
