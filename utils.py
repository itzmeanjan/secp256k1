#!/usr/bin/python3


def bit_count(num: int) -> int:
    """
    Computes bit length of an integer, same as len(bin(num)[2:])
    """
    cnt = 0
    num_ = num

    while num > 0:
        cnt += 1
        num >>= 1

    assert cnt == len(bin(num_)[2:])
    return cnt


def modulo(a: int, b: int) -> int:
    """
    Computes canonical representation of field element `a`, see https://paulmillr.com/posts/noble-secp256k1-fast-ecc/#public-keys
    """
    res = a % b
    if res >= 0:
        return res
    else:
        return b + res


def mul_inv(num: int, mod: int) -> int:
    """
    Computes multiplicative inverse of prime field element a | a * a_inv = 1,
    using extended GCD algorithm, see https://aszepieniec.github.io/stark-anatomy/basic-tools
    """
    if num == 0 or mod <= 0:
        raise Exception("can't invert multiplicative identity 0")

    old_r, r = (mod, modulo(num, mod))
    old_s, s = (1, 0)
    old_t, t = (0, 1)

    while r != 0:
        quotient = old_r // r
        old_r, r = (r, old_r - quotient * r)
        old_s, s = (s, old_s - quotient * s)
        old_t, t = (t, old_t - quotient * t)

    if old_r != 1:
        raise Exception("failed to find multiplicative inverse")

    # (old_s, old_t, old_r) | a, b, g of `ax + by = g`
    return modulo(old_s, mod)


if __name__ == "__main__":
    print("Use `utils` as library module")
