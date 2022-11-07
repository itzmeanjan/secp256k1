#!/usr/bin/python3

import ecdsa

# execute test cases for these many rounds
TEST_CNT: int = 1 << 5


def test_ecdsa():
    """
    Test if ECDSA keygen -> sign -> verify flow works as expected.
    """
    msg = b"this is a message !"

    skey, pkey = ecdsa.keygen()
    (r, s) = ecdsa.sign(skey, msg)
    verified = ecdsa.verify(pkey, msg, (r, s))

    assert verified, "ECDSA signature verification failed"


if __name__ == "__main__":
    print("Use `pytest` to run test cases")
