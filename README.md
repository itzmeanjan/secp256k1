# secp256k1
Elliptic Curve Arithmetic for secp256k1 Curve 

## Prerequisites

- You must have Python3 installed & its path should be added to `$PATH`.


```bash
python3 --version
Python 3.10.8
```

- Download project dependencies from PyPI using `pip`.

```bash
python3 -m pip install --upgrade pip # may not be *always* necessary
python3 -m pip install --user -r requirements.txt
```

## Testing

For ensuring functional correctness of 

- ECDSA keygen/ sign/ verify
- Secp256k1 base & scalar field arithmetic
- Secp256k1 group arithmetic ( point addition, doubling & multiplication )

issue,

```bash
make
```

## Usage

Using ECDSA is fairly easy

- Start by importing `ecdsa` module
- Then generate a random ECDSA keypair
- Now pass message bytes along with ECDSA secret key, for signing

> **Note**

> ECDSA signature consists of two 256 -bit integers (r, s) | r, s ∈ [0, N)

> N = https://github.com/itzmeanjan/secp256k1/blob/b25f4f7/field/scalar_field_consts.py#L19-L21

- Finally use ECDSA public key & message to verify signature

```python3
>>> import ecdsa
>>> (skey, pkey) = ecdsa.keygen()
>>> msg = b'Just a message'
>>> sig = ecdsa.sign(skey, msg)
>>> verified = ecdsa.verify(pkey, msg, sig)
>>> assert verified
```
