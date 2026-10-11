"""Integer factorization for numbers up to roughly 10**18.

The pipeline is:
    1. Strip powers of 2.
    2. Trial-divide by a handful of tiny primes (cost is trivial and it
       shrinks the input before the heavier machinery runs).
    3. Miller-Rabin for remaining composites.
    4. Pollard's rho with Brent's cycle-finding for the actual splitting.

All arithmetic uses Python's big integers, so nothing here can overflow;
the 10**18 ceiling is about keeping Pollard's rho runtime sane on a single
thread rather than a mathematical hard limit.
"""

from __future__ import annotations

import math
import random

# Deterministic Miller-Rabin witnesses for n < 3,317,044,064,679,887,385,961,981
# (Jim Sinclair / the well-known small set). This covers our entire stated
# range and makes is_prime a true primality test, not a probabilistic one.
_MR_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

# Small primes for the trial-division pre-pass. Tuned so that we dispose of
# the vast majority of "easy" factors cheaply before invoking rho.
_SMALL_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                 61, 67, 71, 73, 79, 83, 89, 97)


def is_prime(n: int) -> bool:
    """Return True if *n* is prime, False otherwise.

    Uses the deterministic 12-witness Miller-Rabin set, which is exact for
    every integer below ~3.3 * 10**24 — comfortably beyond our 10**18 target.
    Input is normalized so that callers may pass any integer without
    surprises.
    """
    n = int(n)
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1

    for a in _MR_WITNESSES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n: int) -> int:
    """Return a non-trivial factor of composite *n* via Brent's rho.

    *n* must be composite and greater than 1 (the caller guarantees this).
    The function always returns a proper divisor in the range (1, n); it
    never returns 1 or *n* itself.
    """
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3

    while True:
        # Fresh random parameters each outer iteration until we find a factor.
        # Using math.gcd instead of tracking remainders keeps the inner loop
        # free of conditional branches that would otherwise dominate.
        y = random.randrange(2, n - 1)
        c = random.randrange(1, n - 1)
        m = random.randrange(1, n - 1)
        g = 1
        r = 1
        q = 1
        x = ys = 0

        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                bound = min(m, r - k)
                for _ in range(bound):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2

        if g == n:
            # Backtrack to the last saved state and step linearly. This is
            # the textbook recovery path when the batch GCD collapsed to n.
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)

        if 1 < g < n:
            return g
        # else: g == n again — extremely rare; retry with new parameters.


def factorize(n: int) -> dict[int, int]:
    """Return the prime factorization of *n* as a mapping {prime: exponent}.

    Returns an empty dict for inputs in {-1, 0, 1}. Negative inputs are
    factorized by their absolute value, because a sign is not a prime factor.
    """
    n = int(n)
    if n < 0:
        n = -n
    if n <= 1:
        return {}

    factors: dict[int, int] = {}

    # 1) Strip the factor of 2 explicitly so the rho loop can assume odd n.
    if n % 2 == 0:
        e = 0
        while n % 2 == 0:
            n //= 2
            e += 1
        factors[2] = e

    # 2) Cheap trial division for small odd primes.
    for p in _SMALL_PRIMES:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            factors[p] = e

    # 3) Anything left that is prime needs no further work.
    if n == 1:
        return factors
    if is_prime(n):
        factors[n] = factors.get(n, 0) + 1
        return factors

    # 4) Otherwise split with rho until every piece is prime.
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            factors[m] = factors.get(m, 0) + 1
            continue
        d = _pollard_rho(m)
        stack.append(d)
        stack.append(m // d)

    return factors
