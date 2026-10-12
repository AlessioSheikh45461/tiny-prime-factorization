# tiny-prime-factorization

Factors integers up to around 10**18 into their prime factors. Uses trial division for small primes, Miller-Rabin for primality testing, and Pollard's rho (Brent variant) for splitting larger composites.

```python
from tiny_prime_factorization import factorize, is_prime

print(factorize(360))    # {2: 3, 3: 2, 5: 1}
print(is_prime(7919))   # True
print(factorize(-12))   # {2: 2, 3: 1}
```

## Why this exists

When you need exact prime factorization for integers up to about 10**18 — common in competitive programming, number-theory exercises, and debugging — you typically reach for SymPy. But SymPy is a large dependency if factorization is all you need. This library does exactly one thing with zero third-party dependencies, and the algorithms it uses (deterministic Miller-Rabin with the 12-witness set, Pollard's rho with Brent's cycle detection) are fast enough for inputs in the stated range.

The trade-off: Pollard's rho is not the fastest known factorization method for large semiprimes with two similarly-sized factors. For numbers with two 18-digit prime factors, runtime can reach a few seconds. If you need to factor arbitrary 60-digit numbers, use GMP-ECM or a quadratic sieve instead. This library targets the 99% case: factoring numbers that are mostly small primes with one or two large factors, or moderately-sized semiprimes.

## Edge cases

`factorize(0)`, `factorize(1)`, and `factorize(-1)` all return an empty dict, because none of those values have prime factors. Negative inputs are factorized by their absolute value — a sign is not a prime factor, so it does not appear in the output. The input is passed through `int()`, so numeric strings or other integer-coercible types work, but floats that would lose precision (e.g., 1e18) are your responsibility to convert correctly beforehand.
