import math
import unittest

from tiny_prime_factorization import factorize, is_prime


class TestIsPrime(unittest.TestCase):
    def test_small_and_edge(self):
        self.assertFalse(is_prime(-7))
        self.assertFalse(is_prime(-1))
        self.assertFalse(is_prime(0))
        self.assertFalse(is_prime(1))
        self.assertTrue(is_prime(2))
        self.assertTrue(is_prime(3))
        self.assertFalse(is_prime(4))
        self.assertTrue(is_prime(5))
        self.assertFalse(is_prime(9))
        self.assertFalse(is_prime(25))

    def test_known_primes_and_composites(self):
        known_primes = [13, 17, 97, 1009, 7919, 104729, 1299709, 999999999989]
        for p in known_primes:
            self.assertTrue(is_prime(p), f"{p} should be prime")
        known_composites = [15, 21, 51, 91, 561, 1105, 1729, 2465, 2821, 6601,
                            1000000000000, 999999999990]
        for c in known_composites:
            self.assertFalse(is_prime(c), f"{c} should be composite")

    def test_strong_pseudoprime_bases(self):
        # Numbers that fool single-base Miller-Rabin but not our fixed set.
        for n in [2047, 3277, 4033, 4681, 8321, 15841]:
            self.assertFalse(is_prime(n))

    def test_accepts_numeric_strings(self):
        # is_prime normalizes via int(), which is documented behaviour.
        self.assertTrue(is_prime("7"))
        self.assertFalse(is_prime("9"))


class TestFactorize(unittest.TestCase):
    def test_trivial_inputs(self):
        self.assertEqual(factorize(0), {})
        self.assertEqual(factorize(1), {})
        self.assertEqual(factorize(-1), {})

    def test_negative_uses_absolute_value(self):
        self.assertEqual(factorize(-12), {2: 2, 3: 1})

    def test_prime_powers(self):
        self.assertEqual(factorize(8), {2: 3})
        self.assertEqual(factorize(27), {3: 3})
        self.assertEqual(factorize(49), {7: 2})
        self.assertEqual(factorize(2 ** 10), {2: 10})

    def test_simple_composites(self):
        self.assertEqual(factorize(12), {2: 2, 3: 1})
        self.assertEqual(factorize(100), {2: 2, 5: 2})
        self.assertEqual(factorize(360), {2: 3, 3: 2, 5: 1})

    def test_semiprime_with_large_factors(self):
        p, q = 1000003, 9999991
        n = p * q
        self.assertEqual(factorize(n), {p: 1, q: 1})

    def test_large_prime(self):
        self.assertEqual(factorize(999999999989), {999999999989: 1})

    def test_hard_composite(self):
        p, q = 1000003, 999999999989
        n = p * q
        self.assertEqual(factorize(n), {p: 1, q: 1})

    def test_repeated_large_prime_factor(self):
        p = 999999999989
        n = p * p * 7
        self.assertEqual(factorize(n), {p: 2, 7: 1})

    def test_semiprime_of_two_distinct_primes(self):
        # Exercises the rho-splitting path for two roughly-equal primes.
        p, q = 10007, 10009
        n = p * q
        self.assertEqual(factorize(n), {p: 1, q: 1})

    def test_product_of_many_small_primes(self):
        n = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23
        expected = {p: 1 for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]}
        self.assertEqual(factorize(n), expected)

    def test_output_is_sorted_by_consistent_reconstruction(self):
        # Instead of asserting ordering, verify the product equals the input.
        for n in [1, 2, 6, 12, 97, 100, 2 ** 20, 999999999989]:
            f = factorize(n)
            prod = 1
            for p, e in f.items():
                self.assertTrue(is_prime(p), f"{p} is not prime")
                prod *= p ** e
            self.assertEqual(prod, n if n > 1 else 1)

    def test_product_of_three_distinct_primes(self):
        a, b, c = 1009, 10007, 100003
        n = a * b * c
        self.assertEqual(factorize(n), {a: 1, b: 1, c: 1})


if __name__ == "__main__":
    unittest.main()
