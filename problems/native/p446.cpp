// Project Euler 446 — Retractions B, F(10^7) mod 10^9+7.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using u32 = uint32_t;
using u64 = uint64_t;
using i64 = int64_t;

static constexpr i64 MOD = 1000000007LL;

static i64 mod_pow(i64 a, i64 e, i64 mod) {
  i64 r = 1 % mod;
  a %= mod;
  while (e > 0) {
    if (e & 1) r = r * a % mod;
    a = a * a % mod;
    e >>= 1;
  }
  return r;
}

static std::vector<int> primes_upto(int limit) {
  if (limit < 2) return {};
  std::vector<char> sieve(limit / 2 + 1, 1);
  sieve[0] = 0;
  int r = (int)std::sqrt((double)limit);
  for (int x = 3; x <= r; x += 2)
    if (sieve[x / 2]) {
      int start = x * x / 2;
      for (int i = start; i < (int)sieve.size(); i += x) sieve[i] = 0;
    }
  std::vector<int> primes = {2};
  for (int i = 1; i < (int)sieve.size(); ++i)
    if (sieve[i]) primes.push_back(2 * i + 1);
  return primes;
}

static const int NONRES[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37, 41,
                             43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                             103, 107, 109, 113};

static int sqrt_minus_one(int p) {
  i64 leg = (p - 1) / 2;
  i64 quarter = (p - 1) / 4;
  for (int g : NONRES) {
    if (g >= p) break;
    if (mod_pow(g, leg, p) == p - 1) return (int)mod_pow(g, quarter, p);
  }
  for (int g = 115;; g += 2)
    if (mod_pow(g, leg, p) == p - 1) return (int)mod_pow(g, quarter, p);
}

static i64 solve(int N) {
  if (N < 1) return 0;
  int k_max = N + 1;
  auto primes = primes_upto(k_max);
  std::vector<int> p_list = {2};
  std::vector<int> r_list = {1};
  for (size_t i = 1; i < primes.size(); ++i) {
    int p = primes[i];
    if ((p & 3) == 1) {
      p_list.push_back(p);
      r_list.push_back(sqrt_minus_one(p));
    }
  }
  i64 correction_even = 5 * mod_pow(9, MOD - 2, MOD) % MOD;
  const int block_size = 1000000;
  i64 total = 0;
  i64 prev2_P = -1, prev1_P = -1, prev2_C = -1, prev1_C = -1;

  for (int L = 0; L <= k_max; L += block_size) {
    int R = std::min(k_max + 1, L + block_size);
    int size = R - L;
    std::vector<u64> rem(size);
    std::vector<u32> prod(size, 1), cmod(size);
    i64 k = L;
    i64 v = k * k + 1;
    for (int i = 0; i < size; ++i) {
      rem[i] = (u64)v;
      cmod[i] = (u32)(v % MOD);
      v += 2 * k + 1;
      ++k;
    }
    for (size_t j = 0; j < p_list.size(); ++j) {
      int p = p_list[j];
      int r = r_list[j];
      if (p == 2) {
        int start = (1 - L) & 1;
        for (int idx = start; idx < size; idx += 2) {
          rem[idx] /= 2;
          prod[idx] = (u32)((i64)prod[idx] * 3 % MOD);
        }
        continue;
      }
      auto process_root = [&](int root) {
        int start = ((root - L) % p + p) % p;
        for (int idx = start; idx < size; idx += p) {
          u64 x = rem[idx] / p;
          i64 pe = p;
          while (x % (u64)p == 0) {
            x /= p;
            pe *= p;
          }
          rem[idx] = x;
          i64 t = pe + 1;
          if (t >= MOD) t %= MOD;
          prod[idx] = (u32)((i64)prod[idx] * t % MOD);
        }
      };
      process_root(r);
      process_root(p - r);
    }
    for (int i = 0; i < size; ++i) {
      if (rem[i] > 1) {
        i64 t = (i64)rem[i] + 1;
        if (t >= MOD) t %= MOD;
        prod[i] = (u32)((i64)prod[i] * t % MOD);
      }
    }
    for (int i = 0; i < size; ++i) {
      int kk = L + i;
      i64 Pk = prod[i];
      i64 Ck = cmod[i];
      if (prev2_P >= 0) {
        int n = kk - 1;
        if (n <= N) {
          i64 pm = prev2_P * Pk % MOD;
          if ((n & 1) == 0) pm = pm * correction_even % MOD;
          i64 mm = prev2_C * Ck % MOD;
          total = (total + pm - mm) % MOD;
          if (total < 0) total += MOD;
        }
      }
      prev2_P = prev1_P;
      prev1_P = Pk;
      prev2_C = prev1_C;
      prev1_C = Ck;
    }
  }
  return total % MOD;
}

extern "C" long long pe446_answer(void) { return solve(10000000); }
