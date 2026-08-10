// Project Euler 415 — Titanic Sets, N=10^11 mod 10^8.
#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

using i64 = int64_t;

static constexpr i64 MOD = 100000000LL;
static constexpr i64 DEFAULT_N = 100000000000LL;
static constexpr int PRECOMPUTE = 5000000;

static i64 norm(i64 x) {
  x %= MOD;
  if (x < 0) x += MOD;
  return x;
}

static i64 mod_pow(i64 a, __int128 e) {
  i64 r = 1 % MOD;
  a %= MOD;
  while (e > 0) {
    if (e & 1) r = (i64)((__int128)r * a % MOD);
    a = (i64)((__int128)a * a % MOD);
    e >>= 1;
  }
  return r;
}

// Exact n(n+1)/2 mod MOD
static i64 s1(i64 n) {
  if (n <= 0) return 0;
  i64 a = n, b = n + 1;
  if (a % 2 == 0) a /= 2;
  else b /= 2;
  return (a % MOD) * (b % MOD) % MOD;
}

// Exact n(n+1)(2n+1)/6 mod MOD
static i64 s2(i64 n) {
  if (n <= 0) return 0;
  i64 a = n, b = n + 1, c = 2 * n + 1;
  if (a % 2 == 0) a /= 2;
  else b /= 2;
  if (a % 3 == 0) a /= 3;
  else if (b % 3 == 0) b /= 3;
  else c /= 3;
  return (a % MOD) * (b % MOD) % MOD * (c % MOD) % MOD;
}

static i64 s3(i64 n) {
  if (n <= 0) return 0;
  i64 t = s1(n);
  return t * t % MOD;
}

static i64 range_s1(i64 lo, i64 hi) { return norm(s1(hi) - s1(lo - 1)); }
static i64 range_s2(i64 lo, i64 hi) { return norm(s2(hi) - s2(lo - 1)); }
static i64 range_s3(i64 lo, i64 hi) { return norm(s3(hi) - s3(lo - 1)); }

static i64 pref_k_pow2(i64 n, i64 pow2_next) {
  if (n < 0) return 0;
  return norm(((n - 1) % MOD) * pow2_next + 2);
}
static i64 pref_k2_pow2(i64 n, i64 pow2_next) {
  if (n < 0) return 0;
  i64 k = n % MOD;
  return norm(((k * k % MOD - 2 * k % MOD + 3) % MOD) * pow2_next - 6);
}

struct TotientSums {
  int limit;
  std::vector<i64> pref0, pref1, pref2;
  std::unordered_map<i64, i64> cache0, cache1, cache2;

  explicit TotientSums(i64 max_n) {
    limit = (int)std::min<i64>(max_n, PRECOMPUTE);
    if (limit < 1) limit = 1;
    std::vector<i64> phi(limit + 1);
    std::vector<int> primes;
    std::vector<char> composite(limit + 1, 0);
    phi[1] = 1;
    for (int x = 2; x <= limit; ++x) {
      if (!composite[x]) {
        primes.push_back(x);
        phi[x] = x - 1;
      }
      i64 phix = phi[x];
      for (int p : primes) {
        i64 y = (i64)x * p;
        if (y > limit) break;
        composite[(int)y] = 1;
        if (x % p == 0) {
          phi[(int)y] = phix * p;
          break;
        }
        phi[(int)y] = phix * (p - 1);
      }
    }
    pref0.assign(limit + 1, 0);
    pref1.assign(limit + 1, 0);
    pref2.assign(limit + 1, 0);
    for (int x = 1; x <= limit; ++x) {
      i64 ph = phi[x] % MOD;
      i64 xm = x % MOD;
      pref0[x] = (pref0[x - 1] + ph) % MOD;
      pref1[x] = (pref1[x - 1] + xm * ph) % MOD;
      pref2[x] = (pref2[x - 1] + xm * xm % MOD * ph) % MOD;
    }
  }

  i64 phi_sum(i64 n) {
    if (n <= limit) return pref0[(int)n];
    auto it = cache0.find(n);
    if (it != cache0.end()) return it->second;
    i64 total = s1(n);
    for (i64 lo = 2; lo <= n;) {
      i64 q = n / lo;
      i64 hi = n / q;
      total = norm(total - ((hi - lo + 1) % MOD) * phi_sum(q) % MOD);
      lo = hi + 1;
    }
    return cache0[n] = total;
  }
  i64 i_phi(i64 n) {
    if (n <= limit) return pref1[(int)n];
    auto it = cache1.find(n);
    if (it != cache1.end()) return it->second;
    i64 total = s2(n);
    for (i64 lo = 2; lo <= n;) {
      i64 q = n / lo;
      i64 hi = n / q;
      total = norm(total - range_s1(lo, hi) * i_phi(q) % MOD);
      lo = hi + 1;
    }
    return cache1[n] = total;
  }
  i64 i2_phi(i64 n) {
    if (n <= limit) return pref2[(int)n];
    auto it = cache2.find(n);
    if (it != cache2.end()) return it->second;
    i64 total = s3(n);
    for (i64 lo = 2; lo <= n;) {
      i64 q = n / lo;
      i64 hi = n / q;
      total = norm(total - range_s2(lo, hi) * i2_phi(q) % MOD);
      lo = hi + 1;
    }
    return cache2[n] = total;
  }
};

static void direction_stats(TotientSums &sums, i64 m, i64 &count, i64 &coord_sum, i64 &product_sum) {
  if (m <= 0) {
    count = coord_sum = product_sum = 0;
    return;
  }
  count = norm(2 * sums.phi_sum(m) - 1);
  coord_sum = norm(3 * sums.i_phi(m) - 1);
  product_sum = sums.i2_phi(m);
}

static i64 titanic_sets(i64 n, TotientSums &sums) {
  i64 side = n + 1;
  // side^2 overflows i64 for n=1e11; keep exponent in __int128.
  __int128 point_count = (__int128)side * side;
  i64 all_subsets = mod_pow(2, point_count);
  i64 singleton_part = norm(1 + (i64)(point_count % MOD));
  if (n < 2) return norm(all_subsets - singleton_part);

  struct Block {
    i64 lo, hi, a, b;
  };
  std::vector<Block> blocks;
  std::vector<i64> needed;
  for (i64 lo = 2; lo <= n;) {
    i64 a = n / lo;
    i64 b = lo < n ? n / (lo + 1) : 0;
    i64 hi_a = n / a;
    i64 hi_b = b ? n / b - 1 : n;
    i64 hi = std::min(hi_a, hi_b);
    blocks.push_back({lo, hi, a, b});
    needed.push_back(a);
    if (b) needed.push_back(b);
    lo = hi + 1;
  }
  std::sort(needed.begin(), needed.end());
  needed.erase(std::unique(needed.begin(), needed.end()), needed.end());
  std::unordered_map<i64, std::array<i64, 3>> stats;
  for (i64 m : needed) {
    i64 c, xy, pr;
    direction_stats(sums, m, c, xy, pr);
    stats[m] = {c, xy, pr};
  }
  i64 side_mod = side % MOD;
  i64 side2_mod = side_mod * side_mod % MOD;
  i64 collinear = 0;
  for (auto &bl : blocks) {
    auto s1s = stats[bl.a];
    i64 c1 = s1s[0], xy1 = s1s[1], pr1 = s1s[2];
    i64 c2 = 0, xy2 = 0, pr2 = 0;
    if (bl.b) {
      auto s2s = stats[bl.b];
      c2 = s2s[0];
      xy2 = s2s[1];
      pr2 = s2s[2];
    }
    i64 q2 = norm(pr1 - pr2);
    i64 q1 = norm(side_mod * norm(xy2 - xy1) % MOD - 2 * pr2 % MOD);
    i64 q0 = norm(side2_mod * norm(c1 - c2) % MOD + side_mod * xy2 % MOD - pr2);
    i64 p2 = 2 * q2 % MOD;
    i64 p1 = 2 * q1 % MOD;
    i64 p0 = (2 * q0 + 2 * side_mod) % MOD;
    i64 pow_lo = mod_pow(2, bl.lo);
    i64 pow_after_hi = mod_pow(2, bl.hi + 1);
    i64 e0 = norm(pow_after_hi - pow_lo);
    i64 e1 = norm(pref_k_pow2(bl.hi, pow_after_hi) - pref_k_pow2(bl.lo - 1, pow_lo));
    i64 e2 = norm(pref_k2_pow2(bl.hi, pow_after_hi) - pref_k2_pow2(bl.lo - 1, pow_lo));
    i64 poly_exp = (p2 * e2 + p1 * e1 + p0 * e0) % MOD;
    i64 r1 = range_s1(bl.lo, bl.hi);
    i64 r2 = range_s2(bl.lo, bl.hi);
    i64 r3 = range_s3(bl.lo, bl.hi);
    i64 length = (bl.hi - bl.lo + 1) % MOD;
    i64 poly_plain =
        (p2 * ((r3 + r2) % MOD) + p1 * ((r2 + r1) % MOD) + p0 * ((r1 + length) % MOD)) % MOD;
    collinear = norm(collinear + poly_exp - poly_plain);
  }
  return norm(all_subsets - singleton_part - collinear);
}

extern "C" long long pe415_answer(void) {
  TotientSums sums(std::max(DEFAULT_N, (i64)100000));
  return titanic_sets(DEFAULT_N, sums);
}
