// Project Euler 428 — Necklace of Circles T(10^9) via Min_25-style F,G.
#include <cstdint>
#include <cmath>
#include <unordered_map>
#include <vector>

using i64 = int64_t;

static int chi_mod3(i64 n) {
  int r = (int)(n % 3);
  if (r == 0) return 0;
  return r == 1 ? 1 : -1;
}

static i64 chi_prefix_integers(i64 x) {
  return (x + 2) / 3 - (x + 1) / 3;
}

static std::vector<int> primes_up_to(int limit) {
  if (limit < 2) return {};
  std::vector<char> sieve(limit + 1, 1);
  sieve[0] = sieve[1] = 0;
  int r = (int)std::sqrt((double)limit);
  for (int i = 2; i <= r; ++i)
    if (sieve[i])
      for (int j = i * i; j <= limit; j += i) sieve[j] = 0;
  std::vector<int> primes;
  for (int i = 2; i <= limit; ++i)
    if (sieve[i]) primes.push_back(i);
  return primes;
}

struct PrimeTables {
  i64 N;
  int V;
  std::vector<i64> vals;
  std::vector<int> id1, id2;
  std::vector<i64> pi_tbl, chi_tbl;
  std::vector<int> primes;

  explicit PrimeTables(i64 N_) : N(N_) {
    V = (int)std::sqrt((double)N);
    for (i64 i = 1; i <= N;) {
      i64 v = N / i;
      vals.push_back(v);
      i = N / v + 1;
    }
    id1.assign(V + 1, 0);
    id2.assign(V + 1, 0);
    for (size_t idx = 0; idx < vals.size(); ++idx) {
      i64 v = vals[idx];
      if (v <= V) id1[(int)v] = (int)idx;
      else id2[(int)(N / v)] = (int)idx;
    }
    size_t m = vals.size();
    pi_tbl.resize(m);
    chi_tbl.resize(m);
    for (size_t i = 0; i < m; ++i) {
      pi_tbl[i] = vals[i] - 1;
      chi_tbl[i] = chi_prefix_integers(vals[i]) - 1;
    }
    primes = primes_up_to(V);
    std::vector<i64> chi_prime_pref(1, 0);
    i64 s = 0;
    for (int p : primes) {
      s += chi_mod3(p);
      chi_prime_pref.push_back(s);
    }
    for (size_t pi_idx = 0; pi_idx < primes.size(); ++pi_idx) {
      int p = primes[pi_idx];
      i64 p2 = (i64)p * p;
      if (p2 > N) break;
      i64 pi_before = (i64)pi_idx;
      i64 chi_before = chi_prime_pref[pi_idx];
      int cp = chi_mod3(p);
      for (size_t idx = 0; idx < m; ++idx) {
        if (vals[idx] < p2) break;
        i64 t = vals[idx] / p;
        int j = idx_of(t);
        pi_tbl[idx] -= pi_tbl[j] - pi_before;
        chi_tbl[idx] -= (i64)cp * (chi_tbl[j] - chi_before);
      }
    }
  }

  int idx_of(i64 x) const {
    if (x <= V) return id1[(int)x];
    return id2[(int)(N / x)];
  }
  i64 pi(i64 x) const { return pi_tbl[idx_of(x)]; }
  i64 chi_prime_sum(i64 x) const { return chi_tbl[idx_of(x)]; }
};

struct Summatory {
  i64 N;
  PrimeTables pt;
  std::vector<int> rec_primes;
  using Cache = std::unordered_map<i64, i64>;
  // key = n * 100000 + idx (idx small)
  std::unordered_map<uint64_t, i64> cacheF, cacheG;

  explicit Summatory(i64 N_) : N(N_), pt(N_) {
    for (int p : pt.primes)
      if (p >= 5) rec_primes.push_back(p);
  }

  i64 prime_sum_F(i64 x) const {
    i64 pi = pt.pi(x);
    i64 c23 = (x >= 2) + (x >= 3);
    return 3 * (pi - c23);
  }
  i64 prime_sum_G(i64 x) const {
    i64 pi = pt.pi(x);
    i64 c23 = (x >= 2) + (x >= 3);
    i64 chi_p = pt.chi_prime_sum(x);
    i64 chi_excl = chi_p + (x >= 2 ? 1 : 0);
    return (pi - c23) + 2 * chi_excl;
  }

  static i64 val_F(int /*p*/, int e) { return 2 * e + 1; }
  static i64 val_G(int p, int e) {
    if (p % 3 == 1) return 2 * e + 1;
    return (e & 1) ? -1 : 1;
  }

  static uint64_t key(i64 n, int idx) { return ((uint64_t)n << 20) | (uint32_t)idx; }

  i64 H_F(i64 n, int idx) {
    auto it = cacheF.find(key(n, idx));
    if (it != cacheF.end()) return it->second;
    i64 prev = idx > 0 ? rec_primes[idx - 1] : 1;
    if (idx >= (int)rec_primes.size() || (i64)rec_primes[idx] * rec_primes[idx] > n) {
      i64 r = prime_sum_F(n) - prime_sum_F(prev);
      cacheF[key(n, idx)] = r;
      return r;
    }
    i64 res = prime_sum_F(n) - prime_sum_F(prev);
    for (int k = idx; k < (int)rec_primes.size(); ++k) {
      int p = rec_primes[k];
      if ((i64)p * p > n) break;
      i64 pe = p;
      int e = 1;
      while (pe * p <= n) {
        res += val_F(p, e) * H_F(n / pe, k + 1);
        res += val_F(p, e + 1);
        pe *= p;
        ++e;
      }
    }
    cacheF[key(n, idx)] = res;
    return res;
  }

  i64 H_G(i64 n, int idx) {
    auto it = cacheG.find(key(n, idx));
    if (it != cacheG.end()) return it->second;
    i64 prev = idx > 0 ? rec_primes[idx - 1] : 1;
    if (idx >= (int)rec_primes.size() || (i64)rec_primes[idx] * rec_primes[idx] > n) {
      i64 r = prime_sum_G(n) - prime_sum_G(prev);
      cacheG[key(n, idx)] = r;
      return r;
    }
    i64 res = prime_sum_G(n) - prime_sum_G(prev);
    for (int k = idx; k < (int)rec_primes.size(); ++k) {
      int p = rec_primes[k];
      if ((i64)p * p > n) break;
      i64 pe = p;
      int e = 1;
      while (pe * p <= n) {
        res += val_G(p, e) * H_G(n / pe, k + 1);
        res += val_G(p, e + 1);
        pe *= p;
        ++e;
      }
    }
    cacheG[key(n, idx)] = res;
    return res;
  }

  i64 sum_F(i64 x) {
    if (x <= 0) return 0;
    return 1 + H_F(x, 0);
  }
  i64 sum_G(i64 x) {
    if (x <= 0) return 0;
    return 1 + H_G(x, 0);
  }
};

static i64 T(i64 n) {
  Summatory S(n);
  std::unordered_map<i64, i64> F_cache, G_cache;
  auto F = [&](i64 x) {
    auto it = F_cache.find(x);
    if (it != F_cache.end()) return it->second;
    i64 v = S.sum_F(x);
    F_cache[x] = v;
    return v;
  };
  auto G = [&](i64 x) {
    auto it = G_cache.find(x);
    if (it != G_cache.end()) return it->second;
    i64 v = S.sum_G(x);
    G_cache[x] = v;
    return v;
  };

  i64 total = 0;
  i64 pow2 = 1;
  for (int i = 0; pow2 <= n; ++i, pow2 *= 2) {
    i64 pow3 = 1;
    for (int j = 0; pow2 * pow3 <= n; ++j, pow3 *= 3) {
      i64 x = n / (pow2 * pow3);
      i64 Fx = F(x);
      total += (2 * i + 2) * (2 * j + 1) * Fx;
      total += (2 * i + 3) * (2 * j + 2) * Fx;
      if (j >= 1) total += (2 * j - 1) * (2 * i + 3) * Fx;
    }
  }
  pow2 = 1;
  for (int i = 0; pow2 <= n; ++i, pow2 *= 2) {
    i64 x = n / pow2;
    i64 sign = (i % 2 == 0) ? 1 : -1;
    total += ((2 * i + 3) * F(x) - sign * G(x)) / 2;
  }
  return total;
}

extern "C" long long pe428_answer(void) { return T(1000000000LL); }
