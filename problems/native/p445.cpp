// Project Euler 445 — Retractions A, S(10^7) mod 10^9+7.
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using u32 = uint32_t;
using i64 = int64_t;

static constexpr i64 MOD = 1000000007LL;

static i64 mod_pow(i64 a, i64 e) {
  i64 r = 1;
  a %= MOD;
  while (e > 0) {
    if (e & 1) r = r * a % MOD;
    a = a * a % MOD;
    e >>= 1;
  }
  return r;
}

static void sieve_spf(int n, std::vector<u32> &spf, std::vector<u32> &primes,
                      std::vector<u32> &prime_idx) {
  spf.assign(n + 1, 0);
  prime_idx.assign(n + 1, 0);
  primes.clear();
  if (n >= 1) spf[0] = spf[1] = 1;
  if (n >= 2) {
    spf[2] = 2;
    primes.push_back(2);
    prime_idx[2] = 1;
    for (int x = 4; x <= n; x += 2) spf[x] = 2;
  }
  int limit = (int)std::sqrt((double)n);
  for (int i = 3; i <= n; i += 2) {
    if (spf[i] == 0) {
      spf[i] = i;
      primes.push_back(i);
      prime_idx[i] = (u32)primes.size();
      if (i <= limit) {
        int step = i << 1;
        for (i64 j = (i64)i * i; j <= n; j += step)
          if (spf[(int)j] == 0) spf[(int)j] = i;
      }
    }
  }
}

static std::vector<u32> inverses_upto(int n) {
  std::vector<u32> inv(n + 1);
  inv[1] = 1;
  for (int i = 2; i <= n; ++i)
    inv[i] = (u32)(MOD - (MOD / i) * inv[MOD % i] % MOD);
  return inv;
}

static std::vector<u32> batch_inverse(const std::vector<u32> &vals) {
  int n = (int)vals.size();
  std::vector<u32> invs(n, 0);
  std::vector<int> idxs;
  std::vector<u32> prefix;
  i64 prod = 1;
  for (int i = 0; i < n; ++i) {
    if (vals[i]) {
      prod = prod * vals[i] % MOD;
      idxs.push_back(i);
      prefix.push_back((u32)prod);
    }
  }
  if (idxs.empty()) return invs;
  i64 inv_all = mod_pow(prod, MOD - 2);
  for (int j = (int)idxs.size() - 1; j >= 0; --j) {
    int i = idxs[j];
    i64 prev = j ? prefix[j - 1] : 1;
    invs[i] = (u32)(inv_all * prev % MOD);
    inv_all = inv_all * vals[i] % MOD;
  }
  return invs;
}

static i64 solve(int N) {
  std::vector<u32> spf, primes, prime_idx;
  sieve_spf(N, spf, primes, prime_idx);
  auto inv_num = inverses_upto(N);
  int num_primes = (int)primes.size();
  std::vector<u32> max_exp(num_primes), offset(num_primes);
  int total_terms = 0;
  for (int idx = 0; idx < num_primes; ++idx) {
    int p = (int)primes[idx];
    int t = N, e = 0;
    while (t) {
      t /= p;
      e += t;
    }
    max_exp[idx] = e;
    offset[idx] = total_terms;
    total_terms += e;
  }
  std::vector<u32> inv_terms(total_terms);
  const int CHUNK = 1000000;
  int write_pos = 0;
  std::vector<u32> vals;
  vals.reserve(CHUNK);
  for (int idx = 0; idx < num_primes; ++idx) {
    i64 p = primes[idx];
    int m = (int)max_exp[idx];
    i64 pow_val = p % MOD;
    for (int k = 0; k < m; ++k) {
      vals.push_back((u32)((pow_val + 1) % MOD));
      pow_val = pow_val * p % MOD;
      if ((int)vals.size() >= CHUNK) {
        auto invs = batch_inverse(vals);
        for (size_t i = 0; i < vals.size(); ++i) inv_terms[write_pos++] = invs[i];
        vals.clear();
      }
    }
  }
  if (!vals.empty()) {
    auto invs = batch_inverse(vals);
    for (size_t i = 0; i < vals.size(); ++i) inv_terms[write_pos++] = invs[i];
  }

  std::vector<u32> exp(num_primes, 0), p_pow(num_primes, 1);
  i64 prod = 1;
  int zero_count = 0;
  int mid = N / 2;
  bool even = (N % 2 == 0);
  i64 sum_sigma = 0;

  for (int k = 1; k <= mid; ++k) {
    int numer = N - k + 1;
    int denom = k;

    auto apply_factor = [&](int x, int sign) {
      while (x > 1) {
        int p = (int)spf[x];
        int pi = (int)prime_idx[p] - 1;
        int cnt = 0;
        while (x > 1 && (int)spf[x] == p) {
          x /= p;
          ++cnt;
        }
        int old_e = (int)exp[pi];
        if (sign > 0) {
          if (old_e) {
            i64 term_old = (i64)p_pow[pi] + 1;
            if (term_old == MOD) term_old = 0;
            if (term_old) prod = prod * inv_terms[offset[pi] + old_e - 1] % MOD;
            else --zero_count;
          }
          int new_e = old_e + cnt;
          exp[pi] = new_e;
          if (cnt == 1) p_pow[pi] = (u32)((i64)p_pow[pi] * p % MOD);
          else if (cnt == 2) p_pow[pi] = (u32)((i64)p_pow[pi] * p % MOD * p % MOD);
          else p_pow[pi] = (u32)((i64)p_pow[pi] * mod_pow(p, cnt) % MOD);
          i64 term_new = (i64)p_pow[pi] + 1;
          if (term_new == MOD) term_new = 0;
          if (term_new) prod = prod * term_new % MOD;
          else ++zero_count;
        } else {
          i64 term_old = (i64)p_pow[pi] + 1;
          if (term_old == MOD) term_old = 0;
          if (term_old) prod = prod * inv_terms[offset[pi] + old_e - 1] % MOD;
          else --zero_count;
          int new_e = old_e - cnt;
          exp[pi] = new_e;
          i64 invp = inv_num[p];
          if (cnt == 1) p_pow[pi] = (u32)((i64)p_pow[pi] * invp % MOD);
          else if (cnt == 2) p_pow[pi] = (u32)((i64)p_pow[pi] * invp % MOD * invp % MOD);
          else p_pow[pi] = (u32)((i64)p_pow[pi] * mod_pow(invp, cnt) % MOD);
          if (new_e) {
            i64 term_new = (i64)p_pow[pi] + 1;
            if (term_new == MOD) term_new = 0;
            if (term_new) prod = prod * term_new % MOD;
            else ++zero_count;
          }
        }
      }
    };

    apply_factor(numer, +1);
    apply_factor(denom, -1);

    i64 sigma_val = zero_count ? 0 : prod;
    if (even && k == mid) sum_sigma += sigma_val;
    else sum_sigma += 2 * sigma_val;
    if (sum_sigma >= (1LL << 62)) sum_sigma %= MOD;
  }
  sum_sigma %= MOD;
  i64 sum_binom = (mod_pow(2, N) - 2 + MOD) % MOD;
  return (sum_sigma - sum_binom + MOD) % MOD;
}

extern "C" long long pe445_answer(void) { return solve(10000000); }
