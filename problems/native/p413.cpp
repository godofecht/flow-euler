// Project Euler 413 — One-child numbers F(10^19).
#include <algorithm>
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

using i64 = int64_t;

static int v_factor(int n, int p) {
  int cnt = 0;
  while (n % p == 0) {
    n /= p;
    ++cnt;
  }
  return cnt;
}

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

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

static void precompute_counts(int m, std::vector<std::vector<int>> &count_val,
                              std::vector<std::vector<int>> &inc_code) {
  if (m <= 0) m = 1;
  int total = 1;
  for (int i = 0; i < m; ++i) total *= 3;
  std::vector<int> pow3(m + 1, 1);
  for (int i = 1; i <= m; ++i) pow3[i] = pow3[i - 1] * 3;
  count_val.assign(total, std::vector<int>(m));
  inc_code.assign(total, std::vector<int>(m));
  for (int code = 0; code < total; ++code) {
    int tmp = code;
    std::vector<int> counts(m);
    for (int i = 0; i < m; ++i) {
      counts[i] = tmp % 3;
      tmp /= 3;
    }
    for (int i = 0; i < m; ++i) {
      count_val[code][i] = counts[i];
      int neu = counts[i] + 1;
      if (neu > 2) neu = 2;
      inc_code[code][i] = code + (neu - counts[i]) * pow3[i];
    }
  }
}

static i64 count_one_child_coprime(int d) {
  i64 inv10 = mod_inv(10, d);
  std::vector<i64> w(d + 1);
  w[1] = inv10 % d;
  for (int i = 2; i <= d; ++i) w[i] = w[i - 1] * inv10 % d;

  std::vector<std::vector<std::pair<int, int>>> delta_counts(d);
  for (int i = 1; i <= d; ++i) {
    std::map<int, int> counts;
    i64 wi = w[i];
    int lo = (i == 1) ? 1 : 0;
    for (int digit = lo; digit <= 9; ++digit) {
      int delta = (int)((digit * wi) % d);
      counts[delta] += 1;
    }
    for (auto &kv : counts) delta_counts[i - 1].push_back(kv);
  }

  // key = (mask << 5) | r
  std::map<i64, i64> dp0, dp1;
  dp0[(1LL << 0) << 5 | 0] = 1;
  for (int i = 0; i < d; ++i) {
    std::map<i64, i64> next0, next1;
    for (auto &kv : dp0) {
      i64 key = kv.first;
      i64 ways = kv.second;
      i64 mask = key >> 5;
      int r = (int)(key & 31);
      for (auto &dc : delta_counts[i]) {
        int r2 = r + dc.first;
        if (r2 >= d) r2 -= d;
        i64 bit = 1LL << r2;
        if (mask & bit) {
          i64 key2 = (mask << 5) | r2;
          next1[key2] += ways * dc.second;
        } else {
          i64 key2 = ((mask | bit) << 5) | r2;
          next0[key2] += ways * dc.second;
        }
      }
    }
    for (auto &kv : dp1) {
      i64 key = kv.first;
      i64 ways = kv.second;
      i64 mask = key >> 5;
      int r = (int)(key & 31);
      for (auto &dc : delta_counts[i]) {
        int r2 = r + dc.first;
        if (r2 >= d) r2 -= d;
        i64 bit = 1LL << r2;
        if (mask & bit) continue;
        i64 key2 = ((mask | bit) << 5) | r2;
        next1[key2] += ways * dc.second;
      }
    }
    dp0.swap(next0);
    dp1.swap(next1);
  }
  i64 total = 0;
  for (auto &kv : dp1) total += kv.second;
  return total;
}

static i64 count_one_child_with_t(int d) {
  int a = v_factor(d, 2);
  int b = v_factor(d, 5);
  int t = 1;
  for (int i = 0; i < a; ++i) t *= 2;
  for (int i = 0; i < b; ++i) t *= 5;
  int m = d / t;
  int L = std::max(a, b);

  std::vector<std::vector<int>> count_val, inc_code;
  precompute_counts(m > 0 ? m : 1, count_val, inc_code);

  std::vector<i64> w(d + 1, 0);
  if (m > 1) {
    i64 inv10 = mod_inv(10, m);
    w[1] = inv10 % m;
    for (int i = 2; i <= d; ++i) w[i] = w[i - 1] * inv10 % m;
  }

  if (L == 1) {
    std::map<std::pair<int, int>, i64> dp[3];
    int start_code = inc_code[0][0];
    dp[0][{start_code, 0}] = 1;
    for (int pos = 1; pos <= d; ++pos) {
      std::map<std::pair<int, int>, i64> next_dp[3];
      int lo = (pos == 1) ? 1 : 0;
      i64 wi = w[pos];
      for (int div_count = 0; div_count < 3; ++div_count) {
        for (auto &kv : dp[div_count]) {
          int counts_code = kv.first.first;
          int q = kv.first.second;
          i64 ways = kv.second;
          for (int digit = lo; digit <= 9; ++digit) {
            int q_next = m > 1 ? (int)((q + digit * wi) % m) : 0;
            int new_div = 0;
            if (digit % t == 0) new_div = count_val[counts_code][q_next];
            int new_total = div_count + new_div;
            if (new_total > 2) new_total = 2;
            int new_counts_code = inc_code[counts_code][q_next];
            next_dp[new_total][{new_counts_code, q_next}] += ways;
          }
        }
      }
      for (int i = 0; i < 3; ++i) dp[i].swap(next_dp[i]);
    }
    i64 total = 0;
    for (auto &kv : dp[1]) total += kv.second;
    return total;
  }

  int mod10L = 1;
  for (int i = 0; i < L; ++i) mod10L *= 10;
  std::vector<std::vector<int>> next_tail(mod10L, std::vector<int>(10));
  std::vector<int> short_mask(mod10L), pow10(L + 1, 1);
  std::vector<char> long_ok(mod10L);
  for (int i = 1; i <= L; ++i) pow10[i] = pow10[i - 1] * 10;
  for (int tail = 0; tail < mod10L; ++tail) {
    long_ok[tail] = (tail % t == 0);
    int mask = 0;
    for (int l = 1; l < L; ++l)
      if ((tail % pow10[l]) % t == 0) mask |= 1 << (l - 1);
    short_mask[tail] = mask;
    int base = (tail * 10) % mod10L;
    for (int digit = 0; digit < 10; ++digit) next_tail[tail][digit] = (base + digit) % mod10L;
  }

  // key: (counts_code, recent_packed, recent_len, tail)
  // recent has at most L-1 entries each < m (<=19), pack in base 32
  using Key = std::tuple<int, i64, int, int>;
  std::map<Key, i64> dp[3];
  // recent starts as (0,), matching the Python reference.
  dp[0][Key{0, 0, 1, 0}] = 1;

  auto pack_push = [&](i64 packed, int len, int q) -> std::pair<i64, int> {
    packed = packed * 32 + q;
    ++len;
    return {packed, len};
  };
  auto pack_get = [&](i64 packed, int /*len*/, int from_end) -> int {
    i64 p = packed;
    for (int i = 0; i < from_end - 1; ++i) p /= 32;
    return (int)(p % 32);
  };
  auto pack_pop_front = [&](i64 packed, int len) -> std::tuple<int, i64, int> {
    // return (oldest, new_packed, new_len)
    i64 div = 1;
    for (int i = 1; i < len; ++i) div *= 32;
    int oldest = (int)(packed / div);
    i64 neu = packed % div;
    return {oldest, neu, len - 1};
  };

  for (int pos = 1; pos <= d; ++pos) {
    std::map<Key, i64> next_dp[3];
    int lo = (pos == 1) ? 1 : 0;
    i64 wi = w[pos];
    for (int div_count = 0; div_count < 3; ++div_count) {
      for (auto &kv : dp[div_count]) {
        int counts_code = std::get<0>(kv.first);
        i64 recent = std::get<1>(kv.first);
        int rlen = std::get<2>(kv.first);
        int tail = std::get<3>(kv.first);
        i64 ways = kv.second;
        int current_q = rlen ? pack_get(recent, rlen, 1) : 0;
        int max_l = L - 1;
        if (pos < max_l) max_l = pos;
        for (int digit = lo; digit <= 9; ++digit) {
          int q_next = m > 1 ? (int)((current_q + digit * wi) % m) : 0;
          int new_tail = next_tail[tail][digit];
          int mask = short_mask[new_tail];
          int new_div = 0;
          if (mask) {
            for (int l = 1; l <= max_l; ++l) {
              if (mask & (1 << (l - 1))) {
                if (pack_get(recent, rlen, l) == q_next) ++new_div;
              }
            }
          }
          if (pos >= L && long_ok[new_tail]) new_div += count_val[counts_code][q_next];
          if (new_div > 2) new_div = 2;
          int new_total = div_count + new_div;
          if (new_total > 2) new_total = 2;
          auto pr = pack_push(recent, rlen, q_next);
          i64 new_recent = pr.first;
          int new_rlen = pr.second;
          int new_counts_code = counts_code;
          if (new_rlen > L - 1) {
            auto po = pack_pop_front(new_recent, new_rlen);
            int oldest = std::get<0>(po);
            new_recent = std::get<1>(po);
            new_rlen = std::get<2>(po);
            new_counts_code = inc_code[new_counts_code][oldest];
          }
          next_dp[new_total][Key{new_counts_code, new_recent, new_rlen, new_tail}] += ways;
        }
      }
    }
    for (int i = 0; i < 3; ++i) dp[i].swap(next_dp[i]);
  }
  i64 total = 0;
  for (auto &kv : dp[1]) total += kv.second;
  return total;
}

static i64 count_one_child(int d) {
  if (d == 1) return 9;
  if (d % 2 != 0 && d % 5 != 0) return count_one_child_coprime(d);
  return count_one_child_with_t(d);
}

extern "C" long long pe413_answer(void) {
  i64 total = 0;
  for (int d = 1; d < 20; ++d) total += count_one_child(d);
  return total;
}
