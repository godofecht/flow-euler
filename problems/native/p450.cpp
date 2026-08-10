// Project Euler 450 — Hypocycloid lattice points T(10^6).
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <unordered_map>
#include <vector>

using i64 = int64_t;

static std::vector<int> mobius_sieve(int n) {
  std::vector<int> mu(n + 1), primes;
  std::vector<char> comp(n + 1, 0);
  mu[1] = 1;
  for (int i = 2; i <= n; ++i) {
    if (!comp[i]) {
      primes.push_back(i);
      mu[i] = -1;
    }
    for (int p : primes) {
      i64 ip = (i64)i * p;
      if (ip > n) break;
      comp[(int)ip] = 1;
      if (i % p == 0) {
        mu[(int)ip] = 0;
        break;
      }
      mu[(int)ip] = -mu[i];
    }
  }
  return mu;
}

static void build_prefix(const std::vector<int> &mu, std::vector<i64> &pref_mui,
                         std::vector<i64> &pref_mui_odd) {
  int n = (int)mu.size() - 1;
  pref_mui.assign(n + 1, 0);
  pref_mui_odd.assign(n + 1, 0);
  i64 s1 = 0, so1 = 0;
  for (int i = 1; i <= n; ++i) {
    s1 += (i64)mu[i] * i;
    if (i & 1) so1 += (i64)mu[i] * i;
    pref_mui[i] = s1;
    pref_mui_odd[i] = so1;
  }
}

static i64 G0(i64 n) {
  if (n <= 2) return 0;
  i64 m = (n - 1) / 2;
  return m * (n - m - 1);
}
static i64 GA(i64 n) { return (n + 1) * G0(n) / 2; }
static i64 GB(i64 n) {
  if (n <= 2) return 0;
  i64 m = (n - 1) / 2;
  i64 term = m * (m + 1);
  return term * (3 * n - 4 * m - 2) / 6;
}
static i64 HB(i64 n) {
  i64 M = (n - 2) / 2;
  if (M <= 0) return 0;
  i64 m0 = (M - 1) / 2;
  if (m0 < 0) return 0;
  i64 s_y = m0 * (m0 + 1) / 2;
  i64 s_y2 = m0 * (m0 + 1) * (2 * m0 + 1) / 6;
  return M * (m0 + 1) * (m0 + 1) - 4 * s_y2 - 2 * s_y;
}
static i64 KB(i64 n) {
  i64 res = 0;
  i64 M1 = (n - 2) / 4;
  if (M1 > 0) {
    i64 m1 = (M1 - 1) / 2;
    if (m1 >= 0) {
      i64 s_y = m1 * (m1 + 1) / 2;
      i64 s_y2 = m1 * (m1 + 1) * (2 * m1 + 1) / 6;
      i64 sum_4y1 = (m1 + 1) * (2 * m1 + 1);
      res += M1 * sum_4y1 - 8 * s_y2 - 2 * s_y;
    }
  }
  i64 M3 = (n - 6) / 4;
  if (M3 > 0) {
    i64 m3 = (M3 - 1) / 2;
    if (m3 >= 0) {
      i64 s_y = m3 * (m3 + 1) / 2;
      i64 s_y2 = m3 * (m3 + 1) * (2 * m3 + 1) / 6;
      i64 sum_4y3 = (m3 + 1) * (2 * m3 + 3);
      res += M3 * sum_4y3 - 8 * s_y2 - 6 * s_y;
    }
  }
  return res;
}

template <typename F>
static i64 mobius_sum_weighted(i64 m, const std::vector<i64> &pref_mui, F base) {
  i64 res = 0, l = 1;
  while (l <= m) {
    i64 q = m / l;
    i64 r = m / q;
    res += (pref_mui[(size_t)r] - pref_mui[(size_t)l - 1]) * base(q);
    l = r + 1;
  }
  return res;
}

struct AxisCalculator {
  std::vector<i64> pref_mui, pref_mui_odd;
  std::unordered_map<i64, i64> cache;
  explicit AxisCalculator(int n) {
    auto mu = mobius_sieve(n);
    build_prefix(mu, pref_mui, pref_mui_odd);
  }
  i64 f(i64 m) {
    auto it = cache.find(m);
    if (it != cache.end()) return it->second;
    i64 SA = mobius_sum_weighted(m, pref_mui, GA);
    i64 SB = mobius_sum_weighted(m, pref_mui, GB);
    i64 P2 = mobius_sum_weighted(m, pref_mui_odd, HB);
    i64 P4 = mobius_sum_weighted(m, pref_mui_odd, KB);
    i64 val = 4 * SA + 2 * SB + 2 * P2 - 4 * P4;
    cache[m] = val;
    return val;
  }
  i64 total(i64 n) {
    i64 tot = 0, l = 1;
    while (l <= n) {
      i64 q = n / l;
      i64 r = n / q;
      i64 sum_d = (l + r) * (r - l + 1) / 2;
      tot += sum_d * f(q);
      l = r + 1;
    }
    return tot;
  }
};

static i64 igcd(i64 a, i64 b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  while (b) {
    i64 t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static i64 ipow(i64 b, int e) {
  i64 r = 1;
  while (e > 0) {
    if (e & 1) r *= b;
    b *= b;
    e >>= 1;
  }
  return r;
}

struct Triple {
  i64 a, b, c;
};

static std::vector<Triple> primitive_triples(i64 cmax) {
  std::vector<Triple> triples;
  i64 m_limit = (i64)std::sqrt((double)cmax * 2.0) + 3;
  for (i64 m = 2; m <= m_limit; ++m) {
    i64 mm = m * m;
    for (i64 n = 1; n < m; ++n) {
      if (((m - n) & 1) == 0) continue;
      if (igcd(m, n) != 1) continue;
      i64 c = mm + n * n;
      if (c > cmax) break;
      triples.push_back({mm - n * n, 2 * m * n, c});
    }
  }
  return triples;
}

static void gauss_pow(i64 re, i64 im, int exp, i64 &rr, i64 &ri) {
  rr = 1;
  ri = 0;
  i64 br = re, bi = im;
  int e = exp;
  while (e > 0) {
    if (e & 1) {
      i64 nr = rr * br - ri * bi;
      i64 ni = rr * bi + ri * br;
      rr = nr;
      ri = ni;
    }
    i64 nbr = br * br - bi * bi;
    i64 nbi = br * bi + bi * br;
    br = nbr;
    bi = nbi;
    e >>= 1;
  }
}

static i64 non_axis_total(i64 n) {
  i64 total = 0;
  int maxA = 1;
  while (ipow(3, maxA + 1) <= n) ++maxA;

  for (int Ap = 2; Ap <= maxA; ++Ap) {
    i64 cmax = (i64)std::pow((double)n, 1.0 / Ap) + 2;
    while (ipow(cmax, Ap) > n) --cmax;
    auto triples = primitive_triples(cmax);

    for (int Bp = 1; Bp < Ap; ++Bp) {
      if (std::gcd(Ap, Bp) != 1) continue;
      for (auto &t : triples) {
        i64 den = ipow(t.c, Ap);
        if (den * (Ap + Bp) > n) continue;

        i64 variants[16][2];
        int nv = 0;
        i64 pairs[2][2] = {{t.a, t.b}, {t.b, t.a}};
        for (int pi = 0; pi < 2; ++pi) {
          for (int su : {1, -1}) {
            for (int sv : {1, -1}) {
              i64 u = su * pairs[pi][0], v = sv * pairs[pi][1];
              bool dup = false;
              for (int k = 0; k < nv; ++k)
                if (variants[k][0] == u && variants[k][1] == v) {
                  dup = true;
                  break;
                }
              if (!dup) {
                variants[nv][0] = u;
                variants[nv][1] = v;
                ++nv;
              }
            }
          }
        }

        for (int vi = 0; vi < nv; ++vi) {
          i64 re = variants[vi][0], im = variants[vi][1];
          i64 uA, vA, uB, vB;
          gauss_pow(re, im, Ap, uA, vA);
          gauss_pow(re, im, Bp, uB, vB);
          i64 scale = ipow(t.c, Ap - Bp);
          i64 numX = Ap * uB * scale + Bp * uA;
          i64 numY = Ap * vB * scale - Bp * vA;
          i64 g = igcd(den, igcd(std::llabs(numX), std::llabs(numY)));
          i64 d0 = den / g;
          if (d0 * (Ap + Bp) > n) continue;
          i64 x0 = numX / g, y0 = numY / g;
          i64 kmax = n / (d0 * (Ap + Bp));
          total += (std::llabs(x0) + std::llabs(y0)) * (kmax * (kmax + 1) / 2);
        }
      }
    }
  }
  return total;
}

extern "C" long long pe450_answer(void) {
  const i64 n = 1000000;
  AxisCalculator axis((int)n);
  return axis.total(n) + non_axis_total(n);
}
