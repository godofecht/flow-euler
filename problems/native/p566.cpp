#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <thread>
#include <unordered_set>
#include <vector>


static inline long long llgcd(long long a, long long b) {
  while (b != 0) {
    long long t = a % b;
    a = b;
    b = t;
  }
  return a < 0 ? -a : a;
}

static inline long long lllcm(long long a, long long b) {
  return a / llgcd(a, b) * b;
}

struct Point {
  long long P = 0;
  long long Q = 0;
  bool operator==(const Point &o) const { return P == o.P && Q == o.Q; }
};

struct PointHash {
  size_t operator()(const Point &p) const noexcept {
    uint64_t x = static_cast<uint64_t>(p.P) * 0x9e3779b97f4a7c15ULL;
    uint64_t y = static_cast<uint64_t>(p.Q) + 0x9e3779b97f4a7c15ULL;
    x ^= y + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    return static_cast<size_t>(x);
  }
};

static inline int sign_AB_sqrtc(__int128 A, __int128 B, int c) {
  if (B == 0) {
    if (A == 0) return 0;
    return (A > 0) ? 1 : -1;
  }
  if (A == 0) return (B > 0) ? 1 : -1;
  if ((A > 0 && B > 0) || (A < 0 && B < 0)) {
    return (A > 0) ? 1 : -1;
  }
  __int128 A2 = A * A;
  __int128 B2c = B * B * static_cast<__int128>(c);
  if (A < 0 && B > 0) {
    if (B2c == A2) return 0;
    return (B2c > A2) ? 1 : -1;
  }
  // A > 0, B < 0
  if (B2c == A2) return 0;
  return (B2c > A2) ? -1 : 1;
}

static inline int cmpPoints(const Point &a, const Point &b, long long ab, int c) {
  long long dP = a.P - b.P;
  long long dQ = a.Q - b.Q;
  __int128 A = static_cast<__int128>(dP) * c;
  __int128 B = static_cast<__int128>(dQ) * ab;
  return sign_AB_sqrtc(A, B, c);
}

struct FieldContext {
  int a = 0;
  int b = 0;
  int c = 0;
  long long ab = 0;
  long double sqrtc = 0.0L;
};

static inline Point canonical01(long long P, long long Q, const FieldContext &ctx) {
  long double approx = static_cast<long double>(P) / static_cast<long double>(ctx.ab) +
                       static_cast<long double>(Q) * ctx.sqrtc / static_cast<long double>(ctx.c);
  long long n = static_cast<long long>(std::floor(approx));
  P -= n * ctx.ab;

  auto less0 = [&](long long PP, long long QQ) {
    __int128 A = static_cast<__int128>(PP) * ctx.c;
    __int128 B = static_cast<__int128>(QQ) * ctx.ab;
    return sign_AB_sqrtc(A, B, ctx.c) < 0;
  };
  auto ge1 = [&](long long PP, long long QQ) {
    __int128 A = static_cast<__int128>(PP - ctx.ab) * ctx.c;
    __int128 B = static_cast<__int128>(QQ) * ctx.ab;
    return sign_AB_sqrtc(A, B, ctx.c) >= 0;
  };

  while (less0(P, Q)) P += ctx.ab;
  while (ge1(P, Q)) P -= ctx.ab;

  return {P, Q};
}

static inline Point addP(const Point &x, const Point &s, const FieldContext &ctx) {
  return canonical01(x.P + s.P, x.Q + s.Q, ctx);
}

static inline Point subP(const Point &x, const Point &s, const FieldContext &ctx) {
  return canonical01(x.P - s.P, x.Q - s.Q, ctx);
}

static inline Point oneMinus(const Point &x, const FieldContext &ctx) {
  return canonical01(ctx.ab - x.P, -x.Q, ctx);
}

static std::vector<int> kmpPrefix(const std::vector<uint8_t> &pat) {
  int n = static_cast<int>(pat.size());
  std::vector<int> pi(n, 0);
  for (int i = 1, j = 0; i < n; ++i) {
    while (j > 0 && pat[i] != pat[j]) j = pi[j - 1];
    if (pat[i] == pat[j]) ++j;
    pi[i] = j;
  }
  return pi;
}

static std::vector<int> findRotations(const std::vector<uint8_t> &w,
                                      const std::vector<uint8_t> &e) {
  int n = static_cast<int>(w.size());
  std::vector<uint8_t> text(2 * n);
  for (int i = 0; i < n; ++i) {
    text[i] = w[i];
    text[i + n] = w[i];
  }
  auto pi = kmpPrefix(e);
  std::vector<int> shifts;
  for (int i = 0, j = 0; i < 2 * n - 1; ++i) {
    while (j > 0 && text[i] != e[j]) j = pi[j - 1];
    if (text[i] == e[j]) ++j;
    if (j == n) {
      int pos = i - n + 1;
      if (pos < n) {
        int shift = (n - pos) % n;
        shifts.push_back(shift);
      }
      j = pi[j - 1];
    }
  }
  return shifts;
}

static inline long long egcd(long long a, long long b, long long &x, long long &y) {
  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  }
  long long x1, y1;
  long long g = egcd(b, a % b, x1, y1);
  x = y1;
  y = x1 - (a / b) * y1;
  return g;
}

static inline bool crtPair(long long a1, long long m1, long long a2, long long m2,
                           long long &a, long long &m) {
  long long g = llgcd(m1, m2);
  long long diff = a2 - a1;
  if (diff % g != 0) return false;
  long long m1g = m1 / g;
  long long m2g = m2 / g;

  long long x, y;
  egcd(m1g, m2g, x, y);
  x %= m2g;
  if (x < 0) x += m2g;

  __int128 t = static_cast<__int128>(diff / g) * x % m2g;
  if (t < 0) t += m2g;

  __int128 l = static_cast<__int128>(m1) / g * m2;
  __int128 res = static_cast<__int128>(a1) + static_cast<__int128>(m1) * t;
  long long mod = static_cast<long long>(l);
  long long ans = static_cast<long long>(res % l);
  if (ans < 0) ans += mod;
  a = ans;
  m = mod;
  return true;
}

static void mergeCRT(std::vector<long long> &sols, long long &M,
                     const std::vector<long long> &residues, long long P) {
  std::vector<long long> next;
  next.reserve(sols.size() * residues.size());
  long long newM = lllcm(M, P);
  for (long long s : sols) {
    for (long long r : residues) {
      long long a, m;
      if (crtPair(s, M, r, P, a, m)) {
        next.push_back(a);
      }
    }
  }
  std::sort(next.begin(), next.end());
  next.erase(std::unique(next.begin(), next.end()), next.end());
  sols.swap(next);
  M = newM;
}

static std::pair<long long, std::vector<long long>>
cycleConstraint(const std::vector<uint8_t> &f) {
  int L = static_cast<int>(f.size());
  int n = L / 3;

  uint8_t total_f = 0;
  for (uint8_t b : f) total_f ^= b;
  long long P = (total_f == 0) ? L : 2LL * L;

  std::vector<uint8_t> a(n), b(n), c(n);
  for (int i = 0; i < n; ++i) {
    a[i] = f[3 * i];
    b[i] = f[3 * i + 1];
    c[i] = f[3 * i + 2];
  }

  std::vector<uint8_t> w(n);
  uint8_t total_w = 0;
  for (int i = 0; i < n; ++i) {
    w[i] = a[i] ^ b[i] ^ c[i];
    total_w ^= w[i];
  }

  std::vector<uint8_t> pref(n + 1, 0);
  for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] ^ w[i];
  auto H0 = [&](int s) { return pref[s]; };

  std::vector<long long> good;
  for (int R = 0; R < 3; ++R) {
    std::vector<uint8_t> t(n);
    if (R == 0) {
      std::fill(t.begin(), t.end(), 0);
    } else if (R == 1) {
      t = a;
    } else {
      for (int i = 0; i < n; ++i) t[i] = a[i] ^ b[i];
    }

    std::vector<uint8_t> dt(n);
    for (int i = 0; i < n; ++i) dt[i] = t[i] ^ t[(i + 1) % n];
    std::vector<uint8_t> e(n);
    for (int i = 0; i < n; ++i) e[i] = w[i] ^ dt[i];

    auto shifts = findRotations(w, e);
    for (int sShift : shifts) {
      for (int q = 0; q < 2; ++q) {
        int mTrip = sShift + q * n;
        uint8_t h0 = H0(sShift) ^ static_cast<uint8_t>(q * total_w);
        uint8_t target = t[mTrip % n];
        if (h0 != target) continue;
        long long N = 3LL * mTrip + R;
        long long Nmod = N % P;
        good.push_back(Nmod);
      }
    }
  }

  std::sort(good.begin(), good.end());
  good.erase(std::unique(good.begin(), good.end()), good.end());
  return {P, good};
}

static long long F_square(int a, int b, int c) {
  int k = static_cast<int>(std::llround(std::sqrt(static_cast<long double>(c))));
  long long D = 1;
  D = lllcm(D, a);
  D = lllcm(D, b);
  D = lllcm(D, k);

  long long S[3] = {D / a, D / b, D / k};
  int m = static_cast<int>(D);
  int N = 3 * m;
  std::vector<int> nxt(N);
  std::vector<uint8_t> flip(N);

  for (int ph = 0; ph < 3; ++ph) {
    long long Sph = S[ph];
    for (int j = 0; j < m; ++j) {
      int j2;
      uint8_t fbit;
      if (j < Sph) {
        j2 = m - 1 - j;
        fbit = 1;
      } else {
        j2 = j - static_cast<int>(Sph);
        fbit = 0;
      }
      int cur = ph * m + j;
      int nxtph = (ph + 1) % 3;
      nxt[cur] = nxtph * m + j2;
      flip[cur] = fbit;
    }
  }

  std::vector<char> vis(N, 0);
  std::vector<std::pair<long long, std::vector<long long>>> constraints;

  for (int s = 0; s < N; ++s) {
    if (vis[s]) continue;
    std::vector<int> cyc;
    int cur = s;
    while (!vis[cur]) {
      vis[cur] = 1;
      cyc.push_back(cur);
      cur = nxt[cur];
    }
    int rot = 0;
    for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
      if (cyc[i] < m) { rot = i; break; }
    }
    std::rotate(cyc.begin(), cyc.begin() + rot, cyc.end());

    std::vector<uint8_t> f(cyc.size());
    for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
      f[i] = flip[cyc[i]];
    }
    constraints.push_back(cycleConstraint(f));
  }

  std::sort(constraints.begin(), constraints.end(),
            [](const auto &x, const auto &y) { return x.first < y.first; });

  long long M = 1;
  std::vector<long long> sols = {0};
  for (auto &cc : constraints) {
    mergeCRT(sols, M, cc.second, cc.first);
  }

  long long best = std::numeric_limits<long long>::max();
  for (long long r : sols) {
    if (r == 0) continue;
    if (r < best) best = r;
  }
  if (best == std::numeric_limits<long long>::max()) best = M;
  return best;
}

static int findIntervalIndex(const std::vector<Point> &pts, const Point &x,
                             const FieldContext &ctx) {
  int lo = 0;
  int hi = static_cast<int>(pts.size());
  while (lo + 1 < hi) {
    int mid = (lo + hi) / 2;
    if (cmpPoints(pts[mid], x, ctx.ab, ctx.c) <= 0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

static long long F_nonsquare(int a, int b, int c) {
  FieldContext ctx{a, b, c, static_cast<long long>(a) * b, std::sqrt(static_cast<long double>(c))};

  Point s0 = canonical01(ctx.ab / a, 0, ctx);
  Point s1 = canonical01(ctx.ab / b, 0, ctx);
  Point s2 = canonical01(0, 1, ctx);
  Point s[3] = {s0, s1, s2};
  Point thr[3] = {oneMinus(s0, ctx), oneMinus(s1, ctx), oneMinus(s2, ctx)};

  std::array<std::unordered_set<Point, PointHash>, 3> sets;
  std::array<std::vector<Point>, 3> stack;
  for (int ph = 0; ph < 3; ++ph) {
    sets[ph].reserve(1024);
    sets[ph].insert(Point{0, 0});
    sets[ph].insert(s[ph]);
    stack[ph].push_back(Point{0, 0});
    stack[ph].push_back(s[ph]);
  }

  while (true) {
    bool progressed = false;
    for (int ph = 0; ph < 3; ++ph) {
      while (!stack[ph].empty()) {
        progressed = true;
        Point y = stack[ph].back();
        stack[ph].pop_back();
        int prev = (ph + 2) % 3;
        int cmp = cmpPoints(y, thr[prev], ctx.ab, ctx.c);
        if (cmp < 0) {
          Point x = addP(y, s[prev], ctx);
          if (sets[prev].insert(x).second) stack[prev].push_back(x);
        }
        if (cmp > 0) {
          Point x = oneMinus(y, ctx);
          if (sets[prev].insert(x).second) stack[prev].push_back(x);
        }
      }
    }
    if (!progressed) break;
  }

  std::vector<Point> pts[3];
  int m = -1;
  for (int ph = 0; ph < 3; ++ph) {
    pts[ph].assign(sets[ph].begin(), sets[ph].end());
    std::sort(pts[ph].begin(), pts[ph].end(),
              [&](const Point &u, const Point &v) {
                return cmpPoints(u, v, ctx.ab, ctx.c) < 0;
              });
    if (m < 0) m = static_cast<int>(pts[ph].size());
    if (m != static_cast<int>(pts[ph].size())) {
      std::cerr << "Partition size mismatch for (" << a << "," << b << "," << c << ")\n";
      return -1;
    }
  }

  int N = 3 * m;
  std::vector<int> nxt(N);
  std::vector<uint8_t> flip(N);

  for (int ph = 0; ph < 3; ++ph) {
    for (int i = 0; i < m; ++i) {
      Point l = pts[ph][i];
      Point r = (i + 1 < m) ? pts[ph][i + 1] : Point{ctx.ab, 0};
      bool inside = cmpPoints(l, s[ph], ctx.ab, ctx.c) < 0;
      Point yLeft = inside ? oneMinus(r, ctx) : subP(l, s[ph], ctx);
      int nxtph = (ph + 1) % 3;
      int j = findIntervalIndex(pts[nxtph], yLeft, ctx);
      int cur = ph * m + i;
      nxt[cur] = nxtph * m + j;
      flip[cur] = inside ? 1 : 0;
    }
  }

  std::vector<char> vis(N, 0);
  std::vector<std::pair<long long, std::vector<long long>>> constraints;

  for (int sIdx = 0; sIdx < N; ++sIdx) {
    if (vis[sIdx]) continue;
    std::vector<int> cyc;
    int cur = sIdx;
    while (!vis[cur]) {
      vis[cur] = 1;
      cyc.push_back(cur);
      cur = nxt[cur];
    }

    int rot = 0;
    for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
      if (cyc[i] < m) { rot = i; break; }
    }
    std::rotate(cyc.begin(), cyc.begin() + rot, cyc.end());

    std::vector<uint8_t> f(cyc.size());
    for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
      f[i] = flip[cyc[i]];
    }
    constraints.push_back(cycleConstraint(f));
  }

  std::sort(constraints.begin(), constraints.end(),
            [](const auto &x, const auto &y) { return x.first < y.first; });

  long long M = 1;
  std::vector<long long> sols = {0};
  for (auto &cc : constraints) {
    mergeCRT(sols, M, cc.second, cc.first);
  }

  long long best = std::numeric_limits<long long>::max();
  for (long long r : sols) {
    if (r == 0) continue;
    if (r < best) best = r;
  }
  if (best == std::numeric_limits<long long>::max()) best = M;
  return best;
}

static inline bool isSquare(int c) {
  int r = static_cast<int>(std::llround(std::sqrt(static_cast<long double>(c))));
  return r * r == c;
}

static long long F(int a, int b, int c) {
  if (isSquare(c)) return F_square(a, b, c);
  return F_nonsquare(a, b, c);
}

static long long computeG(int n, int threads) {
  std::atomic<int> next_a(9);
  std::vector<long long> partial(threads);
  std::vector<std::thread> workers;
  workers.reserve(threads);

  for (int t = 0; t < threads; ++t) {
    workers.emplace_back([&, t]() {
      long long local = 0;
      while (true) {
        int a = next_a.fetch_add(1);
        if (a > n - 2) break;
        for (int b = a + 1; b <= n - 1; ++b) {
          for (int c = b + 1; c <= n; ++c) {
            local += F(a, b, c);
          }
        }
      }
      partial[t] = local;
    });
  }

  for (auto &th : workers) th.join();
  long long sum = 0;
  for (const auto &v : partial) sum += v;
  return sum;
}

static bool validate_samples() {
  struct SampleF { int a, b, c; long long expect; };
  const SampleF f_samples[] = {
      {9, 10, 11, 60},
      {10, 14, 16, 506},
      {15, 16, 17, 785232},
  };
  for (const auto &s : f_samples) {
    long long got = F(s.a, s.b, s.c);
    if (got != s.expect) {
      std::cerr << "F(" << s.a << "," << s.b << "," << s.c << ") = " << got
                << ", expected " << s.expect << "\n";
      return false;
    }
  }

  struct SampleG { int n; long long expect; };
  const SampleG g_samples[] = {
      {11, 60},
      {14, 58020},
      {17, 1269260},
  };
  for (const auto &s : g_samples) {
    long long got = computeG(s.n, 1);
    if (got != s.expect) {
      std::cerr << "G(" << s.n << ") = " << got << ", expected " << s.expect << "\n";
      return false;
    }
  }
  return true;
}

extern "C" long long pe566_answer(void) {
  if (!validate_samples()) {
    std::cerr << "Validation failed." << std::endl;
    return -1;
  }
  int threads = static_cast<int>(std::thread::hardware_concurrency());
  if (threads <= 0) threads = 1;
  return computeG(53, threads);
}
