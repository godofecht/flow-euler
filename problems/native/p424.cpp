// Project Euler 424 — Cryptic Kakuro CSP with GAC propagation.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static constexpr int ALL10 = (1 << 10) - 1;
static constexpr int DIGITS_19 = (1 << 10) - 2; // bits 1..9

static int popcount(int x) { return __builtin_popcount((unsigned)x); }

static std::vector<int> iter_bits(int mask) {
  std::vector<int> out;
  while (mask) {
    int b = mask & -mask;
    out.push_back(__builtin_ctz((unsigned)b));
    mask ^= b;
  }
  return out;
}

static std::vector<std::string> split_tokens(const std::string &line) {
  std::vector<std::string> out;
  std::string cur;
  int depth = 0;
  for (char ch : line) {
    if (ch == ',' && depth == 0) {
      out.push_back(cur);
      cur.clear();
      continue;
    }
    if (ch == '(') ++depth;
    else if (ch == ')') --depth;
    cur.push_back(ch);
  }
  out.push_back(cur);
  return out;
}

enum CellType { CT_X, CT_O, CT_L, CT_C };

struct Cell {
  CellType typ;
  char letter; // for L
  std::string h, v; // for C
};

static std::vector<Cell> parse_puzzle(const std::string &line, int &n) {
  auto toks = split_tokens(line);
  n = std::stoi(toks[0]);
  std::vector<Cell> grid;
  grid.reserve((size_t)n * n);
  for (size_t ti = 1; ti < toks.size(); ++ti) {
    const std::string &t = toks[ti];
    Cell c{};
    if (t == "X") c.typ = CT_X;
    else if (t == "O") c.typ = CT_O;
    else if (t.size() == 1 && t[0] >= 'A' && t[0] <= 'J') {
      c.typ = CT_L;
      c.letter = t[0];
    } else if (!t.empty() && t.front() == '(' && t.back() == ')') {
      c.typ = CT_C;
      std::string inside = t.substr(1, t.size() - 2);
      size_t pos = 0;
      while (pos < inside.size()) {
        size_t comma = inside.find(',', pos);
        std::string part = inside.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (!part.empty()) {
          if (part[0] == 'h') c.h = part.substr(1);
          else if (part[0] == 'v') c.v = part.substr(1);
        }
        if (comma == std::string::npos) break;
        pos = comma + 1;
      }
    }
    grid.push_back(c);
  }
  return grid;
}

// Precomputed perms by (len,sum)
static std::vector<std::vector<int>> PERMS[7][46]; // L=1..6, sum up to 45

static void init_perms() {
  static bool done = false;
  if (done) return;
  done = true;
  int digits[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (int L = 1; L <= 6; ++L) {
    std::vector<int> idx(9);
    std::iota(idx.begin(), idx.end(), 0);
    // generate combinations then permutations via next_permutation on chosen set
    std::string mask(9, 0);
    std::fill(mask.end() - L, mask.end(), 1);
    do {
      std::vector<int> chosen;
      for (int i = 0; i < 9; ++i)
        if (mask[i]) chosen.push_back(digits[i]);
      std::sort(chosen.begin(), chosen.end());
      do {
        int s = 0;
        for (int d : chosen) s += d;
        PERMS[L][s].push_back(chosen);
      } while (std::next_permutation(chosen.begin(), chosen.end()));
    } while (std::next_permutation(mask.begin(), mask.end()));
  }
}

struct RunConstraint {
  std::vector<int> cell_vars;
  std::vector<int> sum_letters;
  std::vector<int> involved;
  std::unordered_map<int, int> inv_index;
  int length = 0;
  std::vector<int> sums_list;
  bool dup_cells = false;
  bool sum_same = false;
};

struct CSP {
  std::vector<int> domains;
  std::vector<RunConstraint> runs;
  std::vector<std::vector<int>> var_to_runs;
};

static CSP build_csp(const std::string &line) {
  init_perms();
  int n;
  auto grid = parse_puzzle(line, n);
  CSP csp;
  csp.domains.assign(10, ALL10);
  for (auto &cell : grid)
    if (cell.typ == CT_L) csp.domains[cell.letter - 'A'] &= DIGITS_19;

  std::unordered_map<int, int> cell_to_var; // r*n+c -> var
  int ocount = 0;
  for (int r = 0; r < n; ++r)
    for (int c = 0; c < n; ++c) {
      auto &cell = grid[r * n + c];
      if (cell.typ == CT_O) {
        cell_to_var[r * n + c] = 10 + ocount;
        ++ocount;
      } else if (cell.typ == CT_L) {
        cell_to_var[r * n + c] = cell.letter - 'A';
      }
    }
  for (int i = 0; i < ocount; ++i) csp.domains.push_back(DIGITS_19);

  auto add_run = [&](const std::string &code, const std::vector<int> &cell_vars) {
    if (cell_vars.empty()) return;
    RunConstraint run;
    run.cell_vars = cell_vars;
    if (code.size() == 1) run.sum_letters = {code[0] - 'A'};
    else run.sum_letters = {code[0] - 'A', code[1] - 'A'};
    run.length = (int)cell_vars.size();
    for (int s = 1; s <= 45; ++s)
      if (!PERMS[run.length][s].empty()) run.sums_list.push_back(s);
    std::vector<int> inv = cell_vars;
    inv.insert(inv.end(), run.sum_letters.begin(), run.sum_letters.end());
    std::sort(inv.begin(), inv.end());
    inv.erase(std::unique(inv.begin(), inv.end()), inv.end());
    run.involved = inv;
    for (size_t i = 0; i < inv.size(); ++i) run.inv_index[inv[i]] = (int)i;
    {
      std::vector<int> tmp = cell_vars;
      std::sort(tmp.begin(), tmp.end());
      run.dup_cells = std::unique(tmp.begin(), tmp.end()) != tmp.end();
    }
    run.sum_same = run.sum_letters.size() == 2 && run.sum_letters[0] == run.sum_letters[1];
    csp.runs.push_back(std::move(run));
  };

  for (int r = 0; r < n; ++r)
    for (int c = 0; c < n; ++c) {
      auto &cell = grid[r * n + c];
      if (cell.typ != CT_C) continue;
      if (!cell.h.empty()) {
        std::vector<int> cvs;
        for (int cc = c + 1; cc < n; ++cc) {
          auto &g = grid[r * n + cc];
          if (g.typ != CT_O && g.typ != CT_L) break;
          cvs.push_back(cell_to_var[r * n + cc]);
        }
        add_run(cell.h, cvs);
      }
      if (!cell.v.empty()) {
        std::vector<int> cvs;
        for (int rr = r + 1; rr < n; ++rr) {
          auto &g = grid[rr * n + c];
          if (g.typ != CT_O && g.typ != CT_L) break;
          cvs.push_back(cell_to_var[rr * n + c]);
        }
        add_run(cell.v, cvs);
      }
    }

  csp.var_to_runs.assign(csp.domains.size(), {});
  for (size_t i = 0; i < csp.runs.size(); ++i)
    for (int v : csp.runs[i].involved) csp.var_to_runs[v].push_back((int)i);
  return csp;
}

static std::pair<std::vector<int>, bool> enforce_all_diff(std::vector<int> &dom) {
  int fixed_mask = 0;
  bool seen[10] = {};
  for (int i = 0; i < 10; ++i) {
    if (popcount(dom[i]) == 1) {
      int d = __builtin_ctz((unsigned)dom[i]);
      if (seen[d]) return {{}, true};
      seen[d] = true;
      fixed_mask |= 1 << d;
    }
  }
  std::vector<int> changed;
  for (int i = 0; i < 10; ++i) {
    if (popcount(dom[i]) > 1) {
      int neu = dom[i] & ~fixed_mask;
      if (neu == 0) return {{}, true};
      if (neu != dom[i]) {
        dom[i] = neu;
        changed.push_back(i);
      }
    }
  }
  return {changed, false};
}

static std::pair<std::vector<int>, bool> process_run(const RunConstraint &run, std::vector<int> &dom) {
  if (run.dup_cells) return {{}, true};
  std::vector<int> support(run.involved.size(), 0);
  const auto &cv = run.cell_vars;
  int L = run.length;

  auto try_perm = [&](int S, const std::vector<int> &p, auto mark_sum) {
    for (int i = 0; i < L; ++i)
      if (!(dom[cv[i]] & (1 << p[i]))) return;
    for (int i = 0; i < L; ++i) support[run.inv_index.at(cv[i])] |= 1 << p[i];
    mark_sum();
  };

  if (run.sum_letters.size() == 1) {
    int l = run.sum_letters[0];
    int li = run.inv_index.at(l);
    int doml = dom[l];
    for (int S : run.sums_list) {
      if (S > 9) continue;
      if (!(doml & (1 << S))) continue;
      for (auto &p : PERMS[L][S]) {
        try_perm(S, p, [&] { support[li] |= 1 << S; });
      }
    }
  } else {
    int a = run.sum_letters[0], b = run.sum_letters[1];
    int ai = run.inv_index.at(a), bi = run.inv_index.at(b);
    if (run.sum_same) {
      int doma = dom[a];
      for (int S : run.sums_list) {
        if (S < 10) continue;
        int t = S / 10, o = S % 10;
        if (t != o) continue;
        if (!(doma & (1 << t))) continue;
        for (auto &p : PERMS[L][S]) {
          try_perm(S, p, [&] { support[ai] |= 1 << t; });
        }
      }
    } else {
      int doma = dom[a], domb = dom[b];
      for (int S : run.sums_list) {
        if (S < 10) continue;
        int t = S / 10, o = S % 10;
        if (t == o) continue;
        if (!(doma & (1 << t)) || !(domb & (1 << o))) continue;
        for (auto &p : PERMS[L][S]) {
          try_perm(S, p, [&] {
            support[ai] |= 1 << t;
            support[bi] |= 1 << o;
          });
        }
      }
    }
  }

  std::vector<int> changed;
  for (size_t i = 0; i < run.involved.size(); ++i) {
    int v = run.involved[i];
    int sup = support[i];
    if (sup == 0) return {{}, true};
    int neu = dom[v] & sup;
    if (neu == 0) return {{}, true};
    if (neu != dom[v]) {
      dom[v] = neu;
      changed.push_back(v);
    }
  }
  return {changed, false};
}

static bool propagate(std::vector<int> &dom, const std::vector<RunConstraint> &runs,
                      const std::vector<std::vector<int>> &var_to_runs) {
  std::deque<int> q;
  std::vector<char> inq(runs.size(), 1);
  for (size_t i = 0; i < runs.size(); ++i) q.push_back((int)i);
  while (true) {
    auto [chL, fail] = enforce_all_diff(dom);
    if (fail) return false;
    for (int l : chL)
      for (int ri : var_to_runs[l])
        if (!inq[ri]) {
          q.push_back(ri);
          inq[ri] = 1;
        }
    if (q.empty()) return true;
    int ri = q.front();
    q.pop_front();
    inq[ri] = 0;
    auto [chV, fail2] = process_run(runs[ri], dom);
    if (fail2) return false;
    for (int v : chV)
      for (int rj : var_to_runs[v])
        if (!inq[rj]) {
          q.push_back(rj);
          inq[rj] = 1;
        }
  }
}

static int select_var(const std::vector<int> &dom) {
  int best = -1, best_c = 999;
  for (size_t i = 0; i < dom.size(); ++i) {
    int c = popcount(dom[i]);
    if (c > 1 && c < best_c) {
      best_c = c;
      best = (int)i;
      if (c == 2) break;
    }
  }
  return best;
}

static bool solve_rec(std::vector<int> dom, const std::vector<RunConstraint> &runs,
                      const std::vector<std::vector<int>> &var_to_runs, std::vector<int> &sol) {
  if (!propagate(dom, runs, var_to_runs)) return false;
  int v = select_var(dom);
  if (v < 0) {
    sol = std::move(dom);
    return true;
  }
  for (int d : iter_bits(dom[v])) {
    auto dom2 = dom;
    dom2[v] = 1 << d;
    if (solve_rec(std::move(dom2), runs, var_to_runs, sol)) return true;
  }
  return false;
}

static long long mapping_to_number(const std::vector<int> &sol) {
  long long num = 0;
  for (int i = 0; i < 10; ++i) {
    int d = __builtin_ctz((unsigned)sol[i]);
    num = num * 10 + d;
  }
  return num;
}

static long long solve_puzzle(const std::string &line) {
  CSP csp = build_csp(line);
  std::vector<int> sol;
  if (!solve_rec(csp.domains, csp.runs, csp.var_to_runs, sol)) {
    std::fprintf(stderr, "no solution\n");
    std::abort();
  }
  return mapping_to_number(sol);
}

extern "C" long long pe424_answer(void) {
  const char *paths[] = {"data/p424.txt", "problems/../data/p424.txt", nullptr};
  std::ifstream in;
  for (int i = 0; paths[i]; ++i) {
    in.open(paths[i]);
    if (in) break;
  }
  if (!in) {
    std::fprintf(stderr, "cannot open data/p424.txt\n");
    std::abort();
  }
  long long total = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    // trim CR
    if (!line.empty() && line.back() == '\r') line.pop_back();
    total += solve_puzzle(line);
  }
  return total;
}
