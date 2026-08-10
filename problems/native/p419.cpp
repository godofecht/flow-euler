// Project Euler 419 — Look and Say via Conway element-decay matrices mod 2^30.
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

static constexpr uint32_t MOD_MASK = (1u << 30) - 1u;

static std::string say(const std::string &s) {
  if (s.empty()) return "";
  std::string out;
  out.reserve(s.size() * 2);
  size_t i = 0, n = s.size();
  while (i < n) {
    char ch = s[i];
    size_t j = i + 1;
    while (j < n && s[j] == ch) ++j;
    out += std::to_string(j - i);
    out.push_back(ch);
    i = j;
  }
  return out;
}

static bool spl00_at(const std::string &s, size_t j) {
  size_t n = s.size();
  if (j + 2 < n && s[j] == '1' && s[j + 1] == '1' && s[j + 2] == '1') return true;
  if (j == n - 1 && s[j] == '1') return false;
  if (j + 1 < n && s[j] == '1' && s[j + 1] == '1') return false;
  if (j + 2 < n && s[j] == '1' && s[j + 1] == '2' && s[j + 2] == '2') return false;
  if (j + 2 < n && s[j] == '1' && s[j + 1] == '3' && s[j + 2] == '3') return false;
  if (j < n && s[j] == '2') return false;
  if (j + 3 < n && s[j] == '3' && s[j + 1] == '1' && s[j + 2] == '1' && s[j + 3] == '1')
    return false;
  if (j + 3 < n && s[j] == '3' && s[j + 1] == '2' && s[j + 2] == '2' && s[j + 3] == '2')
    return false;
  if (j + 1 < n && s[j] == '3' && s[j + 1] == '3') return false;
  return true;
}

static bool spl0_at(const std::string &s, size_t i) {
  size_t n = s.size();
  char ch = s[i];
  if (ch == '1' && i == n - 1) return true;
  if (ch == '1' && i + 2 < n && s[i + 1] == '2' && s[i + 2] == '2')
    return spl00_at(s, i + 3);
  if (ch == '2') return spl00_at(s, i + 1);
  if (ch == '3' && i == n - 1) return true;
  if (ch == '3' && i + 2 < n && s[i + 1] == '2' && s[i + 2] == '2')
    return spl00_at(s, i + 3);
  return false;
}

static std::vector<std::string> split_elements(const std::string &s) {
  std::vector<std::string> parts;
  if (s.empty()) return parts;
  size_t start = 0, n = s.size();
  for (size_t i = 0; i < n; ++i) {
    if (spl0_at(s, i)) {
      parts.push_back(s.substr(start, i + 1 - start));
      start = i + 1;
    }
  }
  if (start < n) parts.push_back(s.substr(start));
  return parts;
}

static std::vector<std::vector<uint32_t>> mat_mul(const std::vector<std::vector<uint32_t>> &A,
                                                 const std::vector<std::vector<uint32_t>> &B) {
  size_t m = A.size();
  std::vector<std::vector<uint32_t>> res(m, std::vector<uint32_t>(m, 0));
  for (size_t i = 0; i < m; ++i) {
    for (size_t k = 0; k < m; ++k) {
      uint32_t a = A[i][k];
      if (!a) continue;
      for (size_t j = 0; j < m; ++j)
        res[i][j] = (res[i][j] + (uint64_t)a * B[k][j]) & MOD_MASK;
    }
  }
  return res;
}

static std::vector<uint32_t> vec_mul(const std::vector<uint32_t> &v,
                                    const std::vector<std::vector<uint32_t>> &M) {
  size_t m = v.size();
  std::vector<uint32_t> out(m, 0);
  for (size_t i = 0; i < m; ++i) {
    uint32_t a = v[i];
    if (!a) continue;
    for (size_t j = 0; j < m; ++j)
      out[j] = (out[j] + (uint64_t)a * M[i][j]) & MOD_MASK;
  }
  return out;
}

static std::vector<uint32_t> vec_mul_pow(std::vector<uint32_t> v,
                                        std::vector<std::vector<uint32_t>> M, int64_t exp) {
  while (exp > 0) {
    if (exp & 1) v = vec_mul(v, M);
    exp >>= 1;
    if (exp) M = mat_mul(M, M);
  }
  return v;
}

static void counts_in_term(const std::string &term, uint32_t &A, uint32_t &B, uint32_t &C) {
  A = B = C = 0;
  for (char ch : term) {
    if (ch == '1') ++A;
    else if (ch == '2') ++B;
    else if (ch == '3') ++C;
  }
}

extern "C" void pe419_answer_print(void) {
  const int64_t n = 1000000000000LL;
  std::string term = "1";
  int steps = (int)std::min<int64_t>(n - 1, 39);
  for (int i = 0; i < steps; ++i) term = say(term);

  if (n <= 40) {
    uint32_t A, B, C;
    counts_in_term(term, A, B, C);
    std::printf("%u,%u,%u\n", A, B, C);
    return;
  }

  auto seed_elements = split_elements(term);
  std::vector<std::string> elems;
  std::unordered_map<std::string, int> idx;

  auto add = [&](const std::string &e) -> int {
    auto it = idx.find(e);
    if (it != idx.end()) return it->second;
    int j = (int)elems.size();
    idx[e] = j;
    elems.push_back(e);
    return j;
  };
  for (auto &e : seed_elements) add(e);
  for (size_t p = 0; p < elems.size(); ++p) {
    auto d = split_elements(say(elems[p]));
    for (auto &child : d) add(child);
  }

  size_t m = elems.size();
  std::vector<std::vector<uint32_t>> M(m, std::vector<uint32_t>(m, 0));
  for (size_t i = 0; i < m; ++i) {
    auto d = split_elements(say(elems[i]));
    for (auto &child : d) M[i][idx[child]] += 1;
  }

  std::vector<uint32_t> ones(m), twos(m), threes(m), v(m, 0);
  for (size_t i = 0; i < m; ++i) {
    for (char ch : elems[i]) {
      if (ch == '1') ++ones[i];
      else if (ch == '2') ++twos[i];
      else if (ch == '3') ++threes[i];
    }
  }
  for (auto &e : seed_elements) v[idx[e]] += 1;

  v = vec_mul_pow(std::move(v), std::move(M), n - 40);
  uint32_t A = 0, B = 0, C = 0;
  for (size_t i = 0; i < m; ++i) {
    if (!v[i]) continue;
    A = (A + (uint64_t)v[i] * ones[i]) & MOD_MASK;
    B = (B + (uint64_t)v[i] * twos[i]) & MOD_MASK;
    C = (C + (uint64_t)v[i] * threes[i]) & MOD_MASK;
  }
  std::printf("%u,%u,%u\n", A, B, C);
}
