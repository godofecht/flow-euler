#!/usr/bin/env python3
"""Generate a best-known complexity reference database for Project Euler problems.

Reads problem descriptions from the .flow files and categorises each problem
to assign a best-known time complexity, space complexity, and approach label.

Outputs scripts/best-known.json.

Usage:
    python3 scripts/gen-best-known.py
"""
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def extract_description(n: int) -> str:
    label = f"{n:03d}"
    src = ROOT / f"problems/p{label}.flow"
    if not src.exists():
        return ""
    text = src.read_text(errors="replace")
    lines = text.split("\n")
    comments = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            comments.append(stripped[1:].strip())
        elif stripped == "":
            continue
        else:
            break
    # Skip "Project Euler NNN" header, return the description line(s)
    desc_parts = []
    for c in comments[1:]:
        if c.startswith("Project Euler"):
            continue
        if c == "":
            continue
        desc_parts.append(c)
    return " ".join(desc_parts)


# Curated best-known complexities for well-studied problems.
# Format: problem_number -> (best_time, best_space, approach)
CURATED = {
    1: ("O(1)", "O(1)", "Closed-form inclusion-exclusion"),
    2: ("O(log n)", "O(1)", "Iterate Fibonacci, sum evens"),
    3: ("O(sqrt(n))", "O(1)", "Trial division"),
    4: ("O(n^2)", "O(1)", "Brute-force product, check palindrome"),
    5: ("O(n log n)", "O(1)", "LCM via prime factorisation"),
    6: ("O(1)", "O(1)", "Closed-form sum of squares"),
    7: ("O(n log log n)", "O(n)", "Sieve of Eratosthenes"),
    8: ("O(n)", "O(1)", "Sliding window product"),
    9: ("O(n)", "O(1)", "Euclid parametrisation of Pythagorean triples"),
    10: ("O(n log log n)", "O(n)", "Sieve of Eratosthenes"),
    11: ("O(n)", "O(1)", "Linear scan of grid directions"),
    12: ("O(n log n)", "O(n)", "Sieve-based divisor count"),
    13: ("O(n)", "O(1)", "Big-integer addition"),
    14: ("O(n)", "O(n)", "Memoised Collatz lengths"),
    15: ("O(n^2)", "O(n^2)", "Binomial coefficient via DP or combinatorics"),
    16: ("O(n)", "O(n)", "Big-integer multiplication, digit sum"),
    17: ("O(1)", "O(1)", "Direct word count"),
    18: ("O(n^2)", "O(n^2)", "Bottom-up DP on triangle"),
    19: ("O(1)", "O(1)", "Zeller congruence or direct count"),
    20: ("O(n)", "O(n)", "Big-integer factorial, digit sum"),
    21: ("O(n log log n)", "O(n)", "Sieve-based sum of proper divisors"),
    22: ("O(n log n)", "O(n)", "Sort names, compute scores"),
    23: ("O(n log log n)", "O(n)", "Sieve of abundant sums"),
    24: ("O(n)", "O(1)", "Factorial number system"),
    25: ("O(n)", "O(1)", "Iterate Fibonacci, check digit count"),
    26: ("O(n)", "O(n)", "Long division cycle detection"),
    27: ("O(n^2)", "O(n)", "Sieve + quadratic prime search"),
    28: ("O(1)", "O(1)", "Closed-form diagonal sum"),
    29: ("O(n^2 log n)", "O(n^2)", "Distinct powers via set or prime factoring"),
    30: ("O(n * d)", "O(1)", "Brute-force digit power sum"),
    31: ("O(n * m)", "O(n)", "Coin change DP"),
    32: ("O(n^2)", "O(1)", "Pandigital product search"),
    33: ("O(n^2)", "O(1)", "Digit-cancelling fraction search"),
    34: ("O(n * d)", "O(1)", "Brute-force digit factorial sum"),
    35: ("O(n log log n)", "O(n)", "Sieve + circular prime check"),
    36: ("O(n)", "O(1)", "Check palindrome in base 2 and 10"),
    37: ("O(n log log n)", "O(n)", "Sieve + truncatable prime check"),
    38: ("O(n)", "O(1)", "Pandigital concatenation search"),
    39: ("O(n)", "O(1)", "Pythagorean triple perimeter count"),
    40: ("O(1)", "O(1)", "Champernowne digit positions"),
    41: ("O(n!)", "O(n)", "Pandigital prime search"),
    42: ("O(n log n)", "O(n)", "Sort + triangle number check"),
    43: ("O(n!)", "O(n)", "Pandigital substring divisibility"),
    44: ("O(n)", "O(n)", "Pentagonal number difference search"),
    45: ("O(n)", "O(n)", "Triangular-pentagonal-hexagonal search"),
    46: ("O(n log log n)", "O(n)", "Sieve + Goldbach check"),
    47: ("O(n log log n)", "O(n)", "Sieve + distinct prime factor count"),
    48: ("O(n log n)", "O(1)", "Modular exponentiation"),
    49: ("O(n log log n)", "O(n)", "Sieve + prime permutation search"),
    50: ("O(n^2 / log n)", "O(n)", "Consecutive prime sum search"),
    52: ("O(n log n)", "O(1)", "Permutation check on multiples"),
    53: ("O(n^2)", "O(n^2)", "Pascal triangle combinatorics"),
    55: ("O(n * k)", "O(n)", "Lychrel number iteration"),
    56: ("O(n^2)", "O(n)", "Big-integer power, digit sum"),
    58: ("O(n)", "O(1)", "Spiral prime ratio iteration"),
    63: ("O(n)", "O(1)", "Count n-digit nth powers"),
    92: ("O(n)", "O(n)", "Square digit chain memoisation"),
    97: ("O(log n)", "O(1)", "Modular exponentiation"),
    99: ("O(n)", "O(1)", "Compare logarithms of large numbers"),
    100: ("O(log n)", "O(1)", "Pell equation / continued fractions"),
    206: ("O(n)", "O(1)", "Concealed square search"),
    407: ("O(n log n)", "O(n)", "Idempotents via divisor analysis"),
    500: ("O(n log n)", "O(n)", "Smallest number with 2^500500 divisors"),
}


def heuristic_complexity(desc: str) -> tuple:
    """Assign a best-known complexity based on problem description keywords."""
    d = desc.lower()

    # Closed-form / constant
    if any(k in d for k in ["closed-form", "arithmetic series", "inclusion-exclusion",
                            "zeller", "direct count", "direct word count"]):
        return ("O(1)", "O(1)", "Closed-form formula")

    # Big integer / digit sum
    if any(k in d for k in ["sum of digits", "digit sum", "first ten digits",
                            "digital sum", "last five non-zero", "last 5 digits",
                            "last 8 digits", "last 9 digits", "last digits",
                            "last non-zero"]):
        return ("O(n)", "O(n)", "Big-integer arithmetic")

    # Modular exponentiation / tetration
    if any(k in d for k in ["modular", "self powers", "tetration", "↑↑",
                            "mod p_n^2", "mod a^2"]):
        return ("O(log n)", "O(1)", "Modular exponentiation")

    # Prime sieve
    if any(k in d for k in ["prime below", "primes below", "sum of all primes",
                            "nth prime", "10001st prime", "prime permutation",
                            "circular prime", "truncatable prime", "consecutive prime",
                            "prime family", "prime concatenation", "prime chain",
                            "totient chain", "prime-proof", "sqube",
                            "semiprime", "prime generating", "prime pair",
                            "prime triple", "prime summation", "prime power",
                            "count primes", "prime digit"]):
        return ("O(n log log n)", "O(n)", "Sieve of Eratosthenes")

    # Prime factorisation
    if any(k in d for k in ["prime factor", "largest prime factor", "distinct prime",
                            "prime factors", "factorisation", "factorization",
                            "rad(", "sorted radical", "abc-hit",
                            "rad(n)", "radical"]):
        return ("O(sqrt(n))", "O(1)", "Trial division or Pollard rho")

    # Totient / phi
    if any(k in d for k in ["totient", "phi(n)", "phi(", "n/phi",
                            "euler totient", "sum of phi", "fractions a/b",
                            "counting fractions", "fraction immediately",
                            "ambiguous rational"]):
        return ("O(n log log n)", "O(n)", "Sieve-based totient computation")

    # Mobius / squarefree
    if any(k in d for k in ["squarefree", "square-free", "mobius",
                            "mobius inversion", "mu(n)"]):
        return ("O(n log log n)", "O(n)", "Mobius sieve")

    # Collatz
    if "collatz" in d:
        return ("O(n)", "O(n)", "Memoised Collatz iteration")

    # Fibonacci
    if "fibonacci" in d:
        return ("O(log n)", "O(1)", "Matrix exponentiation")

    # Lattice paths / grid
    if any(k in d for k in ["lattice path", "grid path", "lattice",
                            "crack-free wall", "wall"]):
        return ("O(n^2)", "O(n^2)", "Combinatorial or DP counting")

    # Triangle path
    if any(k in d for k in ["path sum", "triangle"]):
        return ("O(n^2)", "O(n^2)", "Bottom-up DP")

    # Palindrome
    if "palindrome" in d or "palindromic" in d:
        return ("O(n^2)", "O(1)", "Brute-force or constructive search")

    # Permutation / pandigital
    if any(k in d for k in ["permutation", "lexicographic", "pandigital",
                            "anagramic", "digit substitution", "number mind",
                            "passcode", "keylog"]):
        return ("O(n!)", "O(n)", "Permutation enumeration or constraint search")

    # Divisor / amicable / abundant
    if any(k in d for k in ["divisor", "divisors", "amicable", "abundant",
                            "proper divisor", "sum of proper", "sigma",
                            "sigma2", "sum of divisors"]):
        return ("O(n log log n)", "O(n)", "Sieve-based divisor sums")

    # Coin change / partition DP
    if any(k in d for k in ["coin", "ways to", "count ways", "how many ways",
                            "sum of primes in", "prime in over",
                            "fill_ways", "partition", "pentagonal recurrence",
                            "p(n)", "product-sum"]):
        return ("O(n * m)", "O(n)", "Dynamic programming or generating function")

    # Pythagorean
    if "pythagorean" in d:
        return ("O(n)", "O(1)", "Euclid parametrisation")

    # Factorial
    if "factorial" in d:
        return ("O(n)", "O(n)", "Big-integer factorial")

    # Digit power
    if any(k in d for k in ["digit power", "fifth power", "fourth power"]):
        return ("O(n * d)", "O(1)", "Brute-force digit power sum")

    # Square digit chain
    if "square digit" in d:
        return ("O(n)", "O(n)", "Memoised chain iteration")

    # Continued fraction / Pell
    if any(k in d for k in ["continued fraction", "pell", "x^2 - d*y^2",
                            "x² − d y²", "minimal solution"]):
        return ("O(sqrt(n))", "O(1)", "Continued fraction convergents")

    # Matrix exponentiation / linear recurrence
    if any(k in d for k in ["matrix exponentiation", "matrix power",
                            "tribonacci", "linear recurrence",
                            "tours on", "t(10^12)"]):
        return ("O(log n)", "O(n^2)", "Matrix exponentiation")

    # Markov chain / expected value / probability
    if any(k in d for k in ["markov", "expected", "probability",
                            "steady state", "flea", "the chase",
                            "the race", "disc game", "expected singles",
                            "empty squares", "paper envelope"]):
        return ("O(n * s^2)", "O(s^2)", "Markov chain or DP over states")

    # Graph / MST / network
    if any(k in d for k in ["mst", "kruskal", "minimum spanning",
                            "network", "graph", "maximum saving"]):
        return ("O(e log v)", "O(v)", "Kruskal or Prim algorithm")

    # Sudoku / backtracking
    if any(k in d for k in ["sudoku", "backtrack", "constraint"]):
        return ("O(n!)", "O(n^2)", "Backtracking search")

    # Digit DP
    if any(k in d for k in ["digit dp", "balanced number", "18-digit",
                            "20-digit", "hex number", "no digit occurs",
                            "consecutive digits", "bouncy",
                            "non-bouncy", "increasing + decreasing",
                            "reversible number"]):
        return ("O(n * d * s)", "O(d * s)", "Digit DP")

    # Subset sum / subset
    if any(k in d for k in ["subset sum", "subset", "special sum set",
                            "uniquely-occurring", "subset pairs"]):
        return ("O(n * sum)", "O(sum)", "Subset sum DP")

    # Binary / hyperbinary / fusc
    if any(k in d for k in ["hyperbinary", "fusc", "stern",
                            "binary expansion", "shortened binary"]):
        return ("O(log n)", "O(log n)", "Stern-Brocot or binary representation")

    # Hamming number
    if "hamming" in d:
        return ("O(n)", "O(n)", "Priority queue generation")

    # Laser / reflection simulation
    if any(k in d for k in ["laser", "reflection", "laserbeam",
                            "elliptical mirror"]):
        return ("O(n)", "O(1)", "Iterative reflection simulation")

    # Pascal's triangle / binomial
    if any(k in d for k in ["pascal", "binomial", "pascal's pyramid",
                            "pascal row"]):
        return ("O(n^2)", "O(n)", "Pascal triangle computation")

    # Capacitance / network
    if any(k in d for k in ["capacitance", "series/parallel"]):
        return ("O(n^2)", "O(n^2)", "Network enumeration")

    # Tiling
    if any(k in d for k in ["tiling", "triomino"]):
        return ("O(n * s)", "O(s)", "Tiling DP with state")

    # Dragon curve
    if any(k in d for k in ["dragon", "heighway"]):
        return ("O(log n)", "O(log n)", "Recursive curve computation")

    # Numerical / curve area
    if any(k in d for k in ["blancmange", "curve", "area fraction",
                            "descartes", "circle iteration",
                            "uncovered area"]):
        return ("O(n)", "O(1)", "Numerical iteration")

    # Minkowski sum
    if "minkowski" in d:
        return ("O(n log n)", "O(n)", "Convex polygon sum")

    # Sum of squares / representable
    if any(k in d for k in ["sum of squares", "representable as",
                            "a^2+b^2", "a^2+2b^2", "a^2+3b^2",
                            "a^2+7b^2"]):
        return ("O(n log log n)", "O(n)", "Sieve-based representation count")

    # Alexandrian integer
    if "alexandrian" in d:
        return ("O(n log n)", "O(n)", "Divisor-based enumeration")

    # Sphere packing
    if any(k in d for k in ["sphere packing", "pipe of radius"]):
        return ("O(n^2)", "O(n)", "Greedy packing simulation")

    # Repunit
    if "repunit" in d:
        return ("O(sqrt(n))", "O(1)", "Repunit divisibility")

    # Roman numerals
    if "roman" in d:
        return ("O(n)", "O(1)", "Direct numeral conversion")

    # XOR / encryption
    if any(k in d for k in ["xor", "encryption", "decrypt"]):
        return ("O(n * k)", "O(n)", "Key search and XOR decryption")

    # Poker / card game
    if any(k in d for k in ["poker", "hand", "player 1", "card",
                            "magic 5-gon", "cube digit"]):
        return ("O(n)", "O(1)", "Direct hand evaluation or enumeration")

    # Monopoly / board game
    if "monopoly" in d:
        return ("O(s^3)", "O(s^2)", "Markov chain steady state")

    # Rectangles in grid
    if any(k in d for k in ["rectangles in", "cross-hatched"]):
        return ("O(n^2)", "O(1)", "Closed-form or combinatorial counting")

    # Cuboid / shortest path
    if any(k in d for k in ["cuboid", "shortest cuboid"]):
        return ("O(n)", "O(n)", "Enumerate and count")

    # Figurate number / cyclic
    if any(k in d for k in ["figurate", "cyclic 4-digit"]):
        return ("O(n^2)", "O(n)", "Figurate number enumeration")

    # Generating function / FIT
    if any(k in d for k in ["generating function", "fit", "fit for"]):
        return ("O(n^2)", "O(n)", "Polynomial interpolation")

    # Multiplication chain
    if any(k in d for k in ["multiplication chain", "multiplication-chain",
                            "m(k)"]):
        return ("O(n log n)", "O(n)", "Addition chain DP")

    # Hexagonal tile
    if any(k in d for k in ["hexagonal tile", "pd(n)"]):
        return ("O(n)", "O(1)", "Tile difference enumeration")

    # Progressive square
    if "progressive square" in d:
        return ("O(n log n)", "O(n)", "Divisor-based enumeration")

    # Diophantine / 1/a + 1/b
    if any(k in d for k in ["1/a + 1/b", "1/x + 1/y", "diophantine"]):
        return ("O(n log n)", "O(1)", "Divisor-based Diophantine solver")

    # Line segment intersection
    if any(k in d for k in ["line segment", "intersection", "bbs line"]):
        return ("O(n^2 log n)", "O(n)", "Sweep line or pairwise check")

    # Grid of digits
    if any(k in d for k in ["4x4 grid", "grid of digits",
                            "row/column/diagonal"]):
        return ("O(n^d)", "O(n^d)", "Brute-force grid search")

    # Laminae / hollow square
    if any(k in d for k in ["laminae", "hollow square", "tile total"]):
        return ("O(n)", "O(1)", "Enumerative counting")

    # Golden triple
    if "golden triple" in d:
        return ("O(n^2)", "O(n)", "Triple enumeration")

    # Non-terminating decimal
    if any(k in d for k in ["non-terminating", "terminating"]):
        return ("O(n log n)", "O(1)", "Prime factor check")

    # Prize string / DP over states
    if any(k in d for k in ["prize string", "no 3 consecutive",
                            "at most one"]):
        return ("O(n * s)", "O(s)", "DP over states")

    # Best approximation / continued fraction
    if any(k in d for k in ["best approximation", "denominator of best",
                            "approximation to sqrt"]):
        return ("O(log n)", "O(1)", "Continued fraction convergents")

    # Iterative sequence / floor
    if any(k in d for k in ["floor(2^", "iterative sequence",
                            "u_{n+1}", "u_n+u_{n+1}"]):
        return ("O(n)", "O(1)", "Iterative sequence simulation")

    # Dice
    if any(k in d for k in ["dice", "d4", "d6", "d12", "9d4", "6d6",
                            "20d12"]):
        return ("O(n * s)", "O(s)", "Dice distribution DP")

    # Robot / closed path
    if any(k in d for k in ["robot", "closed path", "five-degree"]):
        return ("O(n * s)", "O(s)", "Path counting DP")

    # Truth assignment / circular
    if any(k in d for k in ["truth assignment", "circular logic"]):
        return ("O(n)", "O(n)", "Cycle decomposition")

    # Binary code / Huffman
    if any(k in d for k in ["binary code", "huffman", "minimal cost"]):
        return ("O(n log n)", "O(n)", "Huffman coding")

    # Semidivisible
    if "semidivisible" in d:
        return ("O(n log log n)", "O(n)", "Sieve-based computation")

    # Arithmetic-geometric series
    if any(k in d for k in ["arithmetic-geometric", "s(r)"]):
        return ("O(log n)", "O(1)", "Binary search on series")

    # Hysteresis
    if "hysteresis" in d:
        return ("O(n log n)", "O(n)", "Scan with comparison")

    # BBS / digit stream
    if any(k in d for k in ["bbs", "digit stream", "infinite string"]):
        return ("O(n)", "O(n)", "Period finding and digit sum")

    # Sliding block / puzzle
    if any(k in d for k in ["sliding block", "puzzle", "shortest path"]):
        return ("O(n!)", "O(n!)", "BFS or A* search")

    # Hyperbola / packing
    if any(k in d for k in ["hyperbola", "packing", "square index"]):
        return ("O(n log n)", "O(n)", "Geometric enumeration")

    # Rounded square root
    if "rounded square root" in d:
        return ("O(n log n)", "O(1)", "Newton iteration")

    # Tatami
    if "tatami" in d:
        return ("O(n * s)", "O(s)", "Tiling DP")

    # Reachable integers
    if "reachable" in d:
        return ("O(n * d)", "O(n)", "Expression enumeration")

    # Stone game / game theory
    if any(k in d for k in ["stone game", "losing", "game config"]):
        return ("O(n^3)", "O(n^2)", "Game theory DP")

    # Pivot / square sum
    if any(k in d for k in ["pivot", "pivoted square"]):
        return ("O(n log n)", "O(n)", "Search with pruning")

    # Engineers' paradise
    if "engineers' paradise" in d or "engineers paradise" in d:
        return ("O(n)", "O(1)", "Search with primality checks")

    # Binary circle
    if "binary circle" in d:
        return ("O(2^n)", "O(2^n)", "Bitmask enumeration")

    # Pseudo-Fortunate
    if "pseudo-fortunate" in d:
        return ("O(n log log n)", "O(n)", "Sieve-based search")

    # Square cutting
    if any(k in d for k in ["square cutting", "c(30)"]):
        return ("O(n^2)", "O(n^2)", "Combinatorial counting")

    # CRT / modular
    if any(k in d for k in ["x^3 ≡ 1", "crt", "chinese remainder",
                            "mod 13082761331670030"]):
        return ("O(n log n)", "O(n)", "Chinese Remainder Theorem")

    # Divisibility multiplier
    if "divisibility multiplier" in d:
        return ("O(n log log n)", "O(n)", "Sieve-based computation")

    # f(pq,pr,qr) / triangle
    if any(k in d for k in ["f(pq,pr,qr)", "f(pq", "pq,pr"]):
        return ("O(n^2)", "O(n)", "Triple enumeration")

    # f(m,n) / pizza cutting
    if any(k in d for k in ["f(m,n)", "pizza"]):
        return ("O(n^2)", "O(n)", "Double enumeration")

    # Ackermann
    if any(k in d for k in ["a(n,n)", "ackermann"]):
        return ("O(n * a)", "O(n)", "Modular Ackermann computation")

    # Steady square / base-14
    if any(k in d for k in ["steady square", "base-14", "base 14"]):
        return ("O(n * d)", "O(d)", "Digit DP in base b")

    # Probability / P(exactly
    if any(k in d for k in ["p(exactly", "p(exactly 20",
                            "exactly 20 hits"]):
        return ("O(n * s)", "O(s)", "Binomial or DP")

    # Quadtree / encoding
    if any(k in d for k in ["quadtree", "encoding length"]):
        return ("O(n)", "O(log n)", "Recursive tree traversal")

    # Modular computation
    if any(k in d for k in ["mod 14^8", "mod 10^8", "mod 10^9",
                            "mod 61^10", "mod 3^15", "mod 10^10",
                            "s(11^12)", "n(61,"]):
        return ("O(n log n)", "O(n)", "Modular DP or matrix exponentiation")

    # Eulerian / non-crossing
    if any(k in d for k in ["eulerian", "non-crossing"]):
        return ("O(n * s)", "O(s)", "DP over cycle structure")

    # Panaitopol prime
    if "panaitopol" in d:
        return ("O(n log log n)", "O(n)", "Sieve + parametric search")

    # Radius multiplicity
    if "radius multiplicity" in d:
        return ("O(n^2)", "O(n)", "Geometric enumeration")

    # Zeckendorf
    if "zeckendorf" in d:
        return ("O(n log n)", "O(log n)", "Zeckendorf representation sum")

    # Parallelogram / incenter
    if any(k in d for k in ["parallelogram", "incenter"]):
        return ("O(n^2)", "O(n)", "Geometric enumeration")

    # Protein / H-H contact
    if any(k in d for k in ["protein", "h-h contact", "max h-h"]):
        return ("O(n * s)", "O(s)", "Self-avoiding walk simulation")

    # Optimum / optimisation
    if any(k in d for k in ["optimum", "optimal", "minimise", "minimize",
                            "maximise", "maximize", "minimum", "maximum",
                            "least n", "smallest n", "largest"]):
        return ("O(n log n)", "O(n)", "Search with pruning or sieve")

    # Counting / how many
    if any(k in d for k in ["how many", "count n", "count of",
                            "number of", "sum of all", "sum of n",
                            "sum of", "find the", "find x"]):
        return ("O(n log log n)", "O(n)", "Sieve or enumeration")

    # Default: unknown
    return ("?", "?", "Not curated")


def main():
    results = {}

    # Get all problem numbers from answers.txt
    with open(ROOT / "answers.txt") as f:
        for line in f:
            line = line.strip()
            if line and line[0].isdigit():
                parts = line.split(None, 1)
                if len(parts) == 2:
                    n = int(parts[0])
                    desc = extract_description(n)
                    if n in CURATED:
                        bt, bs, approach = CURATED[n]
                    else:
                        bt, bs, approach = heuristic_complexity(desc)
                    results[n] = {
                        "best_time": bt,
                        "best_space": bs,
                        "approach": approach,
                    }

    out = ROOT / "scripts" / "best-known.json"
    out.write_text(json.dumps(results, indent=2, sort_keys=True))
    print(f"Wrote {len(results)} entries to {out}")

    # Summary
    curated = sum(1 for n in results if n in CURATED)
    unknown = sum(1 for n in results if results[n]["best_time"] == "?")
    print(f"  Curated: {curated}")
    print(f"  Heuristic: {len(results) - curated - unknown}")
    print(f"  Unknown: {unknown}")


if __name__ == "__main__":
    main()
