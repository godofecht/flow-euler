# Flow idioms from the Euler corpus

This repository is both a Project Euler corpus and a pressure test for Flow syntax. The conventions here are derived from repeated solution patterns, not from making code look more functional.

## Loop shape follows the bound

Use `for ... in range` when the iteration space is known before the loop starts. Keep `while` when progress or termination depends on values produced by the loop.

```flow
for i in 0..n {
    total = total + values[i]
}
```

Carry propagation, Euclid's algorithm, continued fractions, parser cursors, and searches whose bound changes during iteration remain natural `while` loops.

The goal is to expose iteration intent without hiding mutation or changing the generated class of loop.

## Put representation-specific operations beside the representation

When multiple problems carry the same state tuple and repeat the same operations, promote the tuple to a struct and the operations to a shared module.

Decimal big integers are the current exemplar:

```flow
import euler.digits { Digits, digits_one, digits_mul_u32, digits_sum }
```

`p016`, `p020`, `p055`, and `p056` should share the same LSD-first decimal representation rather than each defining its own `ptr<i32> + len + carry` protocol.

The same rule should be applied next to sieves and other repeated stateful representations when a shared abstraction remains zero-cost.

## Reuse mathematical kernels

Repeated number-theory functions belong in `euler.nt`. A problem should import `gcd`, `lcm`, `isqrt`, `mulmod`, `mod_pow`, or `is_prime` instead of carrying a private copy when the shared function has the required integer semantics.

Do not force reuse when a problem needs a materially different overflow model, integer width, primality algorithm, or performance contract.

## Prefer a live prefix over clearing dead storage

For heap-backed scratch arrays represented by a pointer plus logical length, the logical length is the state boundary. Do not clear the whole capacity between operations unless stale elements can actually be observed.

This removes unnecessary initialization loops in digit-bigint problems and makes the representation invariant explicit.

## Fixed storage should not require giant literals

Repeated zero literals and heap allocations for tiny fixed buffers indicate a missing ergonomic path. Until Flow has a uniformly good fixed-array initialization form, prefer the least noisy representation that is already supported by the compiler and appropriate to the lifetime.

This corpus should continue to feed concrete cases back into Flow language design rather than invent repository-local syntax.

## Resource ownership should become lexical

Hundreds of solutions manually pair `calloc/free` and `fopen/fclose`. That is valid low-level Flow today, but the repetition suggests a language/library idiom for scoped cleanup.

A future Flow facility should make ownership lifetime visible at the binding site while lowering to the same deterministic cleanup. The Euler corpus is a useful acceptance suite for such a feature because early returns and multiple scratch buffers are common.

## Do not replace explicit hot loops for aesthetics

Pipelines, sorting helpers, spans, enums, `match`, and parallel iteration are useful where they expose domain structure. They should not replace compact arithmetic kernels merely to increase feature usage.

An idiomatic rewrite is successful when it improves at least one of these properties without regressing the others: intent, reuse, correctness surface, generated-code quality, or resource safety.

## Adoption rule

A new corpus idiom should have at least two real consumers or eliminate a clearly recurring failure mode. Add the smallest shared abstraction that satisfies those consumers, migrate representative problems, verify their known Euler answers, and only then expand the sweep.

This keeps Flow idioms empirical: repeated code creates a candidate convention; successful reuse turns it into a repository practice; repeated repository practices become evidence for standard-library or language features.
