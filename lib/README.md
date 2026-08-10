# Shared Flow libraries

## `euler.nt` (`lib/nt.flow`)

Common i64 number-theory helpers: `gcd`, `lcm`, `isqrt`, `mulmod`, `mod_pow`, `is_prime`.

### Import

From a problem file under `problems/`:

```flow
import euler.nt { gcd, lcm, is_prime }
```

Resolution comes from the repo-root `flow.toml`:

```toml
[paths]
euler = "lib"   # import euler.nt → lib/nt.flow
```

`flow run` / `scripts/run.sh` walk upward from the source file to find `flow.toml`, so no cwd tricks are required when launching from the repo root:

```bash
FLOW_BIN=/path/to/flow FLOW_HOST=python ./scripts/run.sh 5
# or
FLOW_HOST=python /path/to/flow run problems/p005.flow
```

No `run.sh` changes were needed for import resolution.
