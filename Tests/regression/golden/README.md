# Golden baselines

One JSON file per regression-test case, written by
`run_regression_tests.py --record` and committed here deliberately -- see
`docs/REGRESSION_TESTS.md` for the schema and the comparison rule.

This file exists only so git tracks this directory before the first
`--record` run produces any `*.json` files; it is not itself read by
anything.
