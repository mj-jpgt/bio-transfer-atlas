# Gate A preregistration

The Gate A thresholds, sample list, random seed, replacement policy, controls,
and sensitivity analyses are frozen before cohort results are inspected.
Machine-readable values in `config/thresholds.yaml` and
`config/stage0_selection.yaml` are authoritative.

## Question

After conventional GRCh38 host filtration, what fraction of confidently
human-derived assembled sequence propagates into quality-passing,
apparently novel bacterial or archaeal bins, and is that fraction higher in the
Tanzania cohort than in the Netherlands cohort?

## Primary endpoints

For each sample:

- `p_count`: propagated human-derived contigs divided by all confidently
  human-derived assembled contigs.
- `p_bp`: propagated human-derived base pairs divided by all confidently
  human-derived assembled base pairs.
- Presence of any endpoint bin containing human-derived sequence.
- Human-derived fraction of each endpoint bin at contact, material, and
  dominant tiers.

A sample with no confidently human-derived assembled contigs receives
`p_count = 0`, `p_bp = 0`, and an explicit zero-denominator flag.

## Gate decision

Gate A is `PASS` only when every rule in `thresholds.yaml:gate_a` is satisfied.
Any failed scientific rule is a scientific `FAIL`; broken controls, downloads,
or software are operational failures that must be repaired and rerun without
changing scientific thresholds.

The final report must state the exact first link in the attrition chain at which
human-derived sequence disappears or, on a pass, the exact quantity reaching
the endpoint.

