# Gate A verdict

**Status:** FAIL

First failed criterion: `ratio_ci_excludes_one`

| Criterion | Pass | Observed | Required |
|---|---:|---|---|
| controls_pass | yes | `{"human_n": 2, "human_sensitivity": 1.0, "human_sensitivity_bp": 1.0, "microbial_false_positive_bp_rate": 0.0, "microbial_false_positive_rate": 0.0, "microbial_n": 2, "passed": true}` | sensitivity >= 0.95 and microbial false-positive bp rate <= 0.001 |
| complete_groups | yes | `{"netherlands": 2, "tanzania": 2}` | exactly 2 valid samples per cohort |
| complete_sensitivity_matrix | yes | `[]` | p_count results for ['strict_pair', 'identity_0.90', 'identity_0.95', 'identity_0.98'] |
| minimum_tanzania_propagation | yes | `0.5` | >= 0.01 |
| minimum_mean_ratio | yes | `"infinity"` | >= 1.5 |
| ratio_ci_excludes_one | no | `1.0` | > 1 |
| difference_ci_excludes_zero | no | `0.0` | > 0 |
| permutation_significant | yes | `0.5074626865671642` | < 1.0 |
| minimum_positive_tanzania_samples | yes | `1` | >= 1 |
| p_bp_preserves_direction | yes | `0.5` | > 0 |
| sensitivity_matrix_preserves_direction | yes | `{"identity_0.90": {"netherlands": 0.0, "tanzania": 0.5}, "identity_0.95": {"netherlands": 0.0, "tanzania": 0.5}, "identity_0.98": {"netherlands": 0.0, "tanzania": 0.5}, "strict_pair": {"netherlands": 0.0, "tanzania": 0.5}}` | Tanzania mean > Netherlands mean in every required analysis |
| leave_one_out_preserves_direction | no | `0.0` | all leave-one-sample-out differences > 0 |
| no_dominant_sample | yes | `1.0` | <= 1.0 |
