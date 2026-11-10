use crate::input::{library_name, OLIGO_ROLE};
use crate::profile::Profile;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

fn fraction(n: f64, d: f64) -> Option<f64> {
    (d > 0.0).then(|| (n / d).clamp(0.0, 1.0))
}

#[derive(Deserialize, Serialize)]
pub struct Pair {
    pub cdna: String,
    pub oligo: String,
    pub cdna_sample: String,
    pub oligo_sample: String,
    /// True when both libraries come from the same sample directory.
    pub expected_pair: bool,
    pub shared_barcodes: usize,
    pub barcode_jaccard: Option<f64>,
    pub cdna_barcode_coverage: Option<f64>,
    pub oligo_barcode_coverage: Option<f64>,
    pub cdna_read_coverage: Option<f64>,
    pub oligo_read_coverage: Option<f64>,
    pub count_cosine: Option<f64>,
    pub log_count_cosine: Option<f64>,
    pub abundance_overlaps: Vec<AbundanceOverlap>,
    pub top_barcode_overlap: Option<f64>,
    pub sufficient_barcode_count: bool,
    pub partial_sampling: bool,
    /// Quantity safeguard only, not a calibrated compatibility threshold.
    pub sufficient_shared_barcode_count: bool,
}

#[derive(Deserialize, Serialize)]
pub struct AbundanceOverlap {
    pub requested_top: usize,
    pub cdna_barcodes: usize,
    pub oligo_barcodes: usize,
    pub shared_barcodes: usize,
    pub overlap: Option<f64>,
}

pub fn compare(cdna: &Profile, oligo: &Profile, min_distinct: usize) -> Pair {
    let mut shared = 0;
    let mut log_dot = 0u128;
    let mut tier_shared = vec![0; cdna.abundance_tiers.len()];
    let tier_filters: Vec<_> = cdna
        .abundance_tiers
        .iter()
        .zip(&oligo.abundance_tiers)
        .enumerate()
        .filter(|(_, (c, o))| {
            c.included_barcodes != cdna.counts.len() || o.included_barcodes != oligo.counts.len()
        })
        .map(|(index, (c, o))| (index, c.minimum_count, o.minimum_count))
        .collect();
    let uniform_log_weights = if cdna.log_weights.len() == 1 && oligo.log_weights.len() == 1 {
        Some((
            *cdna.log_weights.values().next().unwrap(),
            *oligo.log_weights.values().next().unwrap(),
        ))
    } else {
        None
    };
    let (mut cdna_reads, mut oligo_reads, mut dot) = (0u64, 0u64, 0u128);
    let cdna_is_smaller = cdna.counts.len() <= oligo.counts.len();
    let (small, large) = if cdna_is_smaller {
        (&cdna.counts, &oligo.counts)
    } else {
        (&oligo.counts, &cdna.counts)
    };
    for (key, &value) in small {
        if let Some(&other) = large.get(key) {
            let (left, right) = if cdna_is_smaller {
                (value, other)
            } else {
                (other, value)
            };
            shared += 1;
            cdna_reads += left;
            oligo_reads += right;
            dot += left as u128 * right as u128;
            if uniform_log_weights.is_none() {
                log_dot += cdna.log_weights[&left] as u128 * oligo.log_weights[&right] as u128;
            }
            for &(index, minimum_left, minimum_right) in &tier_filters {
                if left >= minimum_left && right >= minimum_right {
                    tier_shared[index] += 1;
                }
            }
        }
    }
    if let Some((left, right)) = uniform_log_weights {
        log_dot = shared as u128 * left as u128 * right as u128;
    }
    for (index, (c, o)) in cdna
        .abundance_tiers
        .iter()
        .zip(&oligo.abundance_tiers)
        .enumerate()
    {
        if c.included_barcodes == cdna.counts.len() && o.included_barcodes == oligo.counts.len() {
            tier_shared[index] = shared;
        }
    }
    Pair {
        cdna: cdna.name.clone(),
        oligo: oligo.name.clone(),
        cdna_sample: cdna.sample.clone(),
        oligo_sample: oligo.sample.clone(),
        expected_pair: cdna.sample == oligo.sample && cdna.role != oligo.role,
        shared_barcodes: shared,
        sufficient_shared_barcode_count: shared >= min_distinct,
        barcode_jaccard: fraction(
            shared as f64,
            (cdna.counts.len() + oligo.counts.len() - shared) as f64,
        ),
        cdna_barcode_coverage: fraction(shared as f64, cdna.counts.len() as f64),
        oligo_barcode_coverage: fraction(shared as f64, oligo.counts.len() as f64),
        cdna_read_coverage: fraction(cdna_reads as f64, cdna.supported_reads as f64),
        oligo_read_coverage: fraction(oligo_reads as f64, oligo.supported_reads as f64),
        count_cosine: fraction(dot as f64, cdna.norm * oligo.norm),
        log_count_cosine: fraction(log_dot as f64, cdna.log_norm * oligo.log_norm),
        abundance_overlaps: cdna
            .abundance_tiers
            .iter()
            .zip(&oligo.abundance_tiers)
            .enumerate()
            .map(|(i, (c, o))| AbundanceOverlap {
                requested_top: c.requested_top,
                cdna_barcodes: c.included_barcodes,
                oligo_barcodes: o.included_barcodes,
                shared_barcodes: tier_shared[i],
                overlap: fraction(
                    tier_shared[i] as f64,
                    c.included_barcodes.min(o.included_barcodes) as f64,
                ),
            })
            .collect(),
        top_barcode_overlap: fraction(
            cdna.top.intersection(&oligo.top).count() as f64,
            cdna.top.len().min(oligo.top.len()) as f64,
        ),
        partial_sampling: cdna
            .files
            .iter()
            .chain(&oligo.files)
            .any(|f| !f.input_fully_scanned),
        sufficient_barcode_count: cdna.counts.len() >= min_distinct
            && oligo.counts.len() >= min_distinct,
    }
}

#[derive(Serialize)]
pub struct Ranking {
    pub cdna: String,
    pub status: &'static str,
    pub eligible_candidates: Vec<String>,
    pub metric_orders: BTreeMap<&'static str, Vec<String>>,
    pub metric_leaders: BTreeMap<&'static str, Vec<String>>,
    pub consensus_leaders: Vec<String>,
    pub preferred_candidate: Option<String>,
    /// Descriptive leader, including weak or sampled evidence; not confirmation.
    pub abundance_candidate: Option<String>,
    pub overlap_candidate: Option<String>,
    /// A tie containing the abundance leader is compatible (true), not unique support.
    pub overlap_agrees: Option<bool>,
    pub coverage_agrees: Option<bool>,
    pub diagnostic_warnings: Vec<&'static str>,
    pub coverage_margin: Option<f64>,
    /// Metrics with a non-tied candidate spread; descriptive for one candidate.
    pub informative_metrics: Vec<&'static str>,
    pub ranking_basis: &'static str,
    pub count_cosine_margin: Option<f64>,
    pub log_count_cosine_margin: Option<f64>,
    /// `<sample>/oligo` implied by the sample directory, when it exists.
    pub expected_candidate: Option<String>,
    pub expected_pair_assessment: &'static str,
}

type Metric = (&'static str, fn(&Pair) -> f64);

pub fn rank(cdna_name: &str, sample: &str, pairs: &[Pair]) -> Ranking {
    let mut eligible: Vec<_> = pairs
        .iter()
        .filter(|p| p.cdna == cdna_name && p.sufficient_barcode_count)
        .collect();
    eligible.sort_by(|a, b| a.oligo.cmp(&b.oligo));
    let mut orders = BTreeMap::new();
    let mut leaders = BTreeMap::new();
    let metrics: [Metric; 4] = [
        ("cdna_read_coverage", |p| {
            p.cdna_read_coverage.unwrap_or(0.0)
        }),
        ("barcode_jaccard", |p| p.barcode_jaccard.unwrap_or(0.0)),
        ("count_cosine", |p| p.count_cosine.unwrap_or(0.0)),
        ("log_count_cosine", |p| p.log_count_cosine.unwrap_or(0.0)),
    ];
    let mut coverage_margin = None;
    let mut count_cosine_margin = None;
    let mut log_count_cosine_margin = None;
    let mut informative_metrics = Vec::new();
    for (name, value) in metrics {
        let mut sorted = eligible.clone();
        sorted.sort_by(|a, b| value(b).total_cmp(&value(a)).then(a.oligo.cmp(&b.oligo)));
        let top = sorted.first().map(|p| value(p));
        let tied: Vec<_> = sorted
            .iter()
            .filter(|p| top.is_some_and(|best| (value(p) - best).abs() <= 1e-12))
            .map(|p| p.oligo.clone())
            .collect();
        if name == "cdna_read_coverage" && sorted.len() >= 2 {
            coverage_margin = Some(value(sorted[0]) - value(sorted[1]));
        }
        if sorted.len() >= 2 {
            let margin = Some(value(sorted[0]) - value(sorted[1]));
            if name == "count_cosine" {
                count_cosine_margin = margin;
            }
            if name == "log_count_cosine" {
                log_count_cosine_margin = margin;
            }
            if tied.len() < sorted.len() {
                informative_metrics.push(name);
            }
        }
        orders.insert(name, sorted.iter().map(|p| p.oligo.clone()).collect());
        leaders.insert(name, tied);
    }
    let abundance_metrics = ["count_cosine", "log_count_cosine"];
    let ranking_basis = if eligible.len() < 2 {
        "descriptive_only"
    } else if abundance_metrics
        .iter()
        .all(|name| leaders[name].len() == eligible.len())
    {
        "no_discriminating_abundance_metrics"
    } else {
        "abundance_only"
    };
    // Overlap and coverage remain diagnostics. They never vote on the winner.
    let consensus: Vec<_> = eligible
        .iter()
        .filter(|p| {
            abundance_metrics
                .iter()
                .all(|name| leaders[name].contains(&p.oligo))
        })
        .map(|p| p.oligo.clone())
        .collect();
    let abundance_candidate =
        (eligible.len() >= 2 && consensus.len() == 1).then(|| consensus[0].clone());
    let overlap_candidate =
        (leaders["barcode_jaccard"].len() == 1).then(|| leaders["barcode_jaccard"][0].clone());
    let overlap_agrees = abundance_candidate
        .as_ref()
        .map(|candidate| leaders["barcode_jaccard"].contains(candidate));
    let coverage_agrees = abundance_candidate
        .as_ref()
        .map(|candidate| leaders["cdna_read_coverage"].contains(candidate));
    let mut diagnostic_warnings = Vec::new();
    if overlap_agrees == Some(false) {
        diagnostic_warnings.push("jaccard_disagrees_with_abundance");
    }
    if coverage_agrees == Some(false) {
        diagnostic_warnings.push("read_coverage_disagrees_with_abundance");
    }
    let status = if eligible.is_empty() {
        "insufficient_data"
    } else if !eligible.iter().any(|p| p.shared_barcodes > 0) {
        "no_shared_barcodes"
    } else if eligible.len() == 1 {
        "single_candidate_unresolved"
    } else if consensus.is_empty() {
        "conflicting_abundance_metrics"
    } else if consensus.len() > 1 {
        "tied_candidates"
    } else if !eligible
        .iter()
        .any(|p| p.oligo == consensus[0] && p.sufficient_shared_barcode_count)
    {
        "insufficient_shared_barcodes"
    } else if eligible.iter().any(|p| p.partial_sampling) {
        "abundance_agree_sampled_unconfirmed"
    } else {
        "abundance_agree_unvalidated"
    };
    let preferred = if status == "abundance_agree_unvalidated" {
        consensus.first().cloned()
    } else {
        None
    };
    let expected = library_name(sample, OLIGO_ROLE);
    let expected_present = pairs
        .iter()
        .any(|p| p.cdna == cdna_name && p.oligo == expected);
    let expected_pair_assessment = match (expected_present, preferred.as_deref()) {
        (false, _) => "expected_candidate_missing",
        (true, Some(chosen)) if chosen == expected => "expected_pair_supported",
        (true, Some(_)) => "suspected_cross_sample_pairing",
        (true, None) => "unresolved",
    };
    Ranking {
        cdna: cdna_name.to_owned(),
        status,
        eligible_candidates: eligible.iter().map(|p| p.oligo.clone()).collect(),
        metric_orders: orders,
        metric_leaders: leaders,
        consensus_leaders: consensus,
        preferred_candidate: preferred,
        abundance_candidate,
        overlap_candidate,
        overlap_agrees,
        coverage_agrees,
        diagnostic_warnings,
        coverage_margin,
        informative_metrics,
        ranking_basis,
        count_cosine_margin,
        log_count_cosine_margin,
        expected_candidate: expected_present.then_some(expected),
        expected_pair_assessment,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::CDNA_ROLE;
    use std::collections::HashMap;
    fn rank(cdna: &Profile, pairs: &[Pair]) -> Ranking {
        super::rank(&cdna.name, &cdna.sample, pairs)
    }
    fn profile(name: &str, sample: &str, role: &'static str, items: &[(u64, u64)]) -> Profile {
        let counts: HashMap<_, _> = items.iter().copied().collect();
        Profile {
            name: name.into(),
            sample: sample.into(),
            role,
            files: vec![],
            observed_barcodes: counts.len(),
            supported_barcodes: counts.len(),
            supported_reads: counts.values().sum(),
            top: counts.keys().copied().collect(),
            norm: (counts
                .values()
                .map(|&n| (n as u128) * (n as u128))
                .sum::<u128>() as f64)
                .sqrt(),
            log_norm: crate::profile::log_summary(&counts).1,
            log_weights: crate::profile::log_summary(&counts).0,
            abundance_tiers: vec![],
            count_effective_barcodes: None,
            largest_barcode_read_fraction: None,
            counts,
        }
    }
    #[test]
    fn weighted_coverage_rank_and_ties() {
        let c = profile("c/cDNA", "c", CDNA_ROLE, &[(1, 90), (2, 10)]);
        let good = profile("good", "good", OLIGO_ROLE, &[(1, 20)]);
        let bad = profile("bad", "bad", OLIGO_ROLE, &[(2, 20)]);
        let pairs = vec![compare(&c, &bad, 1), compare(&c, &good, 1)];
        assert_eq!(pairs[1].cdna_read_coverage, Some(0.9));
        assert_eq!(pairs[1].barcode_jaccard, Some(0.5));
        assert!(!pairs[0].expected_pair);
        let ranking = rank(&c, &pairs);
        assert_eq!(ranking.preferred_candidate.as_deref(), Some("good"));
        assert!((ranking.coverage_margin.unwrap() - 0.8).abs() < 1e-12);
        let ties = vec![
            compare(&c, &good, 1),
            compare(
                &c,
                &profile("also_good", "also_good", OLIGO_ROLE, &[(1, 30)]),
                1,
            ),
        ];
        assert_eq!(rank(&c, &ties).consensus_leaders.len(), 2);
    }
    #[test]
    fn broad_background_is_not_automatically_selected() {
        let c = profile(
            "c/cDNA",
            "c",
            CDNA_ROLE,
            &(0..200).map(|i| (i, 3)).collect::<Vec<_>>(),
        );
        let expected = profile(
            "expected",
            "expected",
            OLIGO_ROLE,
            &(0..180).map(|i| (i, 3)).collect::<Vec<_>>(),
        );
        let broad = profile(
            "broad",
            "broad",
            OLIGO_ROLE,
            &(0..2000).map(|i| (i, 3)).collect::<Vec<_>>(),
        );
        let pairs = vec![compare(&c, &expected, 100), compare(&c, &broad, 100)];
        let result = rank(&c, &pairs);
        assert_eq!(result.status, "abundance_agree_unvalidated");
        assert_eq!(result.preferred_candidate.as_deref(), Some("expected"));
        assert_eq!(result.overlap_agrees, Some(true));
        assert_eq!(result.coverage_agrees, Some(false));
        assert_eq!(
            result.diagnostic_warnings,
            ["read_coverage_disagrees_with_abundance"]
        );
        assert_eq!(result.metric_leaders["cdna_read_coverage"], ["broad"]);
        assert_eq!(result.metric_leaders["barcode_jaccard"], ["expected"]);
    }
    #[test]
    fn partial_sampling_does_not_confirm_preference_even_when_metrics_agree() {
        let c = profile("c/cDNA", "c", CDNA_ROLE, &[(1, 3), (2, 3)]);
        let expected = profile("expected", "expected", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let wrong = profile("wrong", "wrong", OLIGO_ROLE, &[(3, 3), (4, 3)]);
        let mut pairs = vec![compare(&c, &expected, 2), compare(&c, &wrong, 2)];
        pairs[0].partial_sampling = true;
        let result = rank(&c, &pairs);
        assert_eq!(result.status, "abundance_agree_sampled_unconfirmed");
        assert!(result.preferred_candidate.is_none());
        assert_eq!(result.consensus_leaders, ["expected"]);
    }
    #[test]
    fn no_overlap_and_ties_never_select_a_preferred_candidate() {
        let c = profile("c/cDNA", "c", CDNA_ROLE, &[(1, 3), (2, 3)]);
        let a = profile("a", "a", OLIGO_ROLE, &[(3, 3), (4, 3)]);
        let b = profile("b", "b", OLIGO_ROLE, &[(5, 3), (6, 3)]);
        let no_overlap = rank(&c, &[compare(&c, &a, 2), compare(&c, &b, 2)]);
        assert_eq!(no_overlap.status, "no_shared_barcodes");
        assert!(no_overlap.preferred_candidate.is_none());
        let identical = profile("a", "a", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let also_identical = profile("b", "b", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let tied = rank(
            &c,
            &[compare(&c, &identical, 2), compare(&c, &also_identical, 2)],
        );
        assert_eq!(tied.status, "tied_candidates");
        assert!(tied.preferred_candidate.is_none());
    }
    #[test]
    fn directory_pairing_is_assessed_against_the_expected_oligo() {
        let c = profile("S1/cDNA", "S1", CDNA_ROLE, &[(1, 3), (2, 3)]);
        let own = profile("S1/oligo", "S1", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let other = profile("S2/oligo", "S2", OLIGO_ROLE, &[(3, 3), (4, 3)]);
        let pairs = vec![compare(&c, &own, 2), compare(&c, &other, 2)];
        assert!(pairs[0].expected_pair);
        assert!(!pairs[1].expected_pair);
        let ranking = rank(&c, &pairs);
        assert_eq!(ranking.expected_candidate.as_deref(), Some("S1/oligo"));
        assert_eq!(ranking.preferred_candidate.as_deref(), Some("S1/oligo"));
        assert_eq!(ranking.expected_pair_assessment, "expected_pair_supported");

        let swapped_own = profile("S1/oligo", "S1", OLIGO_ROLE, &[(3, 3), (4, 3)]);
        let swapped_other = profile("S2/oligo", "S2", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let swapped = vec![compare(&c, &swapped_own, 2), compare(&c, &swapped_other, 2)];
        let ranking = rank(&c, &swapped);
        assert_eq!(ranking.preferred_candidate.as_deref(), Some("S2/oligo"));
        assert_eq!(
            ranking.expected_pair_assessment,
            "suspected_cross_sample_pairing"
        );

        let missing = rank(&c, &[compare(&c, &other, 2)]);
        assert_eq!(missing.expected_candidate, None);
        assert_eq!(
            missing.expected_pair_assessment,
            "expected_candidate_missing"
        );
    }
    #[test]
    fn weak_consensus_is_rejected_and_shared_floor_is_inclusive() {
        let items = |start, count| (start..start + count).map(|i| (i, 3)).collect::<Vec<_>>();
        let c = profile("S1/cDNA", "S1", CDNA_ROLE, &items(0, 200));
        let wrong = profile("S2/oligo", "S2", OLIGO_ROLE, &items(1000, 200));
        for shared in [1, 99, 100] {
            let mut own_items = items(0, shared);
            own_items.extend(items(300, 200 - shared));
            let own = profile("S1/oligo", "S1", OLIGO_ROLE, &own_items);
            let pairs = [compare(&c, &own, 100), compare(&c, &wrong, 100)];
            let ranked = rank(&c, &pairs);
            assert_eq!(ranked.consensus_leaders, ["S1/oligo"]);
            if shared < 100 {
                assert_eq!(ranked.status, "insufficient_shared_barcodes");
                assert!(ranked.preferred_candidate.is_none());
                assert_eq!(ranked.expected_pair_assessment, "unresolved");
            } else {
                assert_eq!(ranked.status, "abundance_agree_unvalidated");
            }
            assert_eq!(rank(&c, &pairs[..1]).status, "single_candidate_unresolved");
        }
    }

    #[test]
    fn smaller_set_iteration_preserves_directional_read_coverages() {
        let c = profile("c/cDNA", "c", CDNA_ROLE, &[(1, 90), (2, 10), (3, 100)]);
        let o = profile("o/oligo", "o", OLIGO_ROLE, &[(1, 4), (4, 6)]);
        let pair = compare(&c, &o, 1);
        assert_eq!(pair.shared_barcodes, 1);
        assert_eq!(pair.cdna_read_coverage, Some(0.45));
        assert_eq!(pair.oligo_read_coverage, Some(0.4));
        assert_eq!(pair.barcode_jaccard, Some(0.25));
        let reversed = compare(&o, &c, 1);
        assert_eq!(reversed.cdna_read_coverage, pair.oligo_read_coverage);
        assert_eq!(reversed.oligo_read_coverage, pair.cdna_read_coverage);
        assert_eq!(reversed.count_cosine, pair.count_cosine);
    }

    #[test]
    fn identical_sets_are_ranked_by_abundance_and_flat_profiles_stay_tied() {
        let weights: Vec<_> = (0..120)
            .map(|i| (i, if i == 0 { 300 } else { 3 }))
            .collect();
        let c = profile("S/cDNA", "S", CDNA_ROLE, &weights);
        let good = profile("S/oligo", "S", OLIGO_ROLE, &weights);
        let mut swapped = weights.clone();
        swapped[0].1 = 3;
        swapped[1].1 = 300;
        let bad = profile("B/oligo", "B", OLIGO_ROLE, &swapped);
        let pairs = [compare(&c, &good, 100), compare(&c, &bad, 100)];
        assert!(pairs
            .iter()
            .all(|p| p.barcode_jaccard == Some(1.0) && p.cdna_read_coverage == Some(1.0)));
        let ranking = rank(&c, &pairs);
        assert_eq!(ranking.ranking_basis, "abundance_only");
        assert_eq!(
            ranking.informative_metrics,
            ["count_cosine", "log_count_cosine"]
        );
        assert_eq!(ranking.preferred_candidate.as_deref(), Some("S/oligo"));
        assert!(ranking.count_cosine_margin.unwrap() > 0.9);
        assert!(ranking.log_count_cosine_margin.unwrap() > 0.0);
        let flat: Vec<_> = (0..120).map(|i| (i, 3)).collect();
        let c = profile("S/cDNA", "S", CDNA_ROLE, &flat);
        let a = profile("S/oligo", "S", OLIGO_ROLE, &flat);
        let b = profile("B/oligo", "B", OLIGO_ROLE, &flat);
        let ranking = rank(&c, &[compare(&c, &a, 100), compare(&c, &b, 100)]);
        assert_eq!(ranking.ranking_basis, "no_discriminating_abundance_metrics");
        assert_eq!(ranking.status, "tied_candidates");
        assert!(ranking.preferred_candidate.is_none());
    }

    #[test]
    fn raw_and_log_abundance_disagreement_remains_unresolved() {
        let c = profile("S/cDNA", "S", CDNA_ROLE, &[(1, 100), (2, 10), (3, 1)]);
        let a = profile("S/oligo", "S", OLIGO_ROLE, &[(1, 1000), (2, 1), (3, 1)]);
        let b = profile("B/oligo", "B", OLIGO_ROLE, &[(1, 10), (2, 4), (3, 1)]);
        let ranking = rank(&c, &[compare(&c, &a, 1), compare(&c, &b, 1)]);
        assert_eq!(ranking.ranking_basis, "abundance_only");
        assert_eq!(ranking.metric_leaders["count_cosine"], ["S/oligo"]);
        assert_eq!(ranking.metric_leaders["log_count_cosine"], ["B/oligo"]);
        assert_eq!(ranking.status, "conflicting_abundance_metrics");
        assert!(ranking.preferred_candidate.is_none());
    }

    #[test]
    fn saved_deep_sequencing_metrics_keep_set_conflicts_as_diagnostics() {
        let fixture: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/deep-sequencing-pairs.json"))
                .unwrap();
        let pairs: Vec<Pair> = serde_json::from_value(fixture["pairs"].clone()).unwrap();
        assert_eq!(pairs.len(), 121);
        for suffix in 11..=21 {
            let sample = format!("deep_{suffix}");
            let name = format!("{sample}/cDNA");
            let expected = format!("{sample}/oligo");
            let ranked = super::rank(&name, &sample, &pairs);
            assert_eq!(ranked.status, "abundance_agree_unvalidated");
            assert_eq!(
                ranked.preferred_candidate.as_deref(),
                Some(expected.as_str())
            );
            assert_eq!(ranked.coverage_agrees, Some(true));
            assert_eq!(ranked.expected_pair_assessment, "expected_pair_supported");
            if [14, 16, 17, 18, 20, 21].contains(&suffix) {
                assert_eq!(ranked.overlap_agrees, Some(false));
                assert_eq!(ranked.overlap_candidate.as_deref(), Some("deep_13/oligo"));
                assert_eq!(
                    ranked.diagnostic_warnings,
                    ["jaccard_disagrees_with_abundance"]
                );
            } else {
                assert_eq!(ranked.overlap_agrees, Some(true));
                assert!(ranked.diagnostic_warnings.is_empty());
            }
        }
    }

    #[test]
    fn overlap_cannot_break_an_abundance_tie_and_tied_diagnostics_are_compatible() {
        let c = profile("S/cDNA", "S", CDNA_ROLE, &[(1, 3), (2, 3)]);
        let a = profile("S/oligo", "S", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let b = profile("B/oligo", "B", OLIGO_ROLE, &[(1, 3), (2, 3)]);
        let mut pairs = [compare(&c, &a, 1), compare(&c, &b, 1)];
        // Isolate ranking behavior: only the diagnostic values disagree.
        pairs[1].barcode_jaccard = Some(0.5);
        let ranked = rank(&c, &pairs);
        assert_eq!(ranked.status, "tied_candidates");
        assert!(ranked.abundance_candidate.is_none());
        assert_eq!(ranked.overlap_agrees, None);
        assert_eq!(ranked.coverage_agrees, None);
        assert!(ranked.preferred_candidate.is_none());
        pairs[1].count_cosine = Some(0.5);
        pairs[1].log_count_cosine = Some(0.5);
        pairs[1].barcode_jaccard = pairs[0].barcode_jaccard;
        let ranked = rank(&c, &pairs);
        assert_eq!(ranked.overlap_candidate, None);
        assert_eq!(ranked.overlap_agrees, Some(true));
        assert_eq!(ranked.coverage_agrees, Some(true));
        assert!(ranked.diagnostic_warnings.is_empty());
    }

    #[test]
    fn log_cosine_matches_direct_formula_and_keeps_empty_null() {
        let c = profile("c/cDNA", "c", CDNA_ROLE, &[(1, 100), (2, 1), (3, 1)]);
        let o = profile("o/oligo", "o", OLIGO_ROLE, &[(1, 1), (2, 100), (3, 1)]);
        let pair = compare(&c, &o, 1);
        let (hi, lo) = (101f64.ln(), 2f64.ln());
        let expected = (2.0 * hi * lo + lo * lo) / (hi * hi + 2.0 * lo * lo);
        assert!((pair.log_count_cosine.unwrap() - expected).abs() < 1e-9);
        let reverse = compare(&o, &c, 1);
        assert_eq!(pair.log_count_cosine, reverse.log_count_cosine);
        let scaled = profile(
            "scaled",
            "scaled",
            OLIGO_ROLE,
            &[(1, 1000), (2, 10), (3, 10)],
        );
        let scaled_pair = compare(&c, &scaled, 1);
        assert!((scaled_pair.count_cosine.unwrap() - 1.0).abs() < 1e-12);
        assert!(scaled_pair.log_count_cosine.unwrap() < 0.99);
        let empty = profile("empty", "empty", OLIGO_ROLE, &[]);
        assert_eq!(compare(&c, &empty, 1).log_count_cosine, None);
    }

    #[test]
    fn empty_and_insufficient_data_are_not_matches() {
        let empty = profile("empty/cDNA", "empty", CDNA_ROLE, &[]);
        let pair = compare(&empty, &empty, 1);
        assert_eq!(pair.cdna_read_coverage, None);
        assert_eq!(pair.count_cosine, None);
        assert_eq!(rank(&empty, &[pair]).status, "insufficient_data");
    }
}
