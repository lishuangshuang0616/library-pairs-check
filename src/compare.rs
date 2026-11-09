use crate::input::{library_name, OLIGO_ROLE};
use crate::profile::Profile;
use serde::Serialize;
use std::collections::BTreeMap;

fn fraction(n: f64, d: f64) -> Option<f64> {
    (d > 0.0).then(|| (n / d).clamp(0.0, 1.0))
}

#[derive(Serialize)]
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
    pub top_barcode_overlap: Option<f64>,
    pub sufficient_barcode_count: bool,
    pub partial_sampling: bool,
}

pub fn compare(cdna: &Profile, oligo: &Profile, min_distinct: usize) -> Pair {
    let mut shared = 0;
    let (mut cdna_reads, mut oligo_reads, mut dot) = (0u64, 0u64, 0u128);
    for (key, &left) in &cdna.counts {
        if let Some(&right) = oligo.counts.get(key) {
            shared += 1;
            cdna_reads += left;
            oligo_reads += right;
            dot += left as u128 * right as u128;
        }
    }
    Pair {
        cdna: cdna.name.clone(),
        oligo: oligo.name.clone(),
        cdna_sample: cdna.sample.clone(),
        oligo_sample: oligo.sample.clone(),
        expected_pair: cdna.sample == oligo.sample && cdna.role != oligo.role,
        shared_barcodes: shared,
        barcode_jaccard: fraction(
            shared as f64,
            (cdna.counts.len() + oligo.counts.len() - shared) as f64,
        ),
        cdna_barcode_coverage: fraction(shared as f64, cdna.counts.len() as f64),
        oligo_barcode_coverage: fraction(shared as f64, oligo.counts.len() as f64),
        cdna_read_coverage: fraction(cdna_reads as f64, cdna.supported_reads as f64),
        oligo_read_coverage: fraction(oligo_reads as f64, oligo.supported_reads as f64),
        count_cosine: fraction(dot as f64, cdna.norm * oligo.norm),
        top_barcode_overlap: fraction(
            cdna.top.intersection(&oligo.top).count() as f64,
            cdna.top.len().min(oligo.top.len()) as f64,
        ),
        partial_sampling: cdna
            .files
            .iter()
            .chain(&oligo.files)
            .any(|f| f.sampled_records < f.scanned_records || !f.input_fully_scanned),
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
    pub coverage_margin: Option<f64>,
    /// `<sample>/oligo` implied by the sample directory, when it exists.
    pub expected_candidate: Option<String>,
    pub expected_pair_assessment: &'static str,
}

type Metric = (&'static str, fn(&Pair) -> f64);

pub fn rank(cdna: &Profile, pairs: &[Pair]) -> Ranking {
    let mut eligible: Vec<_> = pairs
        .iter()
        .filter(|p| p.cdna == cdna.name && p.sufficient_barcode_count)
        .collect();
    eligible.sort_by(|a, b| a.oligo.cmp(&b.oligo));
    let mut orders = BTreeMap::new();
    let mut leaders = BTreeMap::new();
    let metrics: [Metric; 3] = [
        ("cdna_read_coverage", |p| {
            p.cdna_read_coverage.unwrap_or(0.0)
        }),
        ("barcode_jaccard", |p| p.barcode_jaccard.unwrap_or(0.0)),
        ("count_cosine", |p| p.count_cosine.unwrap_or(0.0)),
    ];
    let mut coverage_margin = None;
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
        orders.insert(name, sorted.iter().map(|p| p.oligo.clone()).collect());
        leaders.insert(name, tied);
    }
    let consensus: Vec<_> = eligible
        .iter()
        .filter(|p| leaders.values().all(|names| names.contains(&p.oligo)))
        .map(|p| p.oligo.clone())
        .collect();
    let status = if eligible.is_empty() {
        "insufficient_data"
    } else if !eligible.iter().any(|p| p.shared_barcodes > 0) {
        "no_shared_barcodes"
    } else if eligible.len() == 1 {
        "single_candidate_unresolved"
    } else if consensus.is_empty() {
        "conflicting_metrics"
    } else if consensus.len() > 1 {
        "tied_candidates"
    } else if eligible.iter().any(|p| p.partial_sampling) {
        "metrics_agree_sampled_unconfirmed"
    } else {
        "metrics_agree_unvalidated"
    };
    let preferred = if status == "metrics_agree_unvalidated" {
        consensus.first().cloned()
    } else {
        None
    };
    let expected = library_name(&cdna.sample, OLIGO_ROLE);
    let expected_present = pairs
        .iter()
        .any(|p| p.cdna == cdna.name && p.oligo == expected);
    let expected_pair_assessment = match (expected_present, preferred.as_deref()) {
        (false, _) => "expected_candidate_missing",
        (true, Some(chosen)) if chosen == expected => "expected_pair_supported",
        (true, Some(_)) => "suspected_cross_sample_pairing",
        (true, None) => "unresolved",
    };
    Ranking {
        cdna: cdna.name.clone(),
        status,
        eligible_candidates: eligible.iter().map(|p| p.oligo.clone()).collect(),
        metric_orders: orders,
        metric_leaders: leaders,
        consensus_leaders: consensus,
        preferred_candidate: preferred,
        coverage_margin,
        expected_candidate: expected_present.then_some(expected),
        expected_pair_assessment,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::CDNA_ROLE;
    use std::collections::HashMap;
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
        assert_eq!(result.status, "conflicting_metrics");
        assert!(result.preferred_candidate.is_none());
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
        assert_eq!(result.status, "metrics_agree_sampled_unconfirmed");
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
    fn empty_and_insufficient_data_are_not_matches() {
        let empty = profile("empty/cDNA", "empty", CDNA_ROLE, &[]);
        let pair = compare(&empty, &empty, 1);
        assert_eq!(pair.cdna_read_coverage, None);
        assert_eq!(pair.count_cosine, None);
        assert_eq!(rank(&empty, &[pair]).status, "insufficient_data");
    }
}
