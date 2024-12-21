use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};

use anyhow::{anyhow, bail};
use itertools::Itertools;
use rv::prelude::*;

use crate::dmr::util::DmrInterval;
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::BorrowingMoniod;

#[derive(Debug, Default, Clone)]
pub(super) struct AggregatedCounts {
    mod_code_counts: HashMap<ModCodeRepr, usize>,
    pub(super) total: usize,
}

impl AggregatedCounts {
    pub(super) fn try_new(
        mod_code_counts: HashMap<ModCodeRepr, usize>,
        total: usize,
    ) -> anyhow::Result<Self> {
        let total_modification_counts = mod_code_counts.values().sum::<usize>();
        if total_modification_counts > total {
            bail!(
                "total modification counts ({total_modification_counts}) \
                 cannot be greater than total counts ({total})"
            )
        }
        Ok(Self { mod_code_counts, total })
    }

    fn get_canonical_counts(&self) -> usize {
        // safe because we check at creation, could be more careful if there
        // was a chance that &mut self was available.
        self.total - self.modified_counts()
    }

    pub(super) fn modified_counts(&self) -> usize {
        self.mod_code_counts.values().sum::<usize>()
    }

    fn combine(&self, other: &Self) -> Self {
        let total = self.total + other.total;
        let mut counts = self.mod_code_counts.clone();
        other.mod_code_counts.iter().for_each(|(mod_code, count)| {
            *counts.entry(*mod_code).or_insert(0) += *count;
        });

        Self { mod_code_counts: counts, total }
    }

    fn categorical_trials(
        &self,
        mod_codes_to_index: &HashMap<ModCodeRepr, usize>,
    ) -> anyhow::Result<Vec<usize>> {
        let mut trials = self
            .mod_code_counts
            .iter()
            .try_fold(Vec::new(), |mut acc, (code, count)| {
                let index = *mod_codes_to_index.get(code)?;
                let mut trials = vec![index; *count];
                acc.append(&mut trials);
                Some(acc)
            })
            .ok_or(anyhow!("failed to make categorical trials"))?;
        let canonical_count = self.get_canonical_counts();
        trials.append(&mut vec![0usize; canonical_count]);
        Ok(trials)
    }

    pub(super) fn string_counts(&self) -> String {
        if self.mod_code_counts.is_empty() {
            ".".to_string()
        } else {
            let csv = self
                .mod_code_counts
                .iter()
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
                .fold(String::new(), |mut acc, (code, count)| {
                    acc.push_str(&format!("{}:{},", code, count));
                    acc
                });
            csv.chars().into_iter().take(csv.len() - 1).collect()
        }
    }

    pub(super) fn string_percentages(&self) -> String {
        if self.mod_code_counts.is_empty() {
            ".".to_string()
        } else {
            let csv = self
                .mod_code_counts
                .iter()
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
                .fold(String::new(), |mut acc, (code, count)| {
                    let frac = *count as f32 / self.total as f32;
                    acc.push_str(&format!("{}:{:.2},", code, frac * 100f32));
                    acc
                });
            csv.chars().into_iter().take(csv.len() - 1).collect()
        }
    }

    pub(super) fn iter_mod_fractions(
        &self,
    ) -> impl Iterator<Item = (ModCodeRepr, f32)> + '_ {
        self.mod_code_counts
            .iter()
            .map(|(code, count)| (*code, *count as f32 / self.total as f32))
    }

    pub(super) fn pct_modified(&self) -> f32 {
        self.modified_counts() as f32 / self.total as f32
    }
}

impl BorrowingMoniod for AggregatedCounts {
    fn zero() -> Self {
        Self { mod_code_counts: HashMap::new(), total: 0usize }
    }

    fn op(self, other: &Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: &Self) {
        for (code, count) in other.mod_code_counts.iter() {
            *self.mod_code_counts.entry(*code).or_insert(0usize) += *count;
        }
        self.total += other.total;
    }

    fn len(&self) -> usize {
        todo!()
    }
}

impl Display for AggregatedCounts {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.string_counts())
    }
}

#[derive(Debug)]
pub(super) struct ModificationCounts {
    control_counts: AggregatedCounts,
    exp_counts: AggregatedCounts,
    interval: DmrInterval,
    pub(crate) score: f64,
}

impl ModificationCounts {
    pub(super) fn header(a_name: &str, b_name: &str) -> String {
        let mut s = [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            &format!("{a_name}_counts"),
            &format!("{a_name}_total"),
            &format!("{b_name}_counts"),
            &format!("{b_name}_total"),
            &format!("{a_name}_mod_percentages"),
            &format!("{b_name}_mod_percentages"),
            &format!("{a_name}_pct_modified"),
            &format!("{b_name}_pct_modified"),
        ]
        .join("\t");
        s.push('\n');
        s
    }

    pub(super) fn new(
        control_counts: AggregatedCounts,
        exp_counts: AggregatedCounts,
        interval: DmrInterval,
    ) -> anyhow::Result<Self> {
        let score = llk_ratio(&control_counts, &exp_counts)?;
        Ok(Self { control_counts, exp_counts, interval, score })
    }

    pub(super) fn to_row(&self) -> anyhow::Result<String> {
        let sep = '\t';
        let start = self.interval.start();
        let stop = self.interval.stop();
        let line = format!(
            "\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}\n\
        ",
            self.interval.chrom,
            start,
            stop,
            self.interval.name,
            self.score,
            self.interval.strand.to_string(),
            self.control_counts.string_counts(),
            self.control_counts.total,
            self.exp_counts.string_counts(),
            self.exp_counts.total,
            self.control_counts.string_percentages(),
            self.exp_counts.string_percentages(),
            self.control_counts.pct_modified(),
            self.exp_counts.pct_modified(),
        );
        Ok(line)
    }
}

fn dirichlet_llk(
    counts: &AggregatedCounts,
    prior: &Dirichlet,
    mod_codes_to_index: &HashMap<ModCodeRepr, usize>,
) -> anyhow::Result<f64> {
    // categorical outputs, die rolls, etc.
    let xs = counts.categorical_trials(&mod_codes_to_index)?;
    let control_data = DataOrSuffStat::Data(&xs);
    let control_posterior = prior.posterior(&control_data);
    Ok(control_posterior.ln_m(&control_data))
}

fn llk_dirichlet(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let mods_to_index = control_counts
        .mod_code_counts
        .keys()
        .chain(exp_counts.mod_code_counts.keys())
        .copied()
        .collect::<HashSet<ModCodeRepr>>()
        .into_iter()
        .sorted_by(|a, b| a.cmp(b))
        .enumerate()
        .map(|(i, c)| (c, i + 1))
        .collect::<HashMap<ModCodeRepr, usize>>();

    let k = mods_to_index.len() + 1;
    let prior = Dirichlet::jeffreys(k)?;
    let llk_control = dirichlet_llk(&control_counts, &prior, &mods_to_index)?;
    let llk_exp = dirichlet_llk(&exp_counts, &prior, &mods_to_index)?;

    let combined_counts = control_counts.combine(exp_counts);
    let llk_combined = dirichlet_llk(&combined_counts, &prior, &mods_to_index)?;

    Ok(llk_control + llk_exp - llk_combined)
}

fn counts_to_trials(count_methyl: usize, count_canonical: usize) -> Vec<bool> {
    let mut x = vec![true; count_methyl];
    let mut y = vec![false; count_canonical];
    x.append(&mut y);
    x
}

fn beta_llk(count_methyl: usize, count_canonical: usize) -> f64 {
    let trials = counts_to_trials(count_methyl, count_canonical);
    let data = BernoulliData::Data(&trials);
    let prior = Beta::jeffreys();
    let posterior = prior.posterior(&data);
    posterior.ln_m(&data)
}

fn llk_beta(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let all_mods = control_counts
        .mod_code_counts
        .keys()
        .copied()
        .chain(exp_counts.mod_code_counts.keys().copied())
        .collect::<HashSet<ModCodeRepr>>();
    if all_mods.len() != 1 {
        bail!("should have exactly one modification to use beta llk")
    }
    let raw_mod_code =
        all_mods.into_iter().take(1).collect::<Vec<ModCodeRepr>>()[0];

    let control_methyls =
        *control_counts.mod_code_counts.get(&raw_mod_code).unwrap_or(&0);
    let control_canonicals = control_counts.get_canonical_counts();

    let llk_control = beta_llk(control_methyls, control_canonicals);
    let exp_methyls =
        *exp_counts.mod_code_counts.get(&raw_mod_code).unwrap_or(&0);
    let exp_canonicals = exp_counts.get_canonical_counts();
    let llk_exp = beta_llk(exp_methyls, exp_canonicals);
    let llk_same = beta_llk(
        exp_methyls + control_methyls,
        exp_canonicals + control_canonicals,
    );

    Ok(llk_control + llk_exp - llk_same)
}

pub(super) fn llk_ratio(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let n_categories = std::cmp::max(
        control_counts.mod_code_counts.keys().len(),
        exp_counts.mod_code_counts.keys().len(),
    ) + 1; // plus 1 for canonical
    if n_categories < 2 {
        return Ok(0f64);
    }
    if n_categories == 2 {
        llk_beta(control_counts, exp_counts)
    } else {
        llk_dirichlet(control_counts, exp_counts)
    }
}

#[cfg(test)]
mod dmr_model_tests {
    use std::collections::HashMap;

    use itertools::Itertools;
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rv::dist::Categorical;
    use rv::prelude::{Bernoulli, Rv};

    use crate::dmr::llr_model::{llk_beta, llk_dirichlet, AggregatedCounts};
    use crate::mod_base_code::{
        ModCodeRepr, HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE,
    };

    fn methyl_sample(p: f64, n: usize, rng: &mut StdRng) -> AggregatedCounts {
        let mod_count = Bernoulli::new(p)
            .unwrap()
            .sample(n, rng)
            .into_iter()
            .filter(|b: &bool| *b)
            .count();
        let mod_code_counts = HashMap::from([('m'.into(), mod_count)]);
        AggregatedCounts::try_new(mod_code_counts, n).unwrap()
    }

    fn hydroxy_sample(
        alphas: &[f64],
        n: usize,
        rng: &mut StdRng,
    ) -> AggregatedCounts {
        let mods = [HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE];
        let counts = Categorical::new(alphas)
            .unwrap()
            .sample(n, rng)
            .into_iter()
            .filter_map(|x: usize| match x {
                0 => None,
                _ => Some(mods[x - 1]),
            })
            .collect::<Vec<ModCodeRepr>>();
        let counts = counts.into_iter().counts();
        AggregatedCounts::try_new(counts, n).unwrap()
    }

    #[test]
    fn test_beta_llk() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let control = methyl_sample(0.9, 1000, &mut rng);
        let exp = methyl_sample(0.1, 1000, &mut rng);
        let llk_a = llk_beta(&control, &exp).unwrap();
        let control = methyl_sample(0.9, 1000, &mut rng);
        let exp = methyl_sample(0.92, 1000, &mut rng);
        let llk_b = llk_beta(&control, &exp).unwrap();
        assert!(llk_a > llk_b);
        let control = methyl_sample(0.1, 1000, &mut rng);
        let exp = methyl_sample(0.12, 1000, &mut rng);
        let llk_c = llk_beta(&control, &exp).unwrap();
        assert!(llk_a > llk_c);
    }

    #[test]
    fn test_dir_llk() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let control = hydroxy_sample(&[0.1, 0.3, 0.6], 1000, &mut rng);
        let exp = hydroxy_sample(&[0.1, 0.6, 0.3], 1000, &mut rng);
        let llk_a = llk_dirichlet(&control, &exp).unwrap();
        let control = hydroxy_sample(&[0.1, 0.3, 0.6], 1000, &mut rng);
        let exp = hydroxy_sample(&[0.1, 0.4, 0.5], 1000, &mut rng);
        let llk_b = llk_dirichlet(&control, &exp).unwrap();
        assert!(llk_a > llk_b);
    }
}
