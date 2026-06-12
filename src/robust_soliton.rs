use rand::Rng;
use std::f64::consts::LN_2;

pub struct RobustSolitonDistribution {
    k: usize,
    c: f64,
    delta: f64,
    z: f64,
    probabilities: Vec<f64>,
    cumulative: Vec<f64>,
}

impl RobustSolitonDistribution {
    pub fn new(k: usize, c: f64, delta: f64) -> Self {
        assert!(k > 0, "k must be positive");
        assert!(c > 0.0, "c must be positive");
        assert!(delta > 0.0 && delta < 1.0, "delta must be in (0, 1)");

        let r = Self::compute_r(k, c, delta);
        let probabilities = Self::compute_probabilities(k, r, delta);
        let z = probabilities.iter().sum::<f64>();
        let normalized: Vec<f64> = probabilities.iter().map(|p| p / z).collect();

        let cumulative = Self::compute_cumulative(&normalized);

        RobustSolitonDistribution {
            k,
            c,
            delta,
            z,
            probabilities: normalized,
            cumulative,
        }
    }

    fn compute_r(k: usize, c: f64, delta: f64) -> f64 {
        let k_f = k as f64;
        c * (k_f.ln() / delta).sqrt() * k_f.sqrt()
    }

    fn compute_probabilities(k: usize, r: f64, delta: f64) -> Vec<f64> {
        let mut probabilities = Vec::with_capacity(k);
        let k_f = k as f64;

        probabilities.push(1.0 / k_f + (r / k_f).min(1.0));

        for i in 2..=(k as isize) {
            let i_f = i as f64;
            let rho = 1.0 / (i_f * (i_f - 1.0));

            let tau = if i <= (k_f / r - 1.0) as isize {
                r / (i_f * k_f)
            } else if i == (k_f / r) as isize {
                r * (r / delta).ln() / k_f
            } else {
                0.0
            };

            probabilities.push(rho + tau);
        }

        probabilities
    }

    fn compute_cumulative(probabilities: &[f64]) -> Vec<f64> {
        let mut cumulative = Vec::with_capacity(probabilities.len());
        let mut sum = 0.0;
        for &p in probabilities {
            sum += p;
            cumulative.push(sum);
        }
        cumulative
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let u: f64 = rng.gen();
        let idx = self.cumulative.binary_search_by(|&p| {
            if p < u {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });

        match idx {
            Ok(i) => i + 1,
            Err(i) => (i + 1).min(self.k),
        }
    }

    pub fn get_probability(&self, d: usize) -> Option<f64> {
        if d == 0 || d > self.k {
            None
        } else {
            Some(self.probabilities[d - 1])
        }
    }

    pub fn expected_degree(&self) -> f64 {
        self.probabilities
            .iter()
            .enumerate()
            .map(|(i, p)| (i + 1) as f64 * p)
            .sum()
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn failure_probability(&self, n: usize) -> f64 {
        let k_f = self.k as f64;
        let n_f = n as f64;
        let overhead = n_f - k_f;
        if overhead <= 0.0 {
            return 1.0;
        }
        self.delta * (-overhead * LN_2 / self.z).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;

    #[test]
    fn test_robust_soliton_creation() {
        let dist = RobustSolitonDistribution::new(100, 0.1, 0.01);
        assert_eq!(dist.k(), 100);
        assert!((dist.probabilities.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sampling_range() {
        let dist = RobustSolitonDistribution::new(100, 0.1, 0.01);
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        for _ in 0..1000 {
            let d = dist.sample(&mut rng);
            assert!(d >= 1 && d <= 100);
        }
    }

    #[test]
    fn test_expected_degree() {
        let dist = RobustSolitonDistribution::new(1000, 0.1, 0.01);
        let ed = dist.expected_degree();
        assert!(ed > 1.0 && ed < 10.0);
    }

    #[test]
    fn test_probability_at_degree_1() {
        let dist = RobustSolitonDistribution::new(1000, 0.1, 0.01);
        let p1 = dist.get_probability(1).unwrap();
        assert!(p1 > 0.0);
    }

    #[test]
    fn test_distribution_statistics() {
        let k = 1000;
        let dist = RobustSolitonDistribution::new(k, 0.1, 0.01);
        let mut rng = ChaCha20Rng::seed_from_u64(12345);
        let samples: Vec<usize> = (0..100_000).map(|_| dist.sample(&mut rng)).collect();

        let mean = samples.iter().sum::<usize>() as f64 / samples.len() as f64;
        let expected = dist.expected_degree();
        assert!((mean - expected).abs() / expected < 0.05);

        let count_d1 = samples.iter().filter(|&&d| d == 1).count();
        let ratio_d1 = count_d1 as f64 / samples.len() as f64;
        let p1 = dist.get_probability(1).unwrap();
        assert!((ratio_d1 - p1).abs() / p1 < 0.05);
    }
}
