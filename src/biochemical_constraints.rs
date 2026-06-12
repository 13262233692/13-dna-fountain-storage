#[derive(Debug, Clone, Copy)]
pub struct BiochemicalConstraints {
    pub gc_min: f64,
    pub gc_max: f64,
    pub max_homopolymer_run: usize,
}

impl Default for BiochemicalConstraints {
    fn default() -> Self {
        BiochemicalConstraints {
            gc_min: 0.40,
            gc_max: 0.60,
            max_homopolymer_run: 3,
        }
    }
}

impl BiochemicalConstraints {
    pub fn new(gc_min: f64, gc_max: f64, max_homopolymer: usize) -> Self {
        assert!(gc_min >= 0.0 && gc_min <= 1.0, "gc_min must be in [0, 1]");
        assert!(gc_max >= 0.0 && gc_max <= 1.0, "gc_max must be in [0, 1]");
        assert!(gc_min < gc_max, "gc_min must be less than gc_max");
        assert!(max_homopolymer >= 1, "max_homopolymer must be at least 1");

        BiochemicalConstraints {
            gc_min,
            gc_max,
            max_homopolymer_run: max_homopolymer,
        }
    }

    pub fn strict() -> Self {
        BiochemicalConstraints::new(0.45, 0.55, 2)
    }

    pub fn relaxed() -> Self {
        BiochemicalConstraints::new(0.30, 0.70, 5)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConstraintViolation {
    GcContentTooLow { actual: f64, min: f64 },
    GcContentTooHigh { actual: f64, max: f64 },
    HomopolymerRun { base: u8, position: usize, run_length: usize, max_allowed: usize },
    InvalidBase { base: u8, position: usize },
}

impl ConstraintViolation {
    pub fn description(&self) -> String {
        match self {
            ConstraintViolation::GcContentTooLow { actual, min } => {
                format!(
                    "GC content too low: {:.2}% (minimum {:.2}%)",
                    actual * 100.0,
                    min * 100.0
                )
            }
            ConstraintViolation::GcContentTooHigh { actual, max } => {
                format!(
                    "GC content too high: {:.2}% (maximum {:.2}%)",
                    actual * 100.0,
                    max * 100.0
                )
            }
            ConstraintViolation::HomopolymerRun {
                base,
                position,
                run_length,
                max_allowed,
            } => {
                format!(
                    "Homopolymer run of {} {}'s at position {} (maximum allowed {})",
                    run_length,
                    *base as char,
                    position,
                    max_allowed
                )
            }
            ConstraintViolation::InvalidBase { base, position } => {
                format!(
                    "Invalid base '{}' (0x{:02X}) at position {}",
                    *base as char,
                    base,
                    position
                )
            }
        }
    }
}

pub struct BiochemicalValidator {
    constraints: BiochemicalConstraints,
}

impl BiochemicalValidator {
    pub fn new(constraints: BiochemicalConstraints) -> Self {
        BiochemicalValidator { constraints }
    }

    pub fn constraints(&self) -> &BiochemicalConstraints {
        &self.constraints
    }

    fn is_valid_base(base: u8) -> bool {
        matches!(base, b'A' | b'T' | b'C' | b'G' | b'a' | b't' | b'c' | b'g')
    }

    fn normalize_base(base: u8) -> u8 {
        base.to_ascii_uppercase()
    }

    pub fn validate(&self, sequence: &[u8]) -> Result<(), Vec<ConstraintViolation>> {
        let mut violations = Vec::new();

        if sequence.is_empty() {
            return Ok(());
        }

        for (i, &base) in sequence.iter().enumerate() {
            if !Self::is_valid_base(base) {
                violations.push(ConstraintViolation::InvalidBase { base, position: i });
            }
        }

        if !violations.is_empty() {
            return Err(violations);
        }

        let normalized: Vec<u8> = sequence.iter().map(|&b| Self::normalize_base(b)).collect();

        if let Err(gc_violation) = self.validate_gc_content(&normalized) {
            violations.push(gc_violation);
        }

        violations.extend(self.validate_homopolymers(&normalized));

        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations)
        }
    }

    fn validate_gc_content(&self, sequence: &[u8]) -> Result<(), ConstraintViolation> {
        let gc_count = sequence
            .iter()
            .filter(|&&b| b == b'C' || b == b'G')
            .count();
        let gc_ratio = gc_count as f64 / sequence.len() as f64;

        if gc_ratio < self.constraints.gc_min {
            Err(ConstraintViolation::GcContentTooLow {
                actual: gc_ratio,
                min: self.constraints.gc_min,
            })
        } else if gc_ratio > self.constraints.gc_max {
            Err(ConstraintViolation::GcContentTooHigh {
                actual: gc_ratio,
                max: self.constraints.gc_max,
            })
        } else {
            Ok(())
        }
    }

    fn validate_homopolymers(&self, sequence: &[u8]) -> Vec<ConstraintViolation> {
        let mut violations = Vec::new();
        if sequence.is_empty() {
            return violations;
        }

        let mut current_base = sequence[0];
        let mut current_run_start = 0;
        let mut current_run_length = 1;

        for (i, &base) in sequence.iter().enumerate().skip(1) {
            if base == current_base {
                current_run_length += 1;
            } else {
                if current_run_length > self.constraints.max_homopolymer_run {
                    violations.push(ConstraintViolation::HomopolymerRun {
                        base: current_base,
                        position: current_run_start,
                        run_length: current_run_length,
                        max_allowed: self.constraints.max_homopolymer_run,
                    });
                }
                current_base = base;
                current_run_start = i;
                current_run_length = 1;
            }
        }

        if current_run_length > self.constraints.max_homopolymer_run {
            violations.push(ConstraintViolation::HomopolymerRun {
                base: current_base,
                position: current_run_start,
                run_length: current_run_length,
                max_allowed: self.constraints.max_homopolymer_run,
            });
        }

        violations
    }

    pub fn gc_content(&self, sequence: &[u8]) -> f64 {
        if sequence.is_empty() {
            return 0.0;
        }
        let gc_count = sequence
            .iter()
            .filter(|&&b| {
                let b = Self::normalize_base(b);
                b == b'C' || b == b'G'
            })
            .count();
        gc_count as f64 / sequence.len() as f64
    }

    pub fn longest_homopolymer_run(&self, sequence: &[u8]) -> (u8, usize, usize) {
        if sequence.is_empty() {
            return (0, 0, 0);
        }

        let mut longest_base = sequence[0];
        let mut longest_start = 0;
        let mut longest_length = 1;

        let mut current_base = sequence[0];
        let mut current_start = 0;
        let mut current_length = 1;

        for (i, &base) in sequence.iter().enumerate().skip(1) {
            if base == current_base {
                current_length += 1;
            } else {
                if current_length > longest_length {
                    longest_base = current_base;
                    longest_start = current_start;
                    longest_length = current_length;
                }
                current_base = base;
                current_start = i;
                current_length = 1;
            }
        }

        if current_length > longest_length {
            longest_base = current_base;
            longest_start = current_start;
            longest_length = current_length;
        }

        (longest_base, longest_start, longest_length)
    }

    pub fn validate_windowed(
        &self,
        sequence: &[u8],
        window_size: usize,
    ) -> Result<(), Vec<ConstraintViolation>> {
        let mut violations = Vec::new();

        if sequence.len() <= window_size {
            return self.validate(sequence);
        }

        for i in 0..=sequence.len().saturating_sub(window_size) {
            let window = &sequence[i..i + window_size];
            if let Err(window_violations) = self.validate(window) {
                for v in window_violations {
                    let adjusted_v = match v {
                        ConstraintViolation::HomopolymerRun {
                            base,
                            position,
                            run_length,
                            max_allowed,
                        } => ConstraintViolation::HomopolymerRun {
                            base,
                            position: position + i,
                            run_length,
                            max_allowed,
                        },
                        ConstraintViolation::InvalidBase { base, position } => {
                            ConstraintViolation::InvalidBase {
                                base,
                                position: position + i,
                            }
                        }
                        other => other,
                    };
                    if !violations.contains(&adjusted_v) {
                        violations.push(adjusted_v);
                    }
                }
            }
        }

        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_constraints() {
        let c = BiochemicalConstraints::default();
        assert_eq!(c.gc_min, 0.40);
        assert_eq!(c.gc_max, 0.60);
        assert_eq!(c.max_homopolymer_run, 3);
    }

    #[test]
    fn test_strict_constraints() {
        let c = BiochemicalConstraints::strict();
        assert_eq!(c.gc_min, 0.45);
        assert_eq!(c.gc_max, 0.55);
        assert_eq!(c.max_homopolymer_run, 2);
    }

    #[test]
    fn test_valid_sequence() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGATCGATCGATCG";
        assert!(validator.validate(seq).is_ok());
    }

    #[test]
    fn test_gc_content_too_low() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"AAAAAAAAAAAAATTTTTTTTTTTTT";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(matches!(
            violations[0],
            ConstraintViolation::GcContentTooLow { .. }
        ));
    }

    #[test]
    fn test_gc_content_too_high() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"CCCCCCCCCCCCGGGGGGGGGGGGG";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(matches!(
            violations[0],
            ConstraintViolation::GcContentTooHigh { .. }
        ));
    }

    #[test]
    fn test_homopolymer_violation() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGAAAAATCG";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(violations.iter().any(|v| matches!(
            v,
            ConstraintViolation::HomopolymerRun {
                base: b'A',
                position: 4,
                run_length: 5,
                max_allowed: 3
            }
        )));
    }

    #[test]
    fn test_homopolymer_edge_case() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGAAATCG";
        assert!(validator.validate(seq).is_ok());
    }

    #[test]
    fn test_invalid_base() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGXATCG";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(matches!(
            violations[0],
            ConstraintViolation::InvalidBase { base: b'X', position: 4 }
        ));
    }

    #[test]
    fn test_lowercase_bases() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"atcgatcgatcg";
        assert!(validator.validate(seq).is_ok());
    }

    #[test]
    fn test_mixed_case_with_violation() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGaaaaaTCG";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(violations.iter().any(|v| matches!(
            v,
            ConstraintViolation::HomopolymerRun {
                base: b'A',
                position: 4,
                run_length: 5,
                max_allowed: 3
            }
        )));
    }

    #[test]
    fn test_gc_content_calculation() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCG";
        assert!((validator.gc_content(seq) - 0.5).abs() < 1e-10);

        let seq2 = b"AAAATTTT";
        assert!((validator.gc_content(seq2) - 0.0).abs() < 1e-10);

        let seq3 = b"CCCCGGGG";
        assert!((validator.gc_content(seq3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_longest_homopolymer() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGGGGGATTTTA";
        let (base, pos, len) = validator.longest_homopolymer_run(seq);
        assert_eq!(base, b'G');
        assert_eq!(pos, 3);
        assert_eq!(len, 5);
    }

    #[test]
    fn test_multiple_violations() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"AAAAACCCCCGGGGGT";
        let result = validator.validate(seq);
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert!(violations.len() >= 3);
    }

    #[test]
    fn test_windowed_validation() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        let seq = b"ATCGATCGAAAAATCGATCG";
        let result = validator.validate_windowed(seq, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_constraint_violation_description() {
        let v = ConstraintViolation::GcContentTooLow {
            actual: 0.30,
            min: 0.40,
        };
        assert!(v.description().contains("30.00%"));

        let v = ConstraintViolation::HomopolymerRun {
            base: b'A',
            position: 10,
            run_length: 5,
            max_allowed: 3,
        };
        assert!(v.description().contains("5 A's"));
    }

    #[test]
    fn test_empty_sequence() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());
        assert!(validator.validate(&[]).is_ok());
        assert!((validator.gc_content(&[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_gc_content() {
        let validator = BiochemicalValidator::new(BiochemicalConstraints::default());

        let seq_40 = b"GGCCATATAT";
        assert!(validator.validate(seq_40).is_ok());

        let seq_50 = b"GCCGATATAT";
        assert!(validator.validate(seq_50).is_ok());

        let seq_60 = b"GCGCCGATAT";
        assert!(validator.validate(seq_60).is_ok());
    }
}
