use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

pub struct Scrambler {
    rng: ChaCha8Rng,
    lfsr_state: u64,
    seed: u64,
    position: u64,
}

impl Scrambler {
    const LFSR_POLYNOMIAL: u64 = 0x800000000000000D;
    const LFSR_INIT: u64 = 0xACE1ACE1ACE1ACE1;

    pub fn new(seed: u64) -> Self {
        Scrambler {
            rng: ChaCha8Rng::seed_from_u64(seed),
            lfsr_state: Self::LFSR_INIT ^ seed,
            seed,
            position: 0,
        }
    }

    pub fn reset(&mut self) {
        self.rng = ChaCha8Rng::seed_from_u64(self.seed);
        self.lfsr_state = Self::LFSR_INIT ^ self.seed;
        self.position = 0;
    }

    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    fn step_lfsr(&mut self) -> u8 {
        let feedback = (self.lfsr_state >> 63) & 1;
        self.lfsr_state <<= 1;
        if feedback == 1 {
            self.lfsr_state ^= Self::LFSR_POLYNOMIAL;
        }
        (self.lfsr_state & 0xFF) as u8
    }

    fn get_keystream_byte(&mut self) -> u8 {
        let lfsr_byte = self.step_lfsr();
        let random_byte: u8 = self.rng.gen();
        lfsr_byte ^ random_byte
    }

    pub fn scramble_byte(&mut self, input: u8) -> u8 {
        let key = self.get_keystream_byte();
        self.position += 1;
        input ^ key
    }

    pub fn descramble_byte(&mut self, input: u8) -> u8 {
        self.scramble_byte(input)
    }

    pub fn scramble(&mut self, data: &[u8]) -> Vec<u8> {
        data.iter().map(|&b| self.scramble_byte(b)).collect()
    }

    pub fn descramble(&mut self, data: &[u8]) -> Vec<u8> {
        self.scramble(data)
    }

    pub fn scramble_in_place(&mut self, data: &mut [u8]) {
        for b in data.iter_mut() {
            *b = self.scramble_byte(*b);
        }
    }

    pub fn descramble_in_place(&mut self, data: &mut [u8]) {
        self.scramble_in_place(data)
    }
}

pub struct AdaptiveScrambler {
    base_seed: u64,
    max_attempts: usize,
    gc_min: f64,
    gc_max: f64,
    max_homopolymer: usize,
}

impl AdaptiveScrambler {
    pub fn new(seed: u64) -> Self {
        AdaptiveScrambler {
            base_seed: seed,
            max_attempts: 1000,
            gc_min: 0.40,
            gc_max: 0.60,
            max_homopolymer: 3,
        }
    }

    pub fn with_constraints(
        seed: u64,
        gc_min: f64,
        gc_max: f64,
        max_homopolymer: usize,
    ) -> Self {
        AdaptiveScrambler {
            base_seed: seed,
            max_attempts: 1000,
            gc_min,
            gc_max,
            max_homopolymer,
        }
    }

    pub fn set_max_attempts(&mut self, attempts: usize) {
        self.max_attempts = attempts;
    }

    fn check_gc_content(&self, bases: &[u8]) -> bool {
        if bases.is_empty() {
            return true;
        }
        let gc_count = bases.iter().filter(|&&b| b == b'C' || b == b'G').count();
        let gc_ratio = gc_count as f64 / bases.len() as f64;
        gc_ratio >= self.gc_min && gc_ratio <= self.gc_max
    }

    fn check_homopolymer(&self, bases: &[u8]) -> bool {
        if bases.is_empty() {
            return true;
        }
        let mut current_base = bases[0];
        let mut current_run = 1;

        for &b in &bases[1..] {
            if b == current_base {
                current_run += 1;
                if current_run > self.max_homopolymer {
                    return false;
                }
            } else {
                current_base = b;
                current_run = 1;
            }
        }
        true
    }

    pub fn scramble_until_valid(
        &mut self,
        data: &[u8],
        to_bases_fn: impl Fn(&[u8]) -> Vec<u8>,
    ) -> (Vec<u8>, u64) {
        for attempt in 0..self.max_attempts {
            let seed = self.base_seed.wrapping_add(attempt as u64);
            let mut scrambler = Scrambler::new(seed);
            let scrambled = scrambler.scramble(data);
            let bases = to_bases_fn(&scrambled);

            if self.check_gc_content(&bases) && self.check_homopolymer(&bases) {
                return (scrambled, seed);
            }
        }

        panic!(
            "Failed to find valid scramble after {} attempts. Consider increasing max_attempts or relaxing constraints.",
            self.max_attempts
        );
    }

    pub fn verify_scramble(
        &self,
        data: &[u8],
        seed: u64,
        to_bases_fn: impl Fn(&[u8]) -> Vec<u8>,
    ) -> bool {
        let mut scrambler = Scrambler::new(seed);
        let scrambled = scrambler.scramble(data);
        let bases = to_bases_fn(&scrambled);
        self.check_gc_content(&bases) && self.check_homopolymer(&bases)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scrambler_symmetry() {
        let mut scrambler = Scrambler::new(42);
        let data: Vec<u8> = (0..255).collect();

        let scrambled = scrambler.scramble(&data);
        assert_ne!(scrambled, data);

        scrambler.reset();
        let descrambled = scrambler.descramble(&scrambled);
        assert_eq!(descrambled, data);
    }

    #[test]
    fn test_scrambler_different_seeds() {
        let data: Vec<u8> = (0..100).collect();

        let mut s1 = Scrambler::new(42);
        let mut s2 = Scrambler::new(43);

        let r1 = s1.scramble(&data);
        let r2 = s2.scramble(&data);

        assert_ne!(r1, r2);
    }

    #[test]
    fn test_scrambler_same_seed() {
        let data: Vec<u8> = (0..100).collect();

        let mut s1 = Scrambler::new(42);
        let mut s2 = Scrambler::new(42);

        let r1 = s1.scramble(&data);
        let r2 = s2.scramble(&data);

        assert_eq!(r1, r2);
    }

    #[test]
    fn test_scramble_in_place() {
        let mut data: Vec<u8> = (0..100).collect();
        let original = data.clone();

        let mut scrambler = Scrambler::new(42);
        scrambler.scramble_in_place(&mut data);
        assert_ne!(data, original);

        scrambler.reset();
        scrambler.descramble_in_place(&mut data);
        assert_eq!(data, original);
    }

    #[test]
    fn test_single_byte_scramble() {
        for byte in 0..=255u8 {
            let mut scrambler = Scrambler::new(42);
            let scrambled = scrambler.scramble_byte(byte);
            let mut scrambler2 = Scrambler::new(42);
            let descrambled = scrambler2.descramble_byte(scrambled);
            assert_eq!(byte, descrambled);
        }
    }

    #[test]
    fn test_lfsr_deterministic() {
        let mut s1 = Scrambler::new(42);
        let mut s2 = Scrambler::new(42);

        for _ in 0..100 {
            assert_eq!(s1.step_lfsr(), s2.step_lfsr());
        }
    }

    #[test]
    fn test_position_tracking() {
        let mut scrambler = Scrambler::new(42);
        assert_eq!(scrambler.position(), 0);

        let data: Vec<u8> = (0..50).collect();
        scrambler.scramble(&data);
        assert_eq!(scrambler.position(), 50);

        scrambler.reset();
        assert_eq!(scrambler.position(), 0);
    }

    fn dummy_to_bases(data: &[u8]) -> Vec<u8> {
        let mut result = Vec::new();
        for &b in data {
            for i in (0..4).rev() {
                let bits = (b >> (i * 2)) & 0b11;
                let base = match bits {
                    0b00 => b'A',
                    0b01 => b'T',
                    0b10 => b'C',
                    0b11 => b'G',
                    _ => unreachable!(),
                };
                result.push(base);
            }
        }
        result
    }

    #[test]
    fn test_gc_check() {
        let scrambler = AdaptiveScrambler::new(42);

        let good_gc = vec![b'A', b'T', b'C', b'G', b'A', b'T', b'C', b'G'];
        assert!(scrambler.check_gc_content(&good_gc));

        let bad_gc_low = vec![b'A', b'A', b'A', b'A', b'T', b'T', b'T', b'T'];
        assert!(!scrambler.check_gc_content(&bad_gc_low));

        let bad_gc_high = vec![b'G', b'G', b'G', b'G', b'C', b'C', b'C', b'C'];
        assert!(!scrambler.check_gc_content(&bad_gc_high));
    }

    #[test]
    fn test_homopolymer_check() {
        let scrambler = AdaptiveScrambler::new(42);

        let good = vec![b'A', b'T', b'C', b'G', b'A', b'T', b'C', b'G'];
        assert!(scrambler.check_homopolymer(&good));

        let bad = vec![b'A', b'A', b'A', b'A', b'T', b'C', b'G'];
        assert!(!scrambler.check_homopolymer(&bad));

        let edge_good = vec![b'A', b'A', b'A', b'T', b'C', b'G'];
        assert!(scrambler.check_homopolymer(&edge_good));
    }

    #[test]
    fn test_adaptive_scrambler() {
        let mut adaptive = AdaptiveScrambler::new(42);
        let data = vec![0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF];

        let (scrambled, seed_used) = adaptive.scramble_until_valid(&data, dummy_to_bases);

        assert!(adaptive.verify_scramble(&data, seed_used, dummy_to_bases));

        let mut scrambler = Scrambler::new(seed_used);
        let descrambled = scrambler.descramble(&scrambled);
        assert_eq!(descrambled, data);
    }

    #[test]
    fn test_adaptive_with_custom_constraints() {
        let mut adaptive = AdaptiveScrambler::with_constraints(42, 0.45, 0.55, 2);
        let data = vec![0xAA, 0xAA, 0x55, 0x55, 0xAA, 0xAA, 0x55, 0x55];

        let (scrambled, seed_used) = adaptive.scramble_until_valid(&data, dummy_to_bases);

        assert!(adaptive.verify_scramble(&data, seed_used, dummy_to_bases));
    }
}
