use crate::biochemical_constraints::{BiochemicalConstraints, BiochemicalValidator, ConstraintViolation};
use crate::scrambler::{AdaptiveScrambler, Scrambler};
use std::fmt;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct DnaSequence {
    pub bases: Vec<u8>,
    pub metadata: DnaMetadata,
}

#[derive(Debug, Clone)]
pub struct DnaMetadata {
    pub original_size: usize,
    pub block_size: usize,
    pub k: usize,
    pub seed: u64,
    pub scramble_seed: u64,
    pub packet_id: u64,
    pub block_indices: Vec<usize>,
    pub crc32: u32,
}

impl fmt::Display for DnaMetadata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let indices_str = self
            .block_indices
            .iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join("; ");
        write!(
            f,
            "packet_id={},original_size={},block_size={},k={},seed={},scramble_seed={},indices=[{}],crc32={:08X}",
            self.packet_id,
            self.original_size,
            self.block_size,
            self.k,
            self.seed,
            self.scramble_seed,
            indices_str,
            self.crc32
        )
    }
}

pub struct DnaConverter {
    validator: BiochemicalValidator,
    adaptive_scrambler: AdaptiveScrambler,
}

impl DnaConverter {
    pub fn new(seed: u64) -> Self {
        let constraints = BiochemicalConstraints::default();
        DnaConverter {
            validator: BiochemicalValidator::new(constraints),
            adaptive_scrambler: AdaptiveScrambler::new(seed),
        }
    }

    pub fn with_constraints(seed: u64, constraints: BiochemicalConstraints) -> Self {
        DnaConverter {
            validator: BiochemicalValidator::new(constraints),
            adaptive_scrambler: AdaptiveScrambler::with_constraints(
                seed,
                constraints.gc_min,
                constraints.gc_max,
                constraints.max_homopolymer_run,
            ),
        }
    }

    pub fn bytes_to_bases(data: &[u8]) -> Vec<u8> {
        let mut bases = Vec::with_capacity(data.len() * 4);
        for &byte in data {
            for i in (0..4).rev() {
                let bits = (byte >> (i * 2)) & 0b11;
                let base = match bits {
                    0b00 => b'A',
                    0b01 => b'T',
                    0b10 => b'C',
                    0b11 => b'G',
                    _ => unreachable!(),
                };
                bases.push(base);
            }
        }
        bases
    }

    pub fn bases_to_bytes(bases: &[u8]) -> Result<Vec<u8>, String> {
        if bases.len() % 4 != 0 {
            return Err(format!(
                "Base sequence length must be multiple of 4, got {}",
                bases.len()
            ));
        }

        let mut bytes = Vec::with_capacity(bases.len() / 4);
        for chunk in bases.chunks(4) {
            let mut byte: u8 = 0;
            for (i, &base) in chunk.iter().enumerate() {
                let bits = match base {
                    b'A' | b'a' => 0b00,
                    b'T' | b't' => 0b01,
                    b'C' | b'c' => 0b10,
                    b'G' | b'g' => 0b11,
                    _ => return Err(format!("Invalid base character: {}", base as char)),
                };
                byte |= bits << ((3 - i) * 2);
            }
            bytes.push(byte);
        }
        Ok(bytes)
    }

    pub fn encode_with_scramble(
        &mut self,
        data: &[u8],
        metadata: &DnaMetadata,
    ) -> Result<DnaSequence, Vec<ConstraintViolation>> {
        let (scrambled_data, scramble_seed) = self
            .adaptive_scrambler
            .scramble_until_valid(data, Self::bytes_to_bases);

        let bases = Self::bytes_to_bases(&scrambled_data);

        if let Err(violations) = self.validator.validate(&bases) {
            return Err(violations);
        }

        let mut metadata = metadata.clone();
        metadata.scramble_seed = scramble_seed;

        Ok(DnaSequence { bases, metadata })
    }

    pub fn decode(&self, sequence: &DnaSequence) -> Result<Vec<u8>, String> {
        let scrambled_data = Self::bases_to_bytes(&sequence.bases)?;
        let mut scrambler = Scrambler::new(sequence.metadata.scramble_seed);
        let data = scrambler.descramble(&scrambled_data);
        Ok(data)
    }

    pub fn validate(&self, bases: &[u8]) -> Result<(), Vec<ConstraintViolation>> {
        self.validator.validate(bases)
    }

    pub fn gc_content(&self, bases: &[u8]) -> f64 {
        self.validator.gc_content(bases)
    }

    pub fn to_fasta(&self, sequences: &[DnaSequence], header: Option<&str>) -> String {
        let mut fasta = String::new();

        if let Some(h) = header {
            fasta.push_str(&format!("; {}\n", h));
        }

        for seq in sequences {
            fasta.push_str(&format!(
                "> packet_{} | {}\n",
                seq.metadata.packet_id, seq.metadata
            ));

            for chunk in seq.bases.chunks(80) {
                fasta.push_str(&format!("{}\n", String::from_utf8_lossy(chunk)));
            }
        }

        fasta
    }

    pub fn write_fasta<P: AsRef<Path>>(
        &self,
        sequences: &[DnaSequence],
        path: P,
        header: Option<&str>,
    ) -> std::io::Result<()> {
        let fasta = self.to_fasta(sequences, header);
        let mut file = File::create(path)?;
        file.write_all(fasta.as_bytes())?;
        Ok(())
    }

    pub fn parse_fasta(s: &str) -> Result<Vec<DnaSequence>, String> {
        let mut sequences = Vec::new();
        let mut current_header: Option<String> = None;
        let mut current_bases = Vec::new();

        for line in s.lines() {
            let line = line.trim();

            if line.is_empty() || line.starts_with(';') {
                continue;
            }

            if line.starts_with('>') {
                if let Some(header) = current_header.take() {
                    if !current_bases.is_empty() {
                        let metadata = Self::parse_metadata(&header)?;
                        sequences.push(DnaSequence {
                            bases: current_bases.clone(),
                            metadata,
                        });
                        current_bases.clear();
                    }
                }
                current_header = Some(line[1..].trim().to_string());
            } else {
                current_bases.extend(line.as_bytes().iter().filter(|&&b| !b.is_ascii_whitespace()));
            }
        }

        if let Some(header) = current_header {
            if !current_bases.is_empty() {
                let metadata = Self::parse_metadata(&header)?;
                sequences.push(DnaSequence {
                    bases: current_bases,
                    metadata,
                });
            }
        }

        Ok(sequences)
    }

    fn parse_metadata(header: &str) -> Result<DnaMetadata, String> {
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() < 2 {
            return Err("Invalid FASTA header format".to_string());
        }

        let metadata_str = parts[1].trim();
        let mut metadata = DnaMetadata {
            original_size: 0,
            block_size: 0,
            k: 0,
            seed: 0,
            scramble_seed: 0,
            packet_id: 0,
            block_indices: Vec::new(),
            crc32: 0,
        };

        for pair in metadata_str.split(',') {
            let pair = pair.trim();
            if let Some((key, value)) = pair.split_once('=') {
                match key.trim() {
                    "packet_id" => {
                        metadata.packet_id = value
                            .parse()
                            .map_err(|e| format!("Invalid packet_id: {}", e))?
                    }
                    "original_size" => {
                        metadata.original_size = value
                            .parse()
                            .map_err(|e| format!("Invalid original_size: {}", e))?
                    }
                    "block_size" => {
                        metadata.block_size = value
                            .parse()
                            .map_err(|e| format!("Invalid block_size: {}", e))?
                    }
                    "k" => metadata.k = value.parse().map_err(|e| format!("Invalid k: {}", e))?,
                    "seed" => {
                        metadata.seed = value
                            .parse()
                            .map_err(|e| format!("Invalid seed: {}", e))?
                    }
                    "scramble_seed" => {
                        metadata.scramble_seed = value
                            .parse()
                            .map_err(|e| format!("Invalid scramble_seed: {}", e))?
                    }
                    "indices" => {
                        let indices_str = value.trim_start_matches('[').trim_end_matches(']');
                        if !indices_str.is_empty() {
                            metadata.block_indices = indices_str
                                .split(|c| c == ';' || c == ',')
                                .map(|s| s.trim().parse())
                                .collect::<Result<Vec<usize>, _>>()
                                .map_err(|e| format!("Invalid indices: {}", e))?;
                        }
                    }
                    "crc32" => {
                        metadata.crc32 = u32::from_str_radix(value.trim_start_matches("0x"), 16)
                            .map_err(|e| format!("Invalid crc32: {}", e))?
                    }
                    _ => {}
                }
            }
        }

        Ok(metadata)
    }

    pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<DnaSequence>, String> {
        let mut file = File::open(path).map_err(|e| e.to_string())?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .map_err(|e| e.to_string())?;
        Self::parse_fasta(&contents)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bytes_to_bases() {
        let data = vec![0b00011011];
        let bases = DnaConverter::bytes_to_bases(&data);
        assert_eq!(bases, vec![b'A', b'T', b'C', b'G']);
    }

    #[test]
    fn test_bases_to_bytes() {
        let bases = vec![b'A', b'T', b'C', b'G'];
        let bytes = DnaConverter::bases_to_bytes(&bases).unwrap();
        assert_eq!(bytes, vec![0b00011011]);
    }

    #[test]
    fn test_bidirectional_conversion() {
        let data: Vec<u8> = (0..255).collect();
        let bases = DnaConverter::bytes_to_bases(&data);
        let decoded = DnaConverter::bases_to_bytes(&bases).unwrap();
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_invalid_base() {
        let bases = vec![b'A', b'T', b'X', b'G'];
        let result = DnaConverter::bases_to_bytes(&bases);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_length() {
        let bases = vec![b'A', b'T', b'C'];
        let result = DnaConverter::bases_to_bytes(&bases);
        assert!(result.is_err());
    }

    #[test]
    fn test_lowercase_bases() {
        let bases = vec![b'a', b't', b'c', b'g'];
        let bytes = DnaConverter::bases_to_bytes(&bases).unwrap();
        assert_eq!(bytes, vec![0b00011011]);
    }

    #[test]
    fn test_encode_with_scramble() {
        let mut converter = DnaConverter::new(42);
        let data = vec![0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF];

        let metadata = DnaMetadata {
            original_size: data.len(),
            block_size: 32,
            k: 10,
            seed: 12345,
            scramble_seed: 0,
            packet_id: 1,
            block_indices: vec![0, 2, 5],
            crc32: 0xDEADBEEF,
        };

        let sequence = converter.encode_with_scramble(&data, &metadata).unwrap();

        converter.validate(&sequence.bases).unwrap();

        let decoded = converter.decode(&sequence).unwrap();
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_fasta_formatting() {
        let converter = DnaConverter::new(42);

        let seq1 = DnaSequence {
            bases: vec![b'A'; 100],
            metadata: DnaMetadata {
                original_size: 25,
                block_size: 32,
                k: 10,
                seed: 12345,
                scramble_seed: 42,
                packet_id: 1,
                block_indices: vec![0, 2],
                crc32: 0xDEADBEEF,
            },
        };

        let seq2 = DnaSequence {
            bases: vec![b'T'; 50],
            metadata: DnaMetadata {
                original_size: 12,
                block_size: 32,
                k: 10,
                seed: 12345,
                scramble_seed: 43,
                packet_id: 2,
                block_indices: vec![1, 3, 5],
                crc32: 0xCAFEBABE,
            },
        };

        let fasta = converter.to_fasta(&[seq1, seq2], Some("Test FASTA"));

        assert!(fasta.contains("; Test FASTA"));
        assert!(fasta.contains("> packet_1"));
        assert!(fasta.contains("> packet_2"));
        assert!(fasta.contains("packet_id=1"));
        assert!(fasta.contains("packet_id=2"));
        assert!(fasta.contains("crc32=DEADBEEF"));
        assert!(fasta.contains("crc32=CAFEBABE"));
    }

    #[test]
    fn test_fasta_roundtrip() {
        let converter = DnaConverter::new(42);

        let sequences = vec![DnaSequence {
            bases: DnaConverter::bytes_to_bases(&[0xAA, 0xBB, 0xCC, 0xDD]),
            metadata: DnaMetadata {
                original_size: 4,
                block_size: 32,
                k: 10,
                seed: 12345,
                scramble_seed: 42,
                packet_id: 42,
                block_indices: vec![0, 1, 2, 3],
                crc32: 0x12345678,
            },
        }];

        let fasta = converter.to_fasta(&sequences, None);
        let parsed = DnaConverter::parse_fasta(&fasta).unwrap();

        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].bases, sequences[0].bases);
        assert_eq!(parsed[0].metadata.packet_id, 42);
        assert_eq!(parsed[0].metadata.original_size, 4);
        assert_eq!(parsed[0].metadata.scramble_seed, 42);
        assert_eq!(parsed[0].metadata.block_indices, vec![0, 1, 2, 3]);
        assert_eq!(parsed[0].metadata.crc32, 0x12345678);
    }

    #[test]
    fn test_metadata_display() {
        let metadata = DnaMetadata {
            original_size: 1000,
            block_size: 32,
            k: 32,
            seed: 12345,
            scramble_seed: 67890,
            packet_id: 42,
            block_indices: vec![0, 5, 10],
            crc32: 0xDEADBEEF,
        };

        let s = format!("{}", metadata);
        assert!(s.contains("packet_id=42"));
        assert!(s.contains("original_size=1000"));
        assert!(s.contains("crc32=DEADBEEF"));
        assert!(s.contains("indices=[0; 5; 10]"));
    }

    #[test]
    fn test_gc_content_encoding() {
        let mut converter = DnaConverter::new(42);
        let data = vec![0xFF; 100];

        let metadata = DnaMetadata {
            original_size: data.len(),
            block_size: 32,
            k: 4,
            seed: 12345,
            scramble_seed: 0,
            packet_id: 1,
            block_indices: vec![0],
            crc32: 0,
        };

        let sequence = converter.encode_with_scramble(&data, &metadata).unwrap();
        let gc = converter.gc_content(&sequence.bases);

        assert!(gc >= 0.40 && gc <= 0.60, "GC content was {:.2}%", gc * 100.0);
    }

    #[test]
    fn test_fasta_file_io() {
        let tmp_dir = std::env::temp_dir();
        let tmp_path = tmp_dir.join("test_dna.fasta");

        let converter = DnaConverter::new(42);
        let sequences = vec![DnaSequence {
            bases: DnaConverter::bytes_to_bases(&[0x01, 0x02, 0x03, 0x04]),
            metadata: DnaMetadata {
                original_size: 4,
                block_size: 32,
                k: 1,
                seed: 12345,
                scramble_seed: 42,
                packet_id: 1,
                block_indices: vec![0],
                crc32: 0,
            },
        }];

        converter.write_fasta(&sequences, &tmp_path, Some("Test")).unwrap();
        let loaded = DnaConverter::read_fasta(&tmp_path).unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].bases, sequences[0].bases);

        std::fs::remove_file(&tmp_path).ok();
    }

    #[test]
    fn test_fasta_with_newlines() {
        let fasta = "> packet_1 | packet_id=1,original_size=4,k=1,seed=1,scramble_seed=1,indices=[0],crc32=00000000\nATCG\nATCG\n> packet_2 | packet_id=2,original_size=4,k=1,seed=1,scramble_seed=1,indices=[0],crc32=00000000\nCGTA\nCGTA\n";

        let sequences = DnaConverter::parse_fasta(fasta).unwrap();
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0].bases, vec![b'A', b'T', b'C', b'G', b'A', b'T', b'C', b'G']);
        assert_eq!(sequences[1].bases, vec![b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'A']);
    }
}
