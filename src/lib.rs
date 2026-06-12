pub mod robust_soliton;
pub mod lt_encoder;
pub mod scrambler;
pub mod biochemical_constraints;
pub mod dna_converter;

pub use robust_soliton::RobustSolitonDistribution;
pub use lt_encoder::{EncodedPacket, LTEncoder, LTDecoder, CRC32};
pub use scrambler::{AdaptiveScrambler, Scrambler};
pub use biochemical_constraints::{
    BiochemicalConstraints, BiochemicalValidator, ConstraintViolation,
};
pub use dna_converter::{DnaConverter, DnaMetadata, DnaSequence};

use std::path::Path;

pub struct DnaFountainEncoder {
    lt_encoder: LTEncoder,
    dna_converter: DnaConverter,
    pub scramble_base_seed: u64,
}

impl DnaFountainEncoder {
    pub fn from_file<P: AsRef<Path>>(
        input_path: P,
        block_size: usize,
        c: f64,
        delta: f64,
        lt_seed: u64,
        scramble_seed: u64,
    ) -> std::io::Result<Self> {
        let lt_encoder = LTEncoder::from_file(input_path, block_size, c, delta, lt_seed)?;
        let dna_converter = DnaConverter::new(scramble_seed);

        Ok(DnaFountainEncoder {
            lt_encoder,
            dna_converter,
            scramble_base_seed: scramble_seed,
        })
    }

    pub fn from_bytes(
        data: &[u8],
        block_size: usize,
        c: f64,
        delta: f64,
        lt_seed: u64,
        scramble_seed: u64,
    ) -> Self {
        let lt_encoder = LTEncoder::from_bytes(data, block_size, c, delta, lt_seed);
        let dna_converter = DnaConverter::new(scramble_seed);

        DnaFountainEncoder {
            lt_encoder,
            dna_converter,
            scramble_base_seed: scramble_seed,
        }
    }

    pub fn generate_dna_packet(&mut self) -> Result<DnaSequence, String> {
        let packet = self.lt_encoder.generate_packet();
        let packet_seed = self.scramble_base_seed.wrapping_add(packet.id);
        self.dna_converter = DnaConverter::new(packet_seed);

        let metadata = DnaMetadata {
            original_size: self.lt_encoder.original_size(),
            block_size: self.lt_encoder.block_size(),
            k: self.lt_encoder.k(),
            seed: self.lt_encoder.seed(),
            scramble_seed: packet_seed,
            packet_id: packet.id,
            block_indices: packet.block_indices.clone(),
            crc32: packet.crc32,
        };

        self.dna_converter
            .encode_with_scramble(&packet.data, &metadata)
            .map_err(|e| format!("Biochemical constraint violations: {:?}", e))
    }

    pub fn generate_dna_packets(&mut self, count: usize) -> Result<Vec<DnaSequence>, String> {
        (0..count)
            .map(|_| self.generate_dna_packet())
            .collect()
    }

    pub fn generate_dna_packets_with_overhead(
        &mut self,
        overhead_factor: f64,
    ) -> Result<Vec<DnaSequence>, String> {
        let count = (self.lt_encoder.k() as f64 * overhead_factor).ceil() as usize;
        self.generate_dna_packets(count)
    }

    pub fn write_fasta<P: AsRef<Path>>(
        &mut self,
        output_path: P,
        overhead_factor: f64,
        header: Option<&str>,
    ) -> Result<usize, String> {
        let sequences = self.generate_dna_packets_with_overhead(overhead_factor)?;
        let count = sequences.len();
        self.dna_converter
            .write_fasta(&sequences, output_path, header)
            .map_err(|e| e.to_string())?;
        Ok(count)
    }

    pub fn k(&self) -> usize {
        self.lt_encoder.k()
    }

    pub fn block_size(&self) -> usize {
        self.lt_encoder.block_size()
    }

    pub fn original_size(&self) -> usize {
        self.lt_encoder.original_size()
    }

    pub fn packet_count(&self) -> u64 {
        self.lt_encoder.packet_count()
    }

    pub fn failure_probability(&self) -> f64 {
        self.lt_encoder.failure_probability()
    }
}

pub struct DnaFountainDecoder {
    decoder: LTDecoder,
    dna_converter: DnaConverter,
}

impl DnaFountainDecoder {
    pub fn new(k: usize, block_size: usize, original_size: usize) -> Self {
        DnaFountainDecoder {
            decoder: LTDecoder::new(k, block_size, original_size),
            dna_converter: DnaConverter::new(0),
        }
    }

    pub fn add_dna_sequence(&mut self, sequence: &DnaSequence) -> Result<bool, String> {
        let packet_data = self.dna_converter.decode(sequence)?;

        let packet = EncodedPacket {
            id: sequence.metadata.packet_id,
            block_indices: sequence.metadata.block_indices.clone(),
            data: packet_data,
            crc32: sequence.metadata.crc32,
        };

        Ok(self.decoder.add_packet(packet))
    }

    pub fn add_dna_sequences(&mut self, sequences: &[DnaSequence]) -> Result<usize, String> {
        let mut added = 0;
        for seq in sequences {
            if self.add_dna_sequence(seq)? {
                added += 1;
            }
        }
        Ok(added)
    }

    pub fn read_and_decode_fasta<P: AsRef<Path>>(&mut self, path: P) -> Result<usize, String> {
        let sequences = DnaConverter::read_fasta(path)?;
        self.add_dna_sequences(&sequences)
    }

    pub fn is_complete(&self) -> bool {
        self.decoder.is_complete()
    }

    pub fn decoded_count(&self) -> usize {
        self.decoder.decoded_count()
    }

    pub fn pending_count(&self) -> usize {
        self.decoder.pending_count()
    }

    pub fn get_data(&self) -> Option<Vec<u8>> {
        self.decoder.get_data()
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> std::io::Result<bool> {
        self.decoder.save_to_file(path)
    }
}

pub fn encode_file_to_fasta<P: AsRef<Path>>(
    input_path: P,
    output_path: P,
    block_size: usize,
    overhead_factor: f64,
    c: f64,
    delta: f64,
    lt_seed: u64,
    scramble_seed: u64,
) -> Result<usize, String> {
    let mut encoder = DnaFountainEncoder::from_file(
        input_path,
        block_size,
        c,
        delta,
        lt_seed,
        scramble_seed,
    )
    .map_err(|e| e.to_string())?;

    let header = format!(
        "DNA Fountain Encoding | original_size={} bytes | k={} blocks | block_size={} bytes | overhead={:.2}x",
        encoder.original_size(),
        encoder.k(),
        encoder.block_size(),
        overhead_factor
    );

    encoder.write_fasta(output_path, overhead_factor, Some(&header))
}

pub fn decode_fasta_to_file<P: AsRef<Path>>(
    fasta_path: P,
    output_path: P,
) -> Result<bool, String> {
    let sequences = DnaConverter::read_fasta(&fasta_path)?;

    if sequences.is_empty() {
        return Err("No sequences found in FASTA file".to_string());
    }

    let first = &sequences[0];
    let mut decoder = DnaFountainDecoder::new(
        first.metadata.k,
        first.metadata.block_size,
        first.metadata.original_size,
    );

    for seq in &sequences {
        decoder.add_dna_sequence(seq)?;
        if decoder.is_complete() {
            break;
        }
    }

    if !decoder.is_complete() {
        return Err(format!(
            "Not enough packets to decode. Decoded {}/{} blocks, need more packets.",
            decoder.decoded_count(),
            first.metadata.k
        ));
    }

    decoder.save_to_file(output_path).map_err(|e| e.to_string())
}
