use crate::robust_soliton::RobustSolitonDistribution;
use crc::{Crc, CRC_32_ISCSI};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub const CRC32: Crc<u32> = Crc::<u32>::new(&CRC_32_ISCSI);

#[derive(Debug, Clone)]
pub struct EncodedPacket {
    pub id: u64,
    pub block_indices: Vec<usize>,
    pub data: Vec<u8>,
    pub crc32: u32,
}

pub struct LTEncoder {
    blocks: Vec<Vec<u8>>,
    block_size: usize,
    k: usize,
    distribution: RobustSolitonDistribution,
    seed: u64,
    rng: ChaCha20Rng,
    packet_count: u64,
    original_size: usize,
}

impl LTEncoder {
    pub fn from_file<P: AsRef<Path>>(
        path: P,
        block_size: usize,
        c: f64,
        delta: f64,
        seed: u64,
    ) -> std::io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        let original_size = buffer.len();

        let blocks = Self::split_into_blocks(&buffer, block_size);
        let k = blocks.len();

        let distribution = RobustSolitonDistribution::new(k, c, delta);
        let rng = ChaCha20Rng::seed_from_u64(seed);

        Ok(LTEncoder {
            blocks,
            block_size,
            k,
            distribution,
            seed,
            rng,
            packet_count: 0,
            original_size,
        })
    }

    pub fn from_bytes(
        data: &[u8],
        block_size: usize,
        c: f64,
        delta: f64,
        seed: u64,
    ) -> Self {
        let original_size = data.len();
        let blocks = Self::split_into_blocks(data, block_size);
        let k = blocks.len();
        let distribution = RobustSolitonDistribution::new(k, c, delta);
        let rng = ChaCha20Rng::seed_from_u64(seed);

        LTEncoder {
            blocks,
            block_size,
            k,
            distribution,
            seed,
            rng,
            packet_count: 0,
            original_size,
        }
    }

    fn split_into_blocks(data: &[u8], block_size: usize) -> Vec<Vec<u8>> {
        let mut blocks = Vec::new();
        for chunk in data.chunks(block_size) {
            let mut block = chunk.to_vec();
            if block.len() < block_size {
                block.resize(block_size, 0);
            }
            blocks.push(block);
        }
        blocks
    }

    pub fn generate_packet(&mut self) -> EncodedPacket {
        let d = self.distribution.sample(&mut self.rng);
        let d = d.min(self.k);

        let mut indices: Vec<usize> = (0..self.k).collect();
        indices.shuffle(&mut self.rng);
        let mut selected_indices: Vec<usize> = indices.into_iter().take(d).collect();
        selected_indices.sort();

        let mut xor_result = vec![0u8; self.block_size];
        for &idx in &selected_indices {
            for (res_byte, &block_byte) in xor_result.iter_mut().zip(self.blocks[idx].iter()) {
                *res_byte ^= block_byte;
            }
        }

        let crc = CRC32.checksum(&xor_result);

        self.packet_count += 1;

        EncodedPacket {
            id: self.packet_count,
            block_indices: selected_indices,
            data: xor_result,
            crc32: crc,
        }
    }

    pub fn generate_packets(&mut self, count: usize) -> Vec<EncodedPacket> {
        (0..count).map(|_| self.generate_packet()).collect()
    }

    pub fn generate_packets_with_overhead(&mut self, overhead_factor: f64) -> Vec<EncodedPacket> {
        let count = (self.k as f64 * overhead_factor).ceil() as usize;
        self.generate_packets(count)
    }

    pub fn verify_packet(&self, packet: &EncodedPacket) -> bool {
        let crc = CRC32.checksum(&packet.data);
        crc == packet.crc32
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn block_size(&self) -> usize {
        self.block_size
    }

    pub fn original_size(&self) -> usize {
        self.original_size
    }

    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn packet_count(&self) -> u64 {
        self.packet_count
    }

    pub fn get_block(&self, index: usize) -> Option<&[u8]> {
        self.blocks.get(index).map(|v| v.as_slice())
    }

    pub fn failure_probability(&self) -> f64 {
        self.distribution.failure_probability(self.packet_count as usize)
    }
}

pub struct LTDecoder {
    block_size: usize,
    k: usize,
    original_size: usize,
    decoded_blocks: Vec<Option<Vec<u8>>>,
    pending_packets: Vec<EncodedPacket>,
    degrees: Vec<usize>,
}

impl LTDecoder {
    pub fn new(k: usize, block_size: usize, original_size: usize) -> Self {
        LTDecoder {
            block_size,
            k,
            original_size,
            decoded_blocks: vec![None; k],
            pending_packets: Vec::new(),
            degrees: Vec::new(),
        }
    }

    pub fn add_packet(&mut self, packet: EncodedPacket) -> bool {
        let crc = CRC32.checksum(&packet.data);
        if crc != packet.crc32 {
            return false;
        }

        self.pending_packets.push(packet);
        self.degrees.push(self.pending_packets.last().unwrap().block_indices.len());
        self.belief_propagation();
        true
    }

    fn belief_propagation(&mut self) {
        let mut changed = true;
        while changed {
            changed = false;

            let mut i = 0;
            while i < self.pending_packets.len() {
                let packet = &self.pending_packets[i];
                let mut unknown_indices = Vec::new();

                for &idx in &packet.block_indices {
                    if self.decoded_blocks[idx].is_none() {
                        unknown_indices.push(idx);
                    }
                }

                if unknown_indices.is_empty() {
                    self.pending_packets.remove(i);
                    self.degrees.remove(i);
                    changed = true;
                    continue;
                }

                if unknown_indices.len() == 1 {
                    let target_idx = unknown_indices[0];
                    let mut result = packet.data.clone();

                    for &idx in &packet.block_indices {
                        if idx != target_idx {
                            if let Some(block) = &self.decoded_blocks[idx] {
                                for (res_byte, &block_byte) in result.iter_mut().zip(block.iter()) {
                                    *res_byte ^= block_byte;
                                }
                            }
                        }
                    }

                    self.decoded_blocks[target_idx] = Some(result);
                    self.pending_packets.remove(i);
                    self.degrees.remove(i);
                    changed = true;

                    for j in 0..self.pending_packets.len() {
                        let indices = &self.pending_packets[j].block_indices;
                        if indices.contains(&target_idx) {
                            self.degrees[j] = self.degrees[j].saturating_sub(1);
                        }
                    }
                    break;
                } else {
                    i += 1;
                }
            }
        }
    }

    pub fn is_complete(&self) -> bool {
        self.decoded_blocks.iter().all(|b| b.is_some())
    }

    pub fn decoded_count(&self) -> usize {
        self.decoded_blocks.iter().filter(|b| b.is_some()).count()
    }

    pub fn pending_count(&self) -> usize {
        self.pending_packets.len()
    }

    pub fn get_data(&self) -> Option<Vec<u8>> {
        if !self.is_complete() {
            return None;
        }

        let mut result = Vec::with_capacity(self.k * self.block_size);
        for block in &self.decoded_blocks {
            if let Some(data) = block {
                result.extend_from_slice(data);
            }
        }
        result.truncate(self.original_size);
        Some(result)
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> std::io::Result<bool> {
        if let Some(data) = self.get_data() {
            let mut file = File::create(path)?;
            file.write_all(&data)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_block_splitting() {
        let data: Vec<u8> = (0..200).collect();
        let blocks = LTEncoder::split_into_blocks(&data, 50);
        assert_eq!(blocks.len(), 4);
        assert_eq!(blocks[0].len(), 50);
        assert_eq!(blocks[0][0], 0);
        assert_eq!(blocks[3][49], 199);
    }

    #[test]
    fn test_block_splitting_with_padding() {
        let data: Vec<u8> = (0..130).collect();
        let blocks = LTEncoder::split_into_blocks(&data, 50);
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[2].len(), 50);
        assert_eq!(blocks[2][29], 129);
        assert_eq!(blocks[2][30], 0);
    }

    #[test]
    fn test_encoder_creation_from_bytes() {
        let data: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();
        let encoder = LTEncoder::from_bytes(&data, 32, 0.1, 0.01, 42);
        assert_eq!(encoder.k(), 32);
        assert_eq!(encoder.block_size(), 32);
        assert_eq!(encoder.original_size(), 1000);
    }

    #[test]
    fn test_packet_generation() {
        let data: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();
        let mut encoder = LTEncoder::from_bytes(&data, 32, 0.1, 0.01, 42);

        let packet = encoder.generate_packet();
        assert_eq!(packet.id, 1);
        assert!(!packet.block_indices.is_empty());
        assert!(packet.block_indices.len() <= encoder.k());
        assert_eq!(packet.data.len(), 32);
        assert!(encoder.verify_packet(&packet));
    }

    #[test]
    fn test_multiple_packets() {
        let data: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();
        let mut encoder = LTEncoder::from_bytes(&data, 32, 0.1, 0.01, 42);

        let packets = encoder.generate_packets(10);
        assert_eq!(packets.len(), 10);
        assert_eq!(packets[9].id, 10);
        for packet in &packets {
            assert!(encoder.verify_packet(packet));
        }
    }

    #[test]
    fn test_crc_verification() {
        let data: Vec<u8> = (0..100).map(|i| (i % 256) as u8).collect();
        let mut encoder = LTEncoder::from_bytes(&data, 32, 0.1, 0.01, 42);

        let mut packet = encoder.generate_packet();
        assert!(encoder.verify_packet(&packet));
        packet.data[0] ^= 0xFF;
        assert!(!encoder.verify_packet(&packet));
    }

    #[test]
    fn test_full_encode_decode() {
        let mut rng = rand::thread_rng();
        let data: Vec<u8> = (0..1000).map(|_| rng.gen()).collect();

        let block_size = 32;
        let mut encoder = LTEncoder::from_bytes(&data, block_size, 0.1, 0.01, 12345);
        let k = encoder.k();

        let mut decoder = LTDecoder::new(k, block_size, data.len());

        let mut packets_processed = 0;
        for _ in 0..(k * 2) {
            let packet = encoder.generate_packet();
            packets_processed += 1;
            decoder.add_packet(packet);
            if decoder.is_complete() {
                break;
            }
        }

        assert!(decoder.is_complete(), "Decoder should complete with enough packets");
        assert!(packets_processed <= k * 2, "Should not need excessive packets");

        let decoded = decoder.get_data().unwrap();
        assert_eq!(decoded, data, "Decoded data should match original");
    }

    #[test]
    fn test_decode_with_loss() {
        let mut rng = rand::thread_rng();
        let data: Vec<u8> = (0..500).map(|_| rng.gen()).collect();

        let block_size = 32;
        let mut encoder = LTEncoder::from_bytes(&data, block_size, 0.1, 0.01, 12345);
        let k = encoder.k();

        let mut decoder = LTDecoder::new(k, block_size, data.len());

        let loss_rate = 0.4;
        let mut packets_processed = 0;
        let mut packets_lost = 0;

        for _ in 0..(k * 3) {
            let packet = encoder.generate_packet();
            if rng.gen::<f64>() < loss_rate {
                packets_lost += 1;
                continue;
            }
            packets_processed += 1;
            decoder.add_packet(packet);
            if decoder.is_complete() {
                break;
            }
        }

        assert!(decoder.is_complete(), 
            "Decoder should complete even with 40% loss (processed: {}, lost: {}, k: {})", 
            packets_processed, packets_lost, k);

        let decoded = decoder.get_data().unwrap();
        assert_eq!(decoded, data, "Decoded data should match original");
    }
}
