use clap::{Parser, Subcommand};
use dna_fountain_storage::{
    BiochemicalValidator, DnaConverter, DnaFountainDecoder, DnaFountainEncoder,
    encode_file_to_fasta, decode_fasta_to_file,
};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "dna-fountain",
    version = "0.1.0",
    author = "DNA Fountain Storage Project",
    about = "DNA Fountain Storage System using LT Codes with Robust Soliton Distribution"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Encode a binary file to DNA sequences in FASTA format
    Encode {
        /// Input binary file path
        #[arg(short, long)]
        input: PathBuf,

        /// Output FASTA file path
        #[arg(short, long)]
        output: PathBuf,

        /// Block size in bytes (default: 32)
        #[arg(short = 'b', long, default_value_t = 32)]
        block_size: usize,

        /// Overhead factor (default: 2.0 for 100% overhead)
        #[arg(short = 'v', long, default_value_t = 2.0)]
        overhead: f64,

        /// Robust Soliton parameter c (default: 0.1)
        #[arg(long, default_value_t = 0.1)]
        c: f64,

        /// Robust Soliton parameter delta (default: 0.01)
        #[arg(long, default_value_t = 0.01)]
        delta: f64,

        /// Random seed for LT code
        #[arg(long, default_value_t = 12345)]
        lt_seed: u64,

        /// Base seed for scrambling
        #[arg(long, default_value_t = 42)]
        scramble_seed: u64,
    },

    /// Decode DNA sequences from FASTA format to binary file
    Decode {
        /// Input FASTA file path
        #[arg(short, long)]
        input: PathBuf,

        /// Output binary file path
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Validate DNA sequences in FASTA format against biochemical constraints
    Validate {
        /// Input FASTA file path
        #[arg(short, long)]
        input: PathBuf,

        /// Minimum GC content (default: 0.40)
        #[arg(long, default_value_t = 0.40)]
        gc_min: f64,

        /// Maximum GC content (default: 0.60)
        #[arg(long, default_value_t = 0.60)]
        gc_max: f64,

        /// Maximum homopolymer run length (default: 3)
        #[arg(long, default_value_t = 3)]
        max_homopolymer: usize,
    },

    /// Analyze a file and show encoding statistics
    Analyze {
        /// Input file path (binary or FASTA)
        #[arg(short, long)]
        input: PathBuf,

        /// Block size in bytes (default: 32)
        #[arg(short = 'b', long, default_value_t = 32)]
        block_size: usize,

        /// Overhead factor (default: 2.0)
        #[arg(short = 'v', long, default_value_t = 2.0)]
        overhead: f64,
    },

    /// Run end-to-end test with random data
    Test {
        /// Data size in bytes (default: 1024)
        #[arg(short, long, default_value_t = 1024)]
        size: usize,

        /// Simulate packet loss rate (default: 0.4)
        #[arg(short = 'l', long, default_value_t = 0.4)]
        loss_rate: f64,

        /// Block size in bytes (default: 32)
        #[arg(short = 'b', long, default_value_t = 32)]
        block_size: usize,

        /// Overhead factor (default: 3.0)
        #[arg(short = 'v', long, default_value_t = 3.0)]
        overhead: f64,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Encode {
            input,
            output,
            block_size,
            overhead,
            c,
            delta,
            lt_seed,
            scramble_seed,
        } => {
            encode_command(input, output, block_size, overhead, c, delta, lt_seed, scramble_seed)
        }
        Commands::Decode { input, output } => decode_command(input, output),
        Commands::Validate {
            input,
            gc_min,
            gc_max,
            max_homopolymer,
        } => validate_command(input, gc_min, gc_max, max_homopolymer),
        Commands::Analyze {
            input,
            block_size,
            overhead,
        } => analyze_command(input, block_size, overhead),
        Commands::Test {
            size,
            loss_rate,
            block_size,
            overhead,
        } => test_command(size, loss_rate, block_size, overhead),
    }
}

fn encode_command(
    input: PathBuf,
    output: PathBuf,
    block_size: usize,
    overhead: f64,
    c: f64,
    delta: f64,
    lt_seed: u64,
    scramble_seed: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║           DNA Fountain Encoder                              ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let input_size = fs::metadata(&input)?.len() as usize;
    let k = (input_size + block_size - 1) / block_size;

    println!("Input file:    {}", input.display());
    println!("Output file:   {}", output.display());
    println!("File size:     {} bytes", input_size);
    println!("Block size:    {} bytes", block_size);
    println!("Blocks (k):    {}", k);
    println!("Overhead:      {:.2}x", overhead);
    println!("LT parameters: c={}, delta={}", c, delta);
    println!("Seeds:         lt_seed={}, scramble_seed={}", lt_seed, scramble_seed);
    println!();

    let start = Instant::now();

    let num_packets = encode_file_to_fasta(
        &input,
        &output,
        block_size,
        overhead,
        c,
        delta,
        lt_seed,
        scramble_seed,
    )?;

    let duration = start.elapsed();

    let output_size = fs::metadata(&output)?.len() as usize;

    println!("✓ Encoding complete!");
    println!();
    println!("Generated {} DNA packets", num_packets);
    println!("Output FASTA size: {} bytes", output_size);
    println!("Encoding time:     {:?}", duration);
    println!("Throughput:        {:.2} KB/s", input_size as f64 / 1024.0 / duration.as_secs_f64());
    println!();

    let sequences = DnaConverter::read_fasta(&output)?;
    let total_bases: usize = sequences.iter().map(|s| s.bases.len()).sum();

    let validator = BiochemicalValidator::new(Default::default());
    let mut violations = 0;
    let mut min_gc = 1.0f64;
    let mut max_gc = 0.0f64;

    for seq in &sequences {
        let gc = validator.gc_content(&seq.bases);
        min_gc = min_gc.min(gc);
        max_gc = max_gc.max(gc);

        if validator.validate(&seq.bases).is_err() {
            violations += 1;
        }
    }

    println!("DNA Sequence Statistics:");
    println!("  Total bases:       {}", total_bases);
    println!("  Bases per packet:  {}", sequences[0].bases.len());
    println!("  GC content range:  {:.2}% - {:.2}%", min_gc * 100.0, max_gc * 100.0);
    println!("  Violations:        {}", violations);
    println!();

    Ok(())
}

fn decode_command(input: PathBuf, output: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║           DNA Fountain Decoder                              ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    println!("Input FASTA:   {}", input.display());
    println!("Output file:   {}", output.display());
    println!();

    let sequences = DnaConverter::read_fasta(&input)?;
    let first = &sequences[0];

    println!("Detected parameters:");
    println!("  Original size: {} bytes", first.metadata.original_size);
    println!("  Block size:    {} bytes", first.metadata.block_size);
    println!("  Blocks (k):    {}", first.metadata.k);
    println!("  Total packets: {}", sequences.len());
    println!();

    let start = Instant::now();

    let success = decode_fasta_to_file(&input, &output)?;

    let duration = start.elapsed();

    if success {
        let output_size = fs::metadata(&output)?.len() as usize;
        println!("✓ Decoding complete!");
        println!();
        println!("Output file size: {} bytes", output_size);
        println!("Decoding time:    {:?}", duration);
        println!("Throughput:       {:.2} KB/s", output_size as f64 / 1024.0 / duration.as_secs_f64());
    } else {
        println!("✗ Decoding failed - not enough packets");
    }

    Ok(())
}

fn validate_command(
    input: PathBuf,
    gc_min: f64,
    gc_max: f64,
    max_homopolymer: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use dna_fountain_storage::BiochemicalConstraints;

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║     DNA Sequence Biochemical Constraint Validator           ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    println!("Input file:   {}", input.display());
    println!("Constraints:");
    println!("  GC content:       {:.2}% - {:.2}%", gc_min * 100.0, gc_max * 100.0);
    println!("  Max homopolymer:  {} consecutive bases", max_homopolymer);
    println!();

    let constraints = BiochemicalConstraints::new(gc_min, gc_max, max_homopolymer);
    let validator = BiochemicalValidator::new(constraints);
    let sequences = DnaConverter::read_fasta(&input)?;

    println!("Validating {} sequences...", sequences.len());
    println!();

    let mut total_violations = 0;
    let mut valid_count = 0;
    let mut min_gc = 1.0f64;
    let mut max_gc = 0.0f64;
    let mut max_homopolymer_found = 0;

    for (i, seq) in sequences.iter().enumerate() {
        let gc = validator.gc_content(&seq.bases);
        min_gc = min_gc.min(gc);
        max_gc = max_gc.max(gc);

        let (_, _, hp_len) = validator.longest_homopolymer_run(&seq.bases);
        max_homopolymer_found = max_homopolymer_found.max(hp_len);

        match validator.validate(&seq.bases) {
            Ok(_) => {
                valid_count += 1;
            }
            Err(violations) => {
                total_violations += violations.len();
                println!("Sequence {} (packet_id={}):", i + 1, seq.metadata.packet_id);
                for v in &violations {
                    println!("  ✗ {}", v.description());
                }
            }
        }
    }

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                      Validation Summary                      ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Valid sequences:    {:>38}  ║", format!("{}/{}", valid_count, sequences.len()));
    println!("║  Total violations:   {:>38}  ║", total_violations);
    println!("║  GC content range:   {:>38}  ║", format!("{:.2}% - {:.2}%", min_gc * 100.0, max_gc * 100.0));
    println!("║  Max homopolymer:    {:>38}  ║", max_homopolymer_found);
    println!("╚══════════════════════════════════════════════════════════════╝");

    if total_violations > 0 {
        std::process::exit(1);
    }

    Ok(())
}

fn analyze_command(
    input: PathBuf,
    block_size: usize,
    overhead: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                DNA Fountain File Analyzer                   ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let input_size = fs::metadata(&input)?.len() as usize;
    let k = (input_size + block_size - 1) / block_size;
    let num_packets = (k as f64 * overhead).ceil() as usize;

    println!("File:           {}", input.display());
    println!("Size:           {} bytes ({:.2} KB)", input_size, input_size as f64 / 1024.0);
    println!();

    println!("Encoding Parameters:");
    println!("  Block size:    {} bytes", block_size);
    println!("  Blocks (k):    {}", k);
    println!("  Overhead:      {:.2}x", overhead);
    println!("  Packets:       {} (for {:.2}x overhead)", num_packets, overhead);
    println!();

    let bases_per_packet = block_size * 4;
    let total_bases = num_packets * bases_per_packet;
    let expansion_ratio = total_bases as f64 / input_size as f64;

    println!("DNA Storage Estimates:");
    println!("  Bases per packet:  {}", bases_per_packet);
    println!("  Total bases:       {}", total_bases);
    println!("  Expansion ratio:   {:.2}x (bytes to bases)", expansion_ratio);
    println!();

    let bytes_per_mb = 1024 * 1024;
    let bases_per_gb = 1_000_000_000;
    let mb_input = input_size as f64 / bytes_per_mb as f64;
    let gb_bases = total_bases as f64 / bases_per_gb as f64;

    println!("For 1 GB input file:");
    let gb_k = (bytes_per_mb * 1024 + block_size - 1) / block_size;
    let gb_packets = (gb_k as f64 * overhead).ceil() as usize;
    let gb_total_bases = gb_packets * bases_per_packet;
    println!("  Blocks:       {}", gb_k);
    println!("  Packets:      {}", gb_packets);
    println!("  Total bases:  {} ({:.2} GB)", gb_total_bases, gb_total_bases as f64 / bases_per_gb as f64);
    println!();

    println!("Expected failure probability at {:.2}x overhead:", overhead);
    use dna_fountain_storage::RobustSolitonDistribution;
    let dist = RobustSolitonDistribution::new(k, 0.1, 0.01);
    let fp = dist.failure_probability(num_packets);
    println!("  P(failure) ≈ {:.2e}", fp);
    println!();

    if input.extension().and_then(|e| e.to_str()) == Some("fasta")
        || input.extension().and_then(|e| e.to_str()) == Some("fa")
    {
        let sequences = DnaConverter::read_fasta(&input)?;
        println!("FASTA File Analysis:");
        println!("  Number of sequences: {}", sequences.len());
        if !sequences.is_empty() {
            let total_bases: usize = sequences.iter().map(|s| s.bases.len()).sum();
            println!("  Total bases:         {}", total_bases);
            println!("  First packet ID:     {}", sequences[0].metadata.packet_id);
            println!("  k:                   {}", sequences[0].metadata.k);
            println!("  Block size:          {}", sequences[0].metadata.block_size);
            println!("  Original size:       {}", sequences[0].metadata.original_size);
        }
    }

    Ok(())
}

fn test_command(
    size: usize,
    loss_rate: f64,
    block_size: usize,
    overhead: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    use rand::Rng;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║         DNA Fountain End-to-End Test Suite                  ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    println!("Test parameters:");
    println!("  Data size:      {} bytes", size);
    println!("  Block size:     {} bytes", block_size);
    println!("  Overhead:       {:.2}x", overhead);
    println!("  Loss rate:      {:.0}%", loss_rate * 100.0);
    println!();

    let mut rng = ChaCha20Rng::seed_from_u64(12345);
    let original_data: Vec<u8> = (0..size).map(|_| rng.gen()).collect();

    println!("Step 1: Encoding data to DNA packets...");
    let start = Instant::now();

    let mut encoder = DnaFountainEncoder::from_bytes(
        &original_data,
        block_size,
        0.1,
        0.01,
        12345,
        42,
    );

    let num_packets = (encoder.k() as f64 * overhead).ceil() as usize;
    let dna_sequences = encoder.generate_dna_packets(num_packets)?;

    let encode_time = start.elapsed();
    println!("✓ Generated {} DNA packets in {:?}", dna_sequences.len(), encode_time);
    println!();

    println!("Step 2: Validating biochemical constraints...");
    let validator = BiochemicalValidator::new(Default::default());
    let mut all_valid = true;
    let mut min_gc = 1.0f64;
    let mut max_gc = 0.0f64;

    for seq in &dna_sequences {
        let gc = validator.gc_content(&seq.bases);
        min_gc = min_gc.min(gc);
        max_gc = max_gc.max(gc);

        if let Err(violations) = validator.validate(&seq.bases) {
            all_valid = false;
            println!("  ✗ Packet {} has {} violations", seq.metadata.packet_id, violations.len());
        }
    }

    if all_valid {
        println!("✓ All packets pass biochemical constraints");
        println!("  GC content range: {:.2}% - {:.2}%", min_gc * 100.0, max_gc * 100.0);
    } else {
        println!("✗ Some packets failed validation");
    }
    println!();

    println!("Step 3: Simulating {:.0}% packet loss...", loss_rate * 100.0);
    let mut rng_loss = ChaCha20Rng::seed_from_u64(67890);
    let mut received_packets = Vec::new();
    let mut lost_count = 0;

    for seq in &dna_sequences {
        if rng_loss.gen::<f64>() >= loss_rate {
            received_packets.push(seq.clone());
        } else {
            lost_count += 1;
        }
    }

    println!("✓ Received {}/{} packets ({} lost)",
        received_packets.len(),
        dna_sequences.len(),
        lost_count
    );
    println!();

    println!("Step 4: Decoding received packets...");
    let start = Instant::now();

    let mut decoder = DnaFountainDecoder::new(
        encoder.k(),
        encoder.block_size(),
        encoder.original_size(),
    );

    let mut packets_processed = 0;
    for seq in &received_packets {
        packets_processed += 1;
        decoder.add_dna_sequence(seq)?;
        if decoder.is_complete() {
            break;
        }
    }

    let decode_time = start.elapsed();

    if decoder.is_complete() {
        println!("✓ Decoding successful!");
        println!("  Processed {} packets", packets_processed);
        println!("  Time: {:?}", decode_time);
    } else {
        println!("✗ Decoding incomplete");
        println!("  Decoded {}/{} blocks", decoder.decoded_count(), encoder.k());
        return Err("Decoding failed".into());
    }
    println!();

    println!("Step 5: Verifying data integrity...");
    let decoded_data = decoder.get_data().unwrap();

    if decoded_data == original_data {
        println!("✓ Data integrity verified - original and decoded data match perfectly!");
    } else {
        println!("✗ Data mismatch detected");
        let mismatches: Vec<usize> = original_data
            .iter()
            .zip(decoded_data.iter())
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(i, _)| i)
            .collect();
        println!("  Mismatches at positions: {:?}", &mismatches[..mismatches.len().min(10)]);
        return Err("Data integrity check failed".into());
    }
    println!();

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                    TEST PASSED ✓                            ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Original size:  {:>38} bytes  ║", size);
    println!("║  Packets sent:   {:>38}  ║", dna_sequences.len());
    println!("║  Packets lost:   {:>38}  ║", lost_count);
    println!("║  Packets used:   {:>38}  ║", packets_processed);
    println!("║  GC range:       {:>38}  ║", format!("{:.2}% - {:.2}%", min_gc * 100.0, max_gc * 100.0));
    println!("║  Encoding time:  {:>38?}  ║", encode_time);
    println!("║  Decoding time:  {:>38?}  ║", decode_time);
    println!("╚══════════════════════════════════════════════════════════════╝");

    Ok(())
}
