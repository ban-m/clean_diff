use clap;
use clap::Parser;
#[derive(Parser)]
#[clap(author, version = "0.1", about = "Simulate reads to be aligned. The output format is <ID>\t<Seq1>\t<Seq2>", long_about = None)]
struct Args {
    /// Numbers of reads to be simulated.
    #[clap(short, long, value_parser, value_name = "NUM", default_value_t = 100)]
    num_of_reads: usize,
    /// Error rate of a read. Should be smaller than 0.25.
    #[clap(short, long, value_parser, value_name = "ERR", default_value_t = 0.1)]
    error_rate: f64,
    /// Length of the template
    #[clap(short, long, value_parser, value_name = "LEN", default_value_t = 500)]
    length: usize,
    /// Seeds for the pseudorandom number generator.
    #[clap(short, long, value_parser, value_name = "SEED", default_value_t = 42)]
    seed: u64,
    /// Wether or not to use three-state Markov model. In the gap state, the prob to see gap is doubled.
    #[clap(short, long)]
    to_use_hmm: bool,
}

use kiley::gen_seq::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(args.seed);
    let error_rate = args.error_rate / 3f64;
    let profile = Profile::new(error_rate, error_rate, error_rate);
    let mat = (1f64 - 2f64 * error_rate, error_rate, error_rate);
    let gap = (
        1f64 - 4f64 * error_rate,
        2f64 * error_rate,
        2f64 * error_rate,
    );
    let mism_prob = (error_rate / (1f64 - 2f64 * error_rate)) / 3f64;
    let match_prob = 1f64 - 3f64 * mism_prob;
    let match_prob = vec![
        vec![match_prob, mism_prob, mism_prob, mism_prob],
        vec![mism_prob, match_prob, mism_prob, mism_prob],
        vec![mism_prob, mism_prob, match_prob, mism_prob],
        vec![mism_prob, mism_prob, mism_prob, match_prob],
    ]
    .concat();
    let hmm =
        kiley::hmm::guided::PairHiddenMarkovModel::new(mat, gap, gap, &match_prob, &[0.25; 20]);
    for id in 0..args.num_of_reads {
        let template = kiley::gen_seq::generate_seq(&mut rng, args.length);
        let mutated = if args.to_use_hmm {
            hmm.gen(&template, &mut rng)
        } else {
            introduce_randomness(&template, &mut rng, &profile)
        };
        let template = String::from_utf8_lossy(&template);
        let mutated = String::from_utf8_lossy(&mutated);
        println!("{id}\t{template}\t{mutated}");
    }
    Ok(())
}
