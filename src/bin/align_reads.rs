use std::io::*;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .enumerate()
        .map(|(i, line)| {
            let mut line = line.split('\t');
            let seq1: Vec<_> = line.next().unwrap().as_bytes().to_vec();
            let seq2: Vec<_> = line.next().unwrap().as_bytes().to_vec();
            (i, seq1, seq2)
        })
        .collect();
    println!("ID\tType\tDist\tNumOfGap\tTime");
    for (i, seq1, seq2) in reads {
        use clean_diff::*;
        // Affine gap
        let start = std::time::Instant::now();
        let (_, aln) = affine_gap::align(&seq1, &seq2, 2, -2, -8, -1);
        let time = (std::time::Instant::now() - start).as_millis();
        let (dist, gaps) = aln.dist_and_num_of_gaps();
        println!("{i}\tAffine\t{dist}\t{gaps}\t{time}");
        // Usual edit
        let start = std::time::Instant::now();
        let (_, aln) = edit_dist_usual::edit_dist(&seq1, &seq2);
        let time = (std::time::Instant::now() - start).as_millis();
        let (dist, gaps) = aln.dist_and_num_of_gaps();
        println!("{i}\tNaive\t{dist}\t{gaps}\t{time}");
        // Clean edit
        let start = std::time::Instant::now();
        let (_, aln) = edit_dist_usual_clean::edit_dist(&seq1, &seq2);
        let time = (std::time::Instant::now() - start).as_millis();
        let (dist, gaps) = aln.dist_and_num_of_gaps();
        println!("{i}\tNClean\t{dist}\t{gaps}\t{time}");
        // Diff
        let start = std::time::Instant::now();
        let (_, aln) = edit_dist_diff::edit_dist(&seq1, &seq2);
        let time = (std::time::Instant::now() - start).as_millis();
        let (dist, gaps) = aln.dist_and_num_of_gaps();
        println!("{i}\tDiff\t{dist}\t{gaps}\t{time}");
        // Clean diff
        let start = std::time::Instant::now();
        let (_, aln) = edit_dist_diff_clean::edit_dist(&seq1, &seq2);
        let time = (std::time::Instant::now() - start).as_millis();
        let (dist, gaps) = aln.dist_and_num_of_gaps();
        println!("{i}\tNClean\t{dist}\t{gaps}\t{time}");
    }
    Ok(())
}
