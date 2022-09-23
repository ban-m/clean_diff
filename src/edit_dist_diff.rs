use super::alignments::*;
/// Usual edit distance alignments and its path.
pub fn edit_dist(xs: &[u8], ys: &[u8]) -> (u32, Alignment) {
    if xs == ys {
        return (0, Alignment::new(vec![Op::Match; xs.len()]));
    }
    // h -> d = the furthest reaching point of the d - h diagonal with edit distace h, and the traceback pointer.
    let mut dp: Vec<Vec<(usize, Option<u8>)>> = vec![vec![(match_len(xs, 0, ys, 0), None)]];
    'outer: for dist in 1..xs.len() + ys.len() + 1 {
        let mut new_d = Vec::with_capacity(2 * dist + 1);
        for diag in 0..2 * dist + 1 {
            let from_above = dp[dist - 1].get(diag).map(|&(x, _)| x);
            let from_mat = match diag < 1 {
                true => None,
                false => dp[dist - 1].get(diag - 1).map(|(x, _)| x + 1),
            };
            let from_left = match diag < 2 {
                true => None,
                false => dp[dist - 1].get(diag - 2).map(|(x, _)| x + 1),
            };
            let (max_reach, trace) = max_three(from_above, from_mat, from_left);
            let i = (max_reach + dist) - diag;
            let snake = match_len(xs, i, ys, max_reach);
            new_d.push((max_reach + snake, trace));
            let j = max_reach + snake;
            let i = (j + dist) - diag;
            if i == xs.len() && j == ys.len() {
                dp.push(new_d);
                break 'outer;
            }
        }
        dp.push(new_d);
        assert!(dist < xs.len() + ys.len() + 1);
    }
    let opt_dist = dp.len() - 1;
    let (mut diag, (mut ypos, mut prev)) = dp.last().unwrap().iter().enumerate().last().unwrap();
    let mut ops = vec![];
    let mut dist = opt_dist;
    while let Some(trace) = prev {
        let old_ypos = ypos;
        let op = match trace {
            0 => Op::Del,
            1 => {
                diag -= 1;
                Op::Mismatch
            }
            2 => {
                diag -= 2;
                Op::Ins
            }
            _ => panic!(),
        };
        dist -= 1;
        (ypos, prev) = dp[dist][diag];
        let len = match op {
            Op::Del => old_ypos - ypos,
            Op::Mismatch | Op::Ins => old_ypos - ypos - 1,
            _ => panic!(),
        };
        ops.extend(std::iter::repeat(Op::Match).take(len));
        ops.push(op);
    }
    assert_eq!(dist, 0);
    let xpos = (ypos + dist) - diag;
    assert_eq!(xpos, ypos);
    ops.extend(std::iter::repeat(Op::Match).take(xpos));
    ops.reverse();
    let aln = Alignment::new(ops);
    (opt_dist as u32, aln)
}

fn match_len(xs: &[u8], x_start: usize, ys: &[u8], y_start: usize) -> usize {
    let xs = xs.iter().skip(x_start);
    let ys = ys.iter().skip(y_start);
    std::iter::zip(xs, ys).take_while(|(x, y)| x == y).count()
}

fn max_three(x: Option<usize>, y: Option<usize>, z: Option<usize>) -> (usize, Option<u8>) {
    let xy = match (x, y) {
        (None, None) => None,
        (Some(x), None) => Some((x, Some(0))),
        (None, Some(y)) => Some((y, Some(1))),
        (Some(x), Some(y)) if x < y => Some((y, Some(1))),
        (Some(x), Some(_)) => Some((x, Some(0))),
    };
    match (xy, z) {
        (None, Some(z)) => (z, Some(2)),
        (Some(xy), None) => xy,
        (Some((pos, Some(tr))), Some(z)) if z < pos => (pos, Some(tr)),
        (Some(_), Some(z)) => (z, Some(2)),
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::rand_core::SeedableRng;
    #[test]
    fn edit_dist_calc() {
        let xs = b"AAACCC";
        let (dist, ops) = edit_dist(xs, xs);
        assert_eq!(dist, 0);
        assert_eq!(ops.ops, vec![Op::Match; xs.len()]);
        let ys = b"AAAGCC";
        let (dist, aln) = edit_dist(xs, ys);
        assert_eq!(dist, 1);
        assert_eq!(
            aln.ops,
            vec![vec![Op::Match; 3], vec![Op::Mismatch], vec![Op::Match; 2]].concat()
        );
        let ys = b"AAA";
        let (dist, aln) = edit_dist(xs, ys);
        assert_eq!(dist, 3);
        assert_eq!(aln.ops, vec![vec![Op::Match; 3], vec![Op::Del; 3]].concat());
        let ys = b"";
        let (dist, aln) = edit_dist(xs, ys);
        assert_eq!(dist, 6);
        assert_eq!(aln.ops, vec![Op::Del; 6]);
        let (dist, aln) = edit_dist(ys, xs);
        assert_eq!(dist, 6);
        assert_eq!(aln.ops, vec![Op::Ins; 6]);
        let xs = b"ACGT";
        let ys = b"ACCTG";
        let (dist, aln) = edit_dist(xs, ys);
        assert_eq!(dist, 2);
        assert_eq!(
            aln.ops,
            vec![Op::Match, Op::Match, Op::Mismatch, Op::Match, Op::Ins]
        );
        let xs = b"AAAAA";
        let ys = b"ATG";
        let (dist, ops) = crate::edit_dist_diff::edit_dist(xs, ys);
        assert_eq!(dist, 4);
        let ylen = ops.ops.iter().filter(|&&op| op != Op::Del).count();
        assert_eq!(ylen, 3);
        let xs = b"ACCCAC";
        let ys = b"ACC";
        crate::edit_dist_diff::edit_dist(xs, ys);
    }
    #[test]
    fn edit_dist_random() {
        use rand_xoshiro::Xoroshiro128Plus;
        let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(32909);
        for _ in 0..100 {
            let seq = kiley::gen_seq::generate_seq(&mut rng, 30);
            let prof = kiley::gen_seq::Profile::new(0.05, 0.05, 0.05);
            let seq2 = kiley::gen_seq::introduce_randomness(&seq, &mut rng, &prof);
            let (dist_1, _) = crate::edit_dist_usual::edit_dist(&seq, &seq2);
            assert!(dist_1 < 30);
            let (dist_2, _) = crate::edit_dist_diff::edit_dist(&seq, &seq2);
            assert_eq!(dist_1, dist_2)
        }
    }
}
