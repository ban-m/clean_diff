use super::alignments::*;
/// Usual edit distance alignments and its path.
pub fn edit_dist(xs: &[u8], ys: &[u8]) -> (u32, Alignment) {
    let mut dp = vec![vec![(0, None); ys.len() + 1]; xs.len() + 1];
    for i in 1..xs.len() + 1 {
        dp[i][0] = (i as u32, Some(Op::Del));
    }
    for j in 1..ys.len() + 1 {
        dp[0][j] = (j as u32, Some(Op::Ins));
    }
    for (i, x) in xs.iter().enumerate().map(|(i, &x)| (i + 1, x)) {
        for (j, y) in ys.iter().enumerate().map(|(j, &y)| (j + 1, y)) {
            let mat_score = dp[i - 1][j - 1].0 + (x != y) as u32;
            let ins_score = dp[i][j - 1].0 + 1;
            let del_score = dp[i - 1][j].0 + 1;
            let min = mat_score.min(ins_score).min(del_score);
            dp[i][j] = if mat_score == min {
                if x == y {
                    (min, Some(Op::Match))
                } else {
                    (min, Some(Op::Mismatch))
                }
            } else if ins_score == min {
                (min, Some(Op::Ins))
            } else if del_score == min {
                (min, Some(Op::Del))
            } else {
                unreachable!()
            };
        }
    }
    assert_eq!(dp[0][0], (0, None));
    let (dist, mut operation) = dp[xs.len()][ys.len()];
    let mut aln = vec![];
    let (mut xpos, mut ypos) = (xs.len(), ys.len());
    while let Some(op) = operation {
        aln.push(op);
        match op {
            Op::Match | Op::Mismatch => {
                xpos -= 1;
                ypos -= 1;
            }
            Op::Ins => ypos -= 1,
            Op::Del => xpos -= 1,
        }
        operation = dp[xpos][ypos].1;
    }
    aln.reverse();
    let aln = Alignment::new(aln);
    (dist, aln)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn edit_dist_calc() {
        let xs = b"AAACCC";
        let (dist, ops) = edit_dist(xs, xs);
        assert_eq!(dist, 0);
        assert_eq!(ops.ops, vec![Op::Match; xs.len()]);
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
    }
}
