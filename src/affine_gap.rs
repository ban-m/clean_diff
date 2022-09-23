use crate::alignments::*;
pub fn align(
    xs: &[u8],
    ys: &[u8],
    mat_score: i64,
    mism: i64,
    gap_open: i64,
    gap_extend: i64,
) -> (i64, Alignment) {
    // (mat,del,ins)
    let min = mism.min(gap_open).min(gap_extend) * (xs.len() + ys.len() + 9) as i64;
    let mut dp = vec![vec![[(min, None); 3]; ys.len() + 1]; xs.len() + 1];
    for i in 1..xs.len() + 1 {
        dp[i][0][1] = (gap_open + (i - 1) as i64 * gap_extend, Some(1));
    }
    for j in 1..ys.len() + 1 {
        dp[0][j][2] = (gap_open + (j - 1) as i64 * gap_extend, Some(2));
    }
    dp[0][0][0] = (0, None);
    for (i, x) in xs.iter().enumerate().map(|(i, &x)| (i + 1, x)) {
        for (j, y) in ys.iter().enumerate().map(|(j, &y)| (j + 1, y)) {
            dp[i][j][0] = {
                let mat = if x == y { mat_score } else { mism };
                dp[i - 1][j - 1]
                    .iter()
                    .enumerate()
                    .map(|(state, (score, _))| (score + mat, Some(state)))
                    .max_by_key(|x| x.0)
                    .unwrap()
            };
            dp[i][j][1] = {
                dp[i - 1][j]
                    .iter()
                    .enumerate()
                    .map(|(state, (score, _))| match state {
                        0 | 2 => (score + gap_open, Some(state)),
                        1 => (score + gap_extend, Some(state)),
                        _ => panic!(),
                    })
                    .max_by_key(|x| x.0)
                    .unwrap()
            };
            dp[i][j][2] = {
                dp[i][j - 1]
                    .iter()
                    .enumerate()
                    .map(|(state, (score, _))| match state {
                        0 | 1 => (score + gap_open, Some(state)),
                        2 => (score + gap_extend, Some(state)),
                        _ => panic!(),
                    })
                    .max_by_key(|x| x.0)
                    .unwrap()
            };
        }
    }
    // Traceback
    assert!(dp[0][0].iter().all(|x| x.1.is_none()));
    let (mut xpos, mut ypos) = (xs.len(), ys.len());
    let (state, &(score, _)) = dp[xpos][ypos]
        .iter()
        .enumerate()
        .max_by_key(|x| x.1 .0)
        .unwrap();
    let mut ops = vec![];
    let mut state = Some(state);
    while let Some(st) = state {
        state = dp[xpos][ypos][st].1;
        if state == None {
            break;
        }
        match st {
            0 => {
                xpos -= 1;
                ypos -= 1;
                if xs[xpos] == ys[ypos] {
                    ops.push(Op::Match)
                } else {
                    ops.push(Op::Mismatch)
                }
            }
            1 => {
                xpos -= 1;
                ops.push(Op::Del);
            }
            2 => {
                ypos -= 1;
                ops.push(Op::Ins);
            }
            _ => unreachable!(),
        }
    }
    ops.reverse();
    let aln = Alignment::new(ops);
    (score, aln)
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn affine_gap_test() {
        let xs = b"AAACCC";
        let (dist, ops) = align(xs, xs, 0, -1, -1, -1);
        assert_eq!(dist, 0);
        assert_eq!(ops.ops, vec![Op::Match; xs.len()]);
        let ys = b"AAA";
        let (dist, aln) = align(xs, ys, 0, -1, -1, -1);
        assert_eq!(dist, -3);
        assert_eq!(aln.ops, vec![vec![Op::Match; 3], vec![Op::Del; 3]].concat());
        let ys = b"";
        let (dist, aln) = align(xs, ys, 0, -1, -1, -1);
        assert_eq!(dist, -6);
        assert_eq!(aln.ops, vec![Op::Del; 6]);
        let (dist, aln) = align(ys, xs, 0, -1, -1, -1);
        assert_eq!(dist, -6);
        assert_eq!(aln.ops, vec![Op::Ins; 6]);
        let xs = b"ACGT";
        let ys = b"ACCTG";
        let (dist, aln) = align(xs, ys, 0, -1, -1, -1);
        assert_eq!(dist, -2);
        assert_eq!(
            aln.ops,
            vec![Op::Match, Op::Match, Op::Mismatch, Op::Match, Op::Ins]
        );
        let xs = b"ACCCGCCCA";
        let ys = b"AGA";
        let (dist, aln) = align(xs, ys, 0, -1, -10, -1);
        assert_eq!(dist, -1 - 10 - 5);
        assert_eq!(
            aln.ops,
            vec![
                vec![Op::Match],
                vec![Op::Mismatch],
                vec![Op::Del; 6],
                vec![Op::Match]
            ]
            .concat()
        );
    }
}
