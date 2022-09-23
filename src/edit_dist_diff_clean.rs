use crate::alignments::*;
/// Usual edit distance alignments and its path.
pub fn edit_dist(xs: &[u8], ys: &[u8]) -> (u32, Alignment) {
    if xs == ys {
        return (0, Alignment::new(vec![Op::Match; xs.len()]));
    }
    // h -> d = (the f.r.p of the d-h diagonal with edit distance h, the num of gaps, w. on indel, the tracing).
    // 32bits, 16bits, 8bits, 8bits.
    let init_match = match_len(xs, 0, ys, 0);
    let mut dp = vec![vec![ReachPoint::new(init_match)]];
    for dist in 1..xs.len() + ys.len() + 2 {
        let new_d: Vec<_> = (0..2 * dist + 1)
            .map(|diag| {
                let from_above = dp[dist - 1].get(diag).map(|x| x.from_above());
                let from_mat = match diag < 1 {
                    true => None,
                    false => dp[dist - 1].get(diag - 1).map(|x| x.from_mat()),
                };
                let from_left = match diag < 2 {
                    true => None,
                    false => dp[dist - 1].get(diag - 2).map(|x| x.from_left()),
                };
                let max_reach = max_of_three(from_above, from_mat, from_left);
                let reached_pos = max_reach.position();
                let i = (reached_pos + dist) - diag;
                let snake_len = match_len(xs, i, ys, reached_pos);
                max_reach.add(snake_len)
            })
            .collect();
        let reached = new_d.iter().enumerate().any(|(diag, pos)| {
            let i = (pos.position() + dist) - diag;
            i == xs.len() && pos.position() == ys.len()
        });
        dp.push(new_d);
        if reached {
            break;
        }
        assert!(dist < xs.len() + ys.len() + 1);
    }
    let opt_dist = dp.len() - 1;
    let (mut diag, opt_cell) = dp
        .last()
        .unwrap()
        .iter()
        .enumerate()
        .find(|(diag, pos)| {
            let i = (pos.position() + opt_dist) - diag;
            i == xs.len() && pos.position() == ys.len()
        })
        .unwrap();
    let mut prev = *opt_cell;
    let mut ops = vec![];
    let mut dist = opt_dist;
    while let Some(trace) = prev.trace() {
        let old_ypos = prev.position();
        match trace {
            Op::Del => {}
            Op::Mismatch => diag -= 1,
            Op::Ins => diag -= 2,
            _ => panic!(),
        }
        dist -= 1;
        prev = dp[dist][diag];
        let new_ypos = prev.position();
        let len = match trace {
            Op::Del => old_ypos - new_ypos,
            Op::Mismatch | Op::Ins => old_ypos - new_ypos - 1,
            _ => panic!(),
        };
        ops.extend(std::iter::repeat(Op::Match).take(len));
        ops.push(trace);
    }
    assert_eq!(dist, 0);
    let ypos = prev.position();
    let xpos = (ypos + dist) - diag;
    assert_eq!(xpos, ypos);
    ops.extend(std::iter::repeat(Op::Match).take(xpos));
    ops.reverse();
    let aln = Alignment::new(ops);
    (opt_dist as u32, aln)
}

#[derive(Debug, Clone, Copy)]
struct ReachPoint(u64);

// TODO: Make the mid 16 bits to be encoded as flipped.
impl ReachPoint {
    fn trace(&self) -> Option<Op> {
        match self.0 & PREV_OF_FLAG {
            0 => Some(Op::Match),
            1 => Some(Op::Del),
            2 => Some(Op::Ins),
            3 => Some(Op::Mismatch),
            _ => None,
        }
    }
    fn new(point: usize) -> Self {
        Self(((point as u64) << 32) | 0b1111_0000)
    }
    fn add(&self, len: usize) -> Self {
        let &Self(mut bits) = self;
        bits += (len as u64) << 32;
        if 0 < len {
            bits &= !ON_GAP_FLAG;
        }
        Self(bits)
    }
    fn position(&self) -> usize {
        let &Self(bits) = self;
        (bits >> 32) as usize
    }
    fn from_above(&self) -> Self {
        let Self(mut bits) = self;
        // If this is not on the deletion line, the # of gap should be updated.
        if ((bits & ON_GAP_FLAG) >> 8) != DEL {
            bits += 1 << 16;
            bits = (bits & !ON_GAP_FLAG) | (DEL << 8);
        }
        bits = (bits & !PREV_OF_FLAG) | DEL;
        Self(bits)
    }
    fn from_mat(&self) -> Self {
        let Self(mut bits) = self;
        // It does not increase the # of gaps, but change the condition though...
        bits = bits & !ON_GAP_FLAG;
        bits = (bits & !PREV_OF_FLAG) | 0b11;
        bits += 1 << 32;
        Self(bits)
    }
    fn from_left(&self) -> Self {
        let Self(mut bits) = self;
        // If this is not on the insertion line, the # of gaps should be updated.
        if ((bits & ON_GAP_FLAG) >> 8) != INS {
            bits += 1 << 16;
            bits = (bits & !ON_GAP_FLAG) | (INS << 8);
        }
        bits += 1 << 32;
        bits = (bits & !PREV_OF_FLAG) | INS;
        Self(bits)
    }
    fn max(&self, other: Self) -> Self {
        let &Self(me) = self;
        let Self(other) = other;
        if (me ^ NUM_GAP) < (other ^ NUM_GAP) {
            Self(other)
        } else {
            Self(me)
        }
    }
}
const DEL: u64 = 0b01;
const INS: u64 = 0b10;
const NUM_GAP: u64 = ((1 << 32) - 1) - ((1 << 16) - 1);
const ON_GAP_FLAG: u64 = 0b1111_1111 << 8;
const PREV_OF_FLAG: u64 = 0b1111_1111;

fn max_of_three(x: Option<ReachPoint>, y: Option<ReachPoint>, z: Option<ReachPoint>) -> ReachPoint {
    match (x, y, z) {
        (None, None, Some(z)) => z,
        (None, Some(y), None) => y,
        (Some(x), None, None) => x,
        (None, Some(y), Some(z)) => z.max(y),
        (Some(x), None, Some(z)) => x.max(z),
        (Some(x), Some(y), None) => x.max(y),
        (Some(x), Some(y), Some(z)) => x.max(y).max(z),
        (None, None, None) => panic!(),
    }
}

fn match_len(xs: &[u8], x_start: usize, ys: &[u8], y_start: usize) -> usize {
    let xs = xs.iter().skip(x_start);
    let ys = ys.iter().skip(y_start);
    std::iter::zip(xs, ys).take_while(|(x, y)| x == y).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    #[test]
    fn edit_dist() {
        let xs = b"AAAAA";
        let ys = b"ATG";
        let (dist, ops) = crate::edit_dist_diff_clean::edit_dist(xs, ys);
        assert_eq!(dist, 4);
        let ylen = ops.ops.iter().filter(|&&op| op != Op::Del).count();
        assert_eq!(ylen, 3);
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
            let (dist_2, _) = crate::edit_dist_diff_clean::edit_dist(&seq, &seq2);
            assert_eq!(dist_1, dist_2)
        }
    }
}
