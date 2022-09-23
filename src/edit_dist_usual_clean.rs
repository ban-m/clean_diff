use super::alignments::*;

#[derive(Debug, Clone, Copy)]
struct OpDist(u64);

fn to_u64(op: Op) -> u64 {
    match op {
        Op::Match => 0,
        Op::Del => 1,
        Op::Ins => 2,
        Op::Mismatch => 3,
    }
}

impl OpDist {
    fn new(dist: u32, num_gap: u16, on_op: Op, traceback: Op) -> Self {
        let code = ((dist as u64) << 32)
            | ((num_gap as u64) << 16)
            | (to_u64(on_op) << 8)
            | to_u64(traceback);
        Self(code)
    }
    fn init(init: u32) -> Self {
        Self(((init as u64) << 32) | 0b1111_0000)
    }
    fn mat_move(self, mat: bool) -> Self {
        let Self(mut bits) = self;
        bits += (!mat as u64) << 32;
        // Set the ON_GAP_FLAG as zero (match).
        bits = bits & !ON_GAP_FLAG;
        // Set the prev score as match
        let op = if mat { 0 } else { 3 };
        bits = (bits & !PREV_OF_FLAG) | op;
        Self(bits)
    }
    fn del_move(self) -> Self {
        let Self(mut bits) = self;
        bits += 1 << 32;
        // Check if this is consective deletion.
        if ((bits & ON_GAP_FLAG) >> 8) != 0b01 {
            bits += 1 << 16;
            bits = (bits & !ON_GAP_FLAG) | (1 << 8);
        }
        bits = (bits & !PREV_OF_FLAG) | 0b01;
        Self(bits)
    }
    fn ins_move(self) -> Self {
        let Self(mut bits) = self;
        bits += 1 << 32;
        if ((bits & ON_GAP_FLAG) >> 8) != 0b10 {
            bits += 1 << 16;
            bits = (bits & !ON_GAP_FLAG) | (0b10 << 8);
        }
        bits = (bits & !PREV_OF_FLAG) | 0b10;
        Self(bits)
    }
    fn min(&self, other: Self) -> Self {
        let min = (self.0).min(other.0);
        Self(min)
    }
    fn num_gaps(&self) -> u32 {
        ((self.0 << 32) >> 48) as u32
    }
    fn score(&self) -> u32 {
        (self.0 >> 32) as u32
    }
    fn prev_op(&self) -> Option<Op> {
        match self.0 & PREV_OF_FLAG {
            0 => Some(Op::Match),
            1 => Some(Op::Del),
            2 => Some(Op::Ins),
            3 => Some(Op::Mismatch),
            _ => None,
        }
    }
}

const ON_GAP_FLAG: u64 = 0b1111_1111 << 8;
const PREV_OF_FLAG: u64 = 0b1111_1111;

/// Usual edit distance alignments and its path.
pub fn edit_dist(xs: &[u8], ys: &[u8]) -> (u32, Alignment) {
    let max = (xs.len() + ys.len() + 3) as u32;
    let mut dp = vec![vec![OpDist::init(max); ys.len() + 1]; xs.len() + 1];
    for i in 1..xs.len() + 1 {
        dp[i][0] = OpDist::new(i as u32, 1, Op::Del, Op::Del);
    }
    for j in 1..ys.len() + 1 {
        dp[0][j] = OpDist::new(j as u32, 1, Op::Ins, Op::Ins);
    }
    dp[0][0] = OpDist::init(0);
    for (i, x) in xs.iter().enumerate().map(|(i, &x)| (i + 1, x)) {
        for (j, y) in ys.iter().enumerate().map(|(j, &y)| (j + 1, y)) {
            let mat_score = dp[i - 1][j - 1].mat_move(x == y);
            let del_score = dp[i - 1][j].del_move();
            let ins_score = dp[i][j - 1].ins_move();
            dp[i][j] = mat_score.min(del_score).min(ins_score);
        }
    }
    // Traceback.
    let (mut xpos, mut ypos) = (xs.len(), ys.len());
    let mut last = dp[xs.len()][ys.len()];
    let dist = last.score();
    let num_gaps = last.num_gaps();
    let mut ops = vec![];
    while let Some(op) = last.prev_op() {
        ops.push(op);
        match op {
            Op::Match | Op::Mismatch => {
                xpos -= 1;
                ypos -= 1;
            }
            Op::Ins => ypos -= 1,
            Op::Del => xpos -= 1,
        }
        last = dp[xpos][ypos];
    }
    ops.reverse();
    let aln = Alignment::new(ops);
    let (_, gaps) = aln.dist_and_num_of_gaps();
    assert_eq!(gaps, num_gaps, "{:?}", aln);
    (dist, aln)
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
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
            let (dist_2, _) = crate::edit_dist_usual_clean::edit_dist(&seq, &seq2);
            assert_eq!(dist_1, dist_2)
        }
    }
}
