//! This module defines the Alignment.
//!

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub enum Op {
    Match,
    Mismatch,
    Ins,
    Del,
}

impl std::convert::TryFrom<char> for Op {
    type Error = ();
    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            '=' => Ok(Op::Match),
            'X' => Ok(Op::Mismatch),
            'I' => Ok(Op::Ins),
            'D' => Ok(Op::Del),
            _ => Err(()),
        }
    }
}

impl std::convert::From<Op> for char {
    fn from(value: Op) -> Self {
        match value {
            Op::Match => '=',
            Op::Mismatch => 'X',
            Op::Ins => 'I',
            Op::Del => 'D',
        }
    }
}

impl Op {}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub ops: Vec<Op>,
}

impl Alignment {
    pub fn new(ops: Vec<Op>) -> Self {
        Self { ops }
    }
    pub fn to_string(&self) -> String {
        self.ops.iter().map(|&x| -> char { x.into() }).collect()
    }
    pub fn from_str(aln: &str) -> Option<Self> {
        let mut ops: Vec<Op> = Vec::with_capacity(aln.len());
        for op in aln.chars() {
            ops.push(op.try_into().ok()?);
        }
        Some(Self { ops })
    }
    pub fn dist_and_num_of_gaps(&self) -> (u32, u32) {
        let dist = self.ops.iter().filter(|&&op| op != Op::Match).count() as u32;
        let mut num_of_gap = 0;
        let mut prev = None;
        for &op in self.ops.iter() {
            // Del gap open.
            if op == Op::Del && prev != Some(Op::Del) {
                num_of_gap += 1;
            }
            if op == Op::Ins && prev != Some(Op::Ins) {
                num_of_gap += 1;
            }
            prev = Some(op);
        }
        (dist, num_of_gap)
    }
    /// xs is the reference, ys is the query.
    pub fn recover(&self, xs: &[u8], ys: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (mut i, mut j) = (0, 0);
        let (mut xr, mut yr, mut aln) = (vec![], vec![], vec![]);
        for &op in self.ops.iter() {
            match op {
                Op::Mismatch | Op::Match => {
                    xr.push(xs[i]);
                    yr.push(ys[j]);
                    if xs[i] == ys[j] {
                        aln.push(b'|');
                    } else {
                        aln.push(b'X');
                    }
                    i += 1;
                    j += 1;
                }
                Op::Del => {
                    xr.push(xs[i]);
                    aln.push(b' ');
                    yr.push(b' ');
                    i += 1;
                }
                Op::Ins => {
                    xr.push(b' ');
                    aln.push(b' ');
                    yr.push(ys[j]);
                    j += 1;
                }
            }
        }
        (xr, aln, yr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn recover() {
        let ops = "=XID";
        let aln = Alignment::from_str(ops).unwrap();
        assert_eq!(vec![Op::Match, Op::Mismatch, Op::Ins, Op::Del], aln.ops);
        let ops = "XX=====IDIDIDIDID====XXXXIDID";
        let aln = Alignment::from_str(ops).unwrap();
        let rec = Alignment::to_string(&aln);
        assert_eq!(ops, rec);
    }
    #[test]
    fn eval() {
        let ops = Alignment::from_str("=========").unwrap();
        let (dist, gaps) = ops.dist_and_num_of_gaps();
        assert_eq!(dist, 0);
        assert_eq!(gaps, 0);
        let ops = Alignment::from_str("XXXXX==").unwrap();
        let dg = ops.dist_and_num_of_gaps();
        assert_eq!(dg, (5, 0));
        let ops = Alignment::from_str("DDDDIIII").unwrap();
        let dg = ops.dist_and_num_of_gaps();
        assert_eq!(dg, (8, 2));
        let ops = Alignment::from_str("IIDD=D=DI=").unwrap();
        let dg = ops.dist_and_num_of_gaps();
        assert_eq!(dg, (7, 5));
    }
}
