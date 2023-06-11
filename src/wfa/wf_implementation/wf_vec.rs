use crate::wfa::wavefront::*;
use bitvec::prelude::*;

/// Wavefront implementation based on a **Vec**; efficient when the number of diagonal actually
/// computed is similar to the maximum hypotheticals number of diagonals.
pub struct WavefrontVec {
    min_diagonal: usize,
    max_diagonal: usize,
    offsets: Vec<usize>,
    exists: BitVec,
    pred_diagonal: Vec<usize>,
}

impl WavefrontVec {
    pub fn new(min_diagonal: isize, max_diagonal: usize) -> Self {
        let min_diagonal: usize = -min_diagonal as usize;
    
        let mut offsets = Vec::with_capacity(max_diagonal + min_diagonal + 1);
        let mut exists = BitVec::with_capacity(max_diagonal + min_diagonal + 1);
        let mut pred_diagonal = Vec::with_capacity(max_diagonal + min_diagonal + 1);
    
        for _i in 0..offsets.capacity() {
            offsets.push(0);
            exists.push(false);
            pred_diagonal.push(0);
        } 
               
        WavefrontVec {
            min_diagonal, 
            max_diagonal,
            offsets,
            exists,
            pred_diagonal,
        }
    }

    #[inline(always)]
    fn get_actual_diagonal(&self, diagonal: isize) -> Option<usize> {
        let min_diagonal = -(self.min_diagonal as isize);

        if diagonal < min_diagonal || diagonal > self.max_diagonal as isize {
            return None;
        }

        Some((diagonal + (self.min_diagonal as isize)) as usize)
    }
}

impl Wavefront for WavefrontVec {

    #[inline(always)]
    fn get_min_diagonal(&self) -> isize {
        -(self.min_diagonal as isize)
    }

    #[inline(always)]
    fn get_max_diagonal(&self) -> isize {
        self.max_diagonal as isize
    }

    fn get_diagonal_offset(&self, diagonal: isize) -> Option<usize> {
       
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            if !self.exists[current_diagonal] {
                None
            }
            else {
                Some(self.offsets[current_diagonal])
            }
        }
        else {
            None
        }
    }

    fn set_diagonal_offset(&mut self, diagonal: isize, offset: usize) -> bool {

        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets[current_diagonal] = offset;
            self.exists.set(current_diagonal, true);
            true
        }
        else {
            false
        }
        
    } 
    
    fn set_pred_diagonal(&mut self, diagonal: isize, predecessor_diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            if let Some(current_pred_diagonal) = self.get_actual_diagonal(predecessor_diagonal) {
                self.pred_diagonal[current_diagonal] = current_pred_diagonal;
                true
            }
            else {
                false
            }
        }
        else {
            false
        }
    }

    fn get_pred_diagonal(&self, diagonal: isize) -> Option<isize> {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            Some(self.pred_diagonal[current_diagonal] as isize - (self.min_diagonal as isize))
        }
        else {
            None
        }
    }

    fn exist(&self, diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.exists[current_diagonal]
        }
        else {
            false
        }
    }
}