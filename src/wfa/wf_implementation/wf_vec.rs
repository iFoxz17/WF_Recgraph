use crate::wfa::wavefront::*;
use bitvec::prelude::*;

/// Wavefront implementation based on a **Vec**; efficient when the number of diagonal actually
/// computed is similar to all the hypotheticals diagonals.
pub struct WavefrontVec<T> {
    min_diagonal: T,
    max_diagonal: T,
    offsets: Vec<T>,
    exists: BitVec,
    pred_diagonal: Vec<T>,
}

impl<T> WavefrontVec<T>
where T: num::NumCast + Copy {

    pub fn new(min_diagonal: isize, max_diagonal: usize) -> Self {
        let new_min_diagonal: usize;
        if min_diagonal < 0 {
            new_min_diagonal = -min_diagonal as usize; 
        } 
        else {
            new_min_diagonal = min_diagonal as usize; 
        }
    
        let mut offsets = Vec::<T>::with_capacity(max_diagonal + new_min_diagonal + 1);
        let mut exists = BitVec::with_capacity(max_diagonal + new_min_diagonal + 1);
        let mut pred_diagonal = Vec::<T>::with_capacity(max_diagonal + new_min_diagonal + 1);
    
        for _i in 0..offsets.capacity() {
            offsets.push(from_usize::<T>(0));
            exists.push(false);
            pred_diagonal.push(from_usize::<T>(0));
        } 
               
        WavefrontVec {
            min_diagonal: from_usize::<T>(new_min_diagonal), 
            max_diagonal: from_usize::<T>(max_diagonal),
            offsets,
            exists,
            pred_diagonal,
        }
    }

    #[inline(always)]
    fn get_actual_diagonal(&self, diagonal: isize) -> Option<usize> {
        let min_diagonal = -(as_usize(self.min_diagonal) as isize);

        if diagonal < min_diagonal || diagonal > as_usize(self.max_diagonal) as isize {
            return None;
        }

        Some((diagonal + as_usize(self.min_diagonal) as isize) as usize)
    }
}

impl<T> Wavefront for WavefrontVec<T> 
where T: num::NumCast + Copy {

    #[inline(always)]
    fn get_min_diagonal(&self) -> isize {
        - (as_usize(self.min_diagonal) as isize)
    }

    #[inline(always)]
    fn get_max_diagonal(&self) -> isize {
        as_usize(self.max_diagonal) as isize
    }

    fn get_diagonal_offset(&self, diagonal: isize) -> Option<usize> {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            if !self.exists[current_diagonal] {
                None
            }
            else {
                Some(as_usize(self.offsets[current_diagonal]))
            }
        }
        else {
            None
        }
    }

    fn set_diagonal_offset(&mut self, diagonal: isize, offset: usize) -> bool {

        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets[current_diagonal] = from_usize::<T>(offset);
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
                self.pred_diagonal[current_diagonal] = from_usize::<T>(current_pred_diagonal);
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
            Some(as_usize(self.pred_diagonal[current_diagonal]) as isize - as_usize(self.min_diagonal) as isize)
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