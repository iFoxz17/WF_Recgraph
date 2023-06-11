use crate::wfa::wavefront::*;
use std::collections::HashMap;

/// Wavefront implementation based on an **Hashmap**; efficient when the number of diagonal actually
/// computed is much less than all the hypotheticals diagonals.
pub struct WavefrontHash {
    min_diagonal: usize,
    max_diagonal: usize,
    offsets: HashMap<usize, usize>,
    pred_diagonal: HashMap<usize, usize>,
}

impl WavefrontHash {
    pub fn new(min_diagonal: isize, max_diagonal: usize) -> Self {
        
        let new_min_diagonal : usize;

        if min_diagonal < 0 {
            new_min_diagonal = (-min_diagonal) as usize;
        } 
        else {
            new_min_diagonal = min_diagonal as usize
        }
               
        WavefrontHash {
            min_diagonal: new_min_diagonal, 
            max_diagonal,
            offsets: HashMap::new(),
            pred_diagonal: HashMap::new(),
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

impl Wavefront for WavefrontHash {

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
            let offset = self.offsets.get(&current_diagonal);

            if offset == None {
                return None;
            }

            Some(*(offset.unwrap()))
        }
        else {
            None
        }
        
    }

    fn set_diagonal_offset(&mut self, diagonal: isize, offset: usize) -> bool {

        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets.insert(current_diagonal, offset);
            true
        }
        else {
            false
        }    
    } 

    fn set_pred_diagonal(&mut self, diagonal: isize, predecessor_diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            if let Some(current_pred_diagonal) = self.get_actual_diagonal(predecessor_diagonal) {
                self.pred_diagonal.insert(current_diagonal, current_pred_diagonal);
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
            let pred_diagonal = self.pred_diagonal.get(&current_diagonal);

            if pred_diagonal == None {
                return None;
            }

            Some(*(pred_diagonal.unwrap()) as isize - self.min_diagonal as isize)
        }
        else {
            None
        }
    }

    fn exist(&self, diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets.contains_key(&current_diagonal)
        }
        else {
            false
        }
    }
}