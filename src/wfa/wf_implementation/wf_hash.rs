use crate::wfa::wavefront::*;
use std::collections::HashMap;
use std::hash::Hash;

/// Wavefront implementation based on an **Hashmap**; efficient when the number of diagonal actually
/// computed is much less than all the hypotheticals diagonals.
pub struct WavefrontHash<T> {
    min_diagonal: T,
    max_diagonal: T,
    low_diagonal: T,
    high_diagonal: T,
    offsets: HashMap<T, T>,
    pred_diagonal: HashMap<T, T>,
}

impl<T> WavefrontHash<T>
where T: num::NumCast + Copy {
    pub fn new(min_diagonal: isize, max_diagonal: usize) -> Self {
        
        let new_min_diagonal : T;

        if min_diagonal < 0 {
            new_min_diagonal = from_usize::<T>((-min_diagonal) as usize);
        } 
        else {
            new_min_diagonal = from_usize::<T>(min_diagonal as usize)
        }
               
        WavefrontHash {
            min_diagonal: new_min_diagonal, 
            max_diagonal: from_usize::<T>(max_diagonal),
            low_diagonal: from_usize::<T>(0), 
            high_diagonal: from_usize::<T>(0),
            offsets: HashMap::new(),
            pred_diagonal: HashMap::new(),
        }
    }

    #[inline(always)]
    fn get_actual_diagonal(&self, diagonal: isize) -> Option<usize> {
        let min_diagonal = - as_isize(self.min_diagonal);

        if diagonal < min_diagonal || diagonal > as_isize(self.max_diagonal) {
            return None;
        }

        Some((diagonal - min_diagonal) as usize)
    }
}

impl<T> Wavefront for WavefrontHash<T> 
where T: num::NumCast + Copy + std::cmp::Eq + Hash + Copy {

    #[inline(always)]
    fn get_min_diagonal(&self) -> isize {
        - (as_isize(self.min_diagonal))
    }

    #[inline(always)]
    fn get_max_diagonal(&self) -> isize {
        as_isize(self.max_diagonal)
    }

    #[inline(always)]
    fn get_low_diagonal(&self) -> isize {
        - (as_isize(self.low_diagonal))
    }

    #[inline(always)]
    fn get_high_diagonal(&self) -> isize {
        as_isize(self.high_diagonal)
    }

    fn set_low_diagonal(&mut self, diagonal: isize) -> bool {
        if diagonal >= - as_isize(self.min_diagonal) && diagonal <= as_isize(self.high_diagonal) {
            self.low_diagonal = from_isize::<T>(-diagonal);
            true
        }
        else {
            false
        }
    }

    fn set_high_diagonal(&mut self, diagonal: isize) -> bool {
        if diagonal >= - as_isize(self.low_diagonal) && diagonal <= as_isize(self.max_diagonal) {
            self.high_diagonal = from_isize::<T>(diagonal);
            true
        }
        else {
            false
        }
    }

    fn get_diagonal_offset(&self, diagonal: isize) -> Option<usize> {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            let offset = self.offsets.get(&from_usize::<T>(current_diagonal));

            if offset == None {
                return None;
            }

            Some(as_usize(*(offset.unwrap())))
        }
        else {
            None
        }
        
    }

    fn set_diagonal_offset(&mut self, diagonal: isize, offset: usize) -> bool {

        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets.insert(
                from_usize::<T>(current_diagonal), 
                from_usize::<T>(offset)
            );
            true
        }
        else {
            false
        }    
    } 

    fn set_pred_diagonal(&mut self, diagonal: isize, predecessor_diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            if let Some(current_pred_diagonal) = self.get_actual_diagonal(predecessor_diagonal) {
                self.pred_diagonal.insert(
                    from_usize::<T>(current_diagonal), 
                    from_usize::<T>(current_pred_diagonal)
                );
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
            let pred_diagonal = self.pred_diagonal.get(&from_usize::<T>(current_diagonal));

            if pred_diagonal == None {
                return None;
            }

            Some(as_isize(*pred_diagonal.unwrap()) - as_isize(self.min_diagonal))
        }
        else {
            None
        }
    }

    fn exist(&self, diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets.contains_key(&from_usize::<T>(current_diagonal))
        }
        else {
            false
        }
    }

    fn remove_diagonal(&mut self, diagonal: isize) -> bool {
        if let Some(current_diagonal) = self.get_actual_diagonal(diagonal) {
            self.offsets.remove(&from_usize::<T>(current_diagonal));
            true
        }
        else {
            false
        }
    }
}