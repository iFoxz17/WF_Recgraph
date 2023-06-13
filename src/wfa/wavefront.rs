/// **Interface** of a **Wavefront**; to mantain the algorithm **time complexity**, 
/// all the operations should run in **O(1)**. 
pub trait Wavefront {

    /// Returns the hypothetical **minor diagonal** <code>(-(path length) + 1)</code>. 
    fn get_min_diagonal(&self) -> isize;

    /// Returns the hypothetical **major diagonal** <code>(sequence.len() - 1)</code>. 
    fn get_max_diagonal(&self) -> isize;

    /// Returns the **offset** of <code>diagonal</code> (if exists).
    fn get_diagonal_offset(&self, diagonal: isize) -> Option<usize>;

    /// Sets the **offset** of <code>diagonal</code> (if exists) to <code>offset</code>.
    fn set_diagonal_offset(&mut self, diagonal: isize, offset: usize) -> bool;

    /// Sets the **previous diagonal** of <code>diagonal</code> to <code>pred_diagonal</code>.
    fn set_pred_diagonal(&mut self, diagonal: isize, pred_diagonal: isize) -> bool;

    /// Returns the **previous diagonal** of <code>diagonal</code> (if exists).
    fn get_pred_diagonal(&self, diagonal: isize) -> Option<isize>;

    /// Returns <code>true</code> if <code>diagonal</code> exists; <code>false</code> otherwise.
    fn exist(&self, diagonal: isize) -> bool;
}

#[inline(always)]
pub fn as_usize<T>(value: T) -> usize 
where T: num::NumCast {
    num::NumCast::from::<T>(value).unwrap()
}

#[inline(always)]
pub fn from_usize<T>(value: usize) -> T 
where T: num::NumCast {
    num::NumCast::from::<usize>(value).unwrap()
}

#[inline(always)]
pub fn as_isize<T>(value: T) -> isize 
where T: num::NumCast {
    num::NumCast::from::<T>(value).unwrap()
}

#[allow(unused)]
#[inline(always)]
pub fn from_isize<T>(value: isize) -> T 
where T: num::NumCast {
    num::NumCast::from::<isize>(value).unwrap()
}