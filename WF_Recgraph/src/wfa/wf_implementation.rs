/// <code>Vec</code> implementation.
pub mod wf_vec;
/// <code>Hashmap</code> implementation.
pub mod wf_hash;

#[allow(warnings)]

#[derive(Clone, Copy, Debug)]
pub enum WavefrontImpl {
    WavefrontHash,
    WavefrontVec,
    /// Mixed implementation: as long as the **percentage** of **diagonals**
    /// in the queue is less than the **float argument**, **WavefrontHash** is used; otherwise, it's chosen
    /// **WavefrontVec**. 
    WavefrontMixed(f64),
}