extern crate queues;

use crate::wfa::wavefront::*;
use crate::wfa::wf_implementation::*;
use crate::wfa::wf_implementation::wf_vec::*;
use crate::wfa::wf_implementation::wf_hash::*;

use std::thread;
use std::sync::{Arc, Mutex};
use std::sync::mpsc;

use std::collections::HashSet;
use std::hash::Hash;
use queues::*;

/// **Multiple string rapresentation** of the **graph**; <code>lnz</code> contains the linearization of 
/// the graph, while every path is mapped in <code>lnz</code> by <code>paths_mapping</code>. 
/// # Example
/// <code>lnz =    ['A', 'T', 'G', 'T', 'G', 'T', 'G', 'C']</code>  
/// <code>path_1 = ['A', 'T', 'G',                'G', 'C']</code>  
/// <code>path_2 = ['A', 'T',           'G', 'T', 'G', 'C']</code>  
/// <code>path_3 = ['A', 'T',      'T',           'G', 'C']</code>  
/// <code>paths_mapping = [[0, 1, 2, 6, 7], [0, 1, 4, 5, 6, 7], [0, 1, 3, 6, 7]]</code>  
/// <code>paths_number = 3</code>
pub struct PathStrings {
    lnz: Vec<char>,
    paths_mapping: Vec<Vec<usize>>,
    paths_number: usize,
}

impl PathStrings {
    pub fn new(lnz: Vec<char>, paths_mapping: Vec<Vec<usize>>, paths_number: usize) -> PathStrings {
        PathStrings{
            lnz,
            paths_mapping,
            paths_number,
        }
    } 

    #[inline(always)]
    fn get_path_node_abs_mapping(&self, path: usize, v: usize) -> usize { 
        self.paths_mapping[path][v]
    }

    #[inline(always)]
    fn get_path_node_label(&self, path: usize, v: usize) -> char { 
        self.lnz[self.get_path_node_abs_mapping(path, v)]
    }

    #[inline(always)]
    pub fn get_min_diagonal(&self, path: usize) -> isize {
        -(self.paths_mapping[path].len() as isize) + 1
    } 

    #[inline(always)]
    pub fn get_path_length(&self, path: usize) -> usize {
        self.paths_mapping[path].len()
    } 
}

/// Struct which stores a **set of wavefront** (all the **wavefronts** of scores within the range 
/// <code>[0, score_wavefront.len()]</code>).
struct WavefrontSet {
    score_wavefront: Vec<Box<dyn Wavefront>>,
}

impl WavefrontSet {
    fn new<T>(
        min_diagonal: isize,
        max_diagonal: usize, 
        wavefront_impl: WavefrontImpl,
        diagonals_number: usize,
        max_penalty: usize
    ) -> WavefrontSet 
    where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {

        let mut wavefront_set = WavefrontSet {
            score_wavefront: Vec::new(),
        };

        wavefront_set.add_wavefront::<T>(
            min_diagonal, max_diagonal, wavefront_impl, diagonals_number, max_penalty
        );
        wavefront_set
    }

    fn add_wavefront<T>(
        &mut self, 
        min_diagonal: isize, 
        max_diagonal: usize,
        wavefront_impl: WavefrontImpl,
        diagonals_number: usize,
        max_penalty: usize
    ) 
    where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {
        self.score_wavefront.push(
            match wavefront_impl {
                WavefrontImpl::WavefrontVec => Box::new(WavefrontVec::<T>::new(min_diagonal, max_diagonal)),
                WavefrontImpl::WavefrontHash => Box::new(WavefrontHash::<T>::new(min_diagonal, max_diagonal)),
                WavefrontImpl::WavefrontMixed(threshold) => {
                    if diagonals_number as f64 / 
                    (((max_diagonal + (-min_diagonal) as usize + 1) * max_penalty) as f64) 
                    < threshold {
                        Box::new(WavefrontVec::<T>::new(min_diagonal, max_diagonal))
                    }
                    else {
                        Box::new(WavefrontHash::<T>::new(min_diagonal, max_diagonal))
                    }
                },
            }
        )
    }

    fn get_wavefront(&mut self, score: usize) -> Option<&mut Box<dyn Wavefront>> {

        if score < self.score_wavefront.len() {
            Some(&mut self.score_wavefront[score])
        }
        else {
            None
        }
    }

    /// Returns the max wavefront score stored.
    #[inline(always)]
    fn get_max_score(&self) -> usize {
        if self.score_wavefront.len() > 0 {
            self.score_wavefront.len() - 1
        }
        else {
            0
        }
    }
}

/// Struct which stores all the **wavefronts** of a certain **path**.
/// # Fields
/// - <code>path</code>: **path** chosen;
/// - <code>wavefronts</code>: all the **wavefronts** of the **path**; 
/// - <code>diagonals_queue</code>: <code>Queue</code> that maintains the diagonals to extend for every score.
pub struct PathWavefronts {
    path: usize,
    wavefronts: WavefrontSet,
    diagonals_queue: Queue<(usize, isize)>,
}

impl PathWavefronts {
    #[inline(always)]
    fn new<T>(
        path: usize, 
        min_diagonal: isize, 
        max_diagonal: usize, 
        wavefront_impl: WavefrontImpl,
        diagonals_number: usize,
        max_penalty: usize
    ) -> PathWavefronts 
    where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {
        PathWavefronts {
            path,
            wavefronts: WavefrontSet::new::<T>(
                min_diagonal, max_diagonal, wavefront_impl, diagonals_number, max_penalty
            ),
            diagonals_queue: Queue::new(),
        }
    }
}

#[allow(unused)]
#[derive(Clone, Copy)]
pub enum AlignmentMod {
    Global,
    Semiglobal,
    StartFreeGlobal,
    EndFreeGlobal,
}

#[derive(Clone, Copy, Debug)]
pub enum AlignmentOperation {
    Match(char),
    Mismatch(char, char),
    Insertion(char),
    Deletion(char),
}

#[derive(Debug)]
/// Struct that stores the result of an **alignment**;
/// # Fields
/// - <code>path</code>: **path** of the **optimal alignment**;
/// - <code>path_nodes</code>: **nodes** of the path actually **aligned**; 
/// - <code>path_length</code>: **path length**; 
/// - <code>path_start</code>: **first node** actually **aligned** of the **path**; 
/// - <code>path_end</code>: **last node** actually **aligned** of the **path**; 
/// - <code>score</code>: **score** of the alignment; 
/// - <code>operations</code>: <code>Vec</code> of the **alignemnt operations performed**.
/// - <code>diagonals_pruned</code>: number of **diagonals pruned** by the algorithm.
pub struct Alignment {
    pub path: usize,
    pub path_nodes: Vec<usize>,
    pub path_length: usize,
    pub path_start: usize,
    pub path_end: usize,
    pub score: usize,
    pub operations: Vec<AlignmentOperation>,
    pub diagonals_pruned: usize
} 

impl Alignment {
    #[inline(always)]
    pub fn new(path: usize, path_length: usize, path_start: usize, 
        path_end: usize, score: usize, diagonals_pruned: usize) -> Alignment {
        Alignment {
            path,
            path_nodes: Vec::with_capacity(path_length),
            path_length,
            path_start,
            path_end,
            score,
            operations: Vec::new(),
            diagonals_pruned,
        }
    }

    fn get_graph_operations_sequence_strings(&self) -> (String, String, String) {
        let mut graph_string = String::new();
        let mut op_string = String::new();
        let mut sequence_string = String::new();

        for op in self.operations.iter() {
            match op {
                AlignmentOperation::Match(match_car) => {
                    graph_string.push(*match_car);
                    op_string.push('M');
                    sequence_string.push(*match_car);
                }
                AlignmentOperation::Insertion(inserted_car) => {
                    graph_string.push('-');
                    op_string.push('I');
                    sequence_string.push(*inserted_car);
                    
                }
                AlignmentOperation::Mismatch(graph_car, sequence_car) => {
                    graph_string.push(*graph_car);
                    op_string.push('X');
                    sequence_string.push(*sequence_car);
                }
                AlignmentOperation::Deletion(deleted_car) => {
                    graph_string.push(*deleted_car);
                    op_string.push('D');
                    sequence_string.push('-');
                }
            }
        }

        (graph_string, op_string, sequence_string)
    }

    #[allow(unused)]
    pub fn get_path_string(&self) -> String {
        let mut path_string = String::new();
        for op in self.operations.iter() {
            match op {
                AlignmentOperation::Match(match_car) => {
                    path_string.push(*match_car);
                }
                AlignmentOperation::Insertion(inserted_car) => {
                    path_string.push(*inserted_car);
                }
                AlignmentOperation::Mismatch(graph_car, _sequence_car) => {
                    path_string.push(*graph_car);
                }
                AlignmentOperation::Deletion(_deleted_car) => {},
            }
        }

        path_string
    }

    fn change_operation(&self, new_op: char, last_op: &mut char, op_count: &mut usize, cigar: &mut String) {
        if *last_op == new_op  {
            *op_count += 1;
        }
        else {
            if *op_count > 0 {
                (*cigar).push_str(&(*op_count).to_string());
                (*cigar).push(*last_op);
            }
            *last_op = new_op;
            *op_count = 1;
        }
    }

    #[allow(unused)]
    pub fn get_cigar(&self) -> String {
        let mut cigar = String::new();
        let mut last_op = '\0';
        let mut op_count = 0;

        for op in self.operations.iter() {
            match op {
                AlignmentOperation::Match(_match_car) =>
                    self.change_operation('M', &mut last_op, &mut op_count, &mut cigar),
                AlignmentOperation::Insertion(_inserted_car) =>
                    self.change_operation('I', &mut last_op, &mut op_count, &mut cigar),
                AlignmentOperation::Mismatch(_graph_car, _sequence_car) =>
                    self.change_operation('X', &mut last_op, &mut op_count, &mut cigar),
                AlignmentOperation::Deletion(_deleted_car) => 
                    self.change_operation('D', &mut last_op, &mut op_count, &mut cigar),
            }
        }
        cigar.push_str(&op_count.to_string());
        cigar.push(last_op);
    
        cigar
    }
}

impl std::fmt::Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {

        let (graph_string, op_string, sequence_string) = self.get_graph_operations_sequence_strings();

        write!(f, "Path:{}\nPath length: {}\nPath start: {}\nPath end: {}\nAlignment penalty: {}\n{}\n{}\n{}", 
        self.path, self.path_length, self.path_start, self.path_end, 
        self.score, graph_string, op_string, sequence_string)
    }
}

#[inline(always)]
fn max(x: usize, y: usize, z: usize) -> usize {
    if x >= y {
        if x >= z {
            x
        }
        else {
            z
        }
    }
    else {
        if y >= z {
            y
        }
        else {
            z
        }
    }
}

fn set_base_case(
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    modality: AlignmentMod
) {
    match modality {
        AlignmentMod::Global | AlignmentMod::EndFreeGlobal => {
            path_wavefronts.wavefronts
            .get_wavefront(0)
            .unwrap()
            .set_diagonal_offset(0, 0);

            path_wavefronts.diagonals_queue.add((0, 0)).unwrap();
        },
        
        AlignmentMod::Semiglobal | AlignmentMod::StartFreeGlobal=> {
            for j in 0..graph.get_path_length(path_wavefronts.path) {
                path_wavefronts.wavefronts
                .get_wavefront(0)
                .unwrap()
                .set_diagonal_offset(-(j as isize), 0);
                path_wavefronts.diagonals_queue.add((0, -(j as isize))).unwrap(); 
            } 
        },
    }
}

fn check_termination(
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    score: usize,
    modality: AlignmentMod,
) -> Option<isize> {

    let mut final_diagonal = None;

    match modality {
        AlignmentMod::Global | AlignmentMod::StartFreeGlobal => {
            if let Some(offset) = path_wavefronts.wavefronts
            .get_wavefront(score)
            .unwrap().
            get_diagonal_offset(sequence.len() as isize - 
            graph.get_path_length(path_wavefronts.path) as isize) {

                if offset == sequence.len() - 1 {
                    final_diagonal = Some(sequence.len() as isize - 
                    graph.get_path_length(path_wavefronts.path) as isize);
                }
            }
        },

        AlignmentMod::Semiglobal | AlignmentMod::EndFreeGlobal => {
            for j in 0..graph.get_path_length(path_wavefronts.path) {
                if let Some(offset) = path_wavefronts.wavefronts
                .get_wavefront(score)
                .unwrap()
                .get_diagonal_offset(sequence.len() as isize - 1 - j as isize) {

                    if offset == sequence.len() - 1 {
                        final_diagonal = Some(sequence.len() as isize - 1 - j as isize);
                        break;
                    }
                }
            }
        },
    }

    final_diagonal
}

#[inline]
fn get_diagonal_match_count(
    sequence: &[char], 
    graph: &PathStrings,
    path: usize,
    k: isize,
    i: usize
) -> usize {
        let v = (i as isize - k) as usize;
        let mut increase = 0;

        while i + increase < sequence.len() - 1 && 
        v + increase < graph.paths_mapping[path].len() - 1 && 
        graph.get_path_node_label(path, v + increase + 1) == sequence[i + increase + 1] {
            increase += 1; 
        }

        increase
}

fn thread_run_condition(
    sequence: &[char], 
    graph: &PathStrings,
    path: usize,
    k: isize,
    i: usize
) -> bool {
        let v = (i as isize - k) as usize;

        i < sequence.len() - 1 && 
        v < graph.paths_mapping[path].len() - 1 && 
        graph.get_path_node_label(path, v + 1) == sequence[i +  1]
}

/// Every diagonal that respect the <code>thread_run_condition</code> is extended in parallel.
fn extend_with_threads(
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    score: usize
) {
    let queue = &mut path_wavefronts.diagonals_queue;
    let mut mem_queue = Queue::new();
    let wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

    thread::scope(|s| {

        let mut diagonals_to_extend = 0;  
        let (tx, rx) = mpsc::channel();
        let path = path_wavefronts.path;

        while (*queue).size() > 0 {    
            let (penalty, k) = queue.remove().unwrap();
            
            if score == penalty {
                let i = wavefront.get_diagonal_offset(k).unwrap();

                if thread_run_condition(sequence, graph, path, k, i) {
                    diagonals_to_extend += 1;
                    let tx = tx.clone();    
                    
                    s.spawn(move || {
                        let increase = get_diagonal_match_count(
                            sequence,
                            graph,
                            path,
                            k,
                            i
                        );
                        tx.send((k, increase)).unwrap();
                        });
                }
                else {
                    let increase = get_diagonal_match_count(sequence, graph, path_wavefronts.path, k, i);
                    wavefront.set_diagonal_offset(k, i + increase); 
                    mem_queue.add((penalty, k)).unwrap();
                }
            }
            else {
                mem_queue.add((penalty, k)).unwrap();
            }
        }

        while diagonals_to_extend > 0 {
            match rx.recv() {
                Ok((k, increase)) => {
                    wavefront.set_diagonal_offset(
                        k, 
                        wavefront.get_diagonal_offset(k).unwrap() + increase
                    );
                    mem_queue.add((score, k)).unwrap();
                },
                Err(_) => (),
            }
            diagonals_to_extend -= 1;
        }

        *queue = mem_queue;
    });
}

fn extend_no_threads(
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    score: usize,
) {

    let queue = &mut path_wavefronts.diagonals_queue;
    let mut mem_queue = Queue::new();
    let wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

    while (*queue).size() > 0 {    
        let (penalty, k) = queue.remove().unwrap();

        if score == penalty {
            let i = wavefront.get_diagonal_offset(k).unwrap();
            let increase = get_diagonal_match_count(sequence, graph, path_wavefronts.path, k, i);
            wavefront.set_diagonal_offset(k, i + increase);
        }

        mem_queue.add((penalty, k)).unwrap();
    }

    *queue = mem_queue;
}

#[inline(always)]
fn extend(
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    score: usize,
    parallelize_match: bool
) {

    if parallelize_match {
        extend_with_threads(
            sequence,
            graph,
            path_wavefronts,
            score
        )   
    }
    else {
        extend_no_threads(
            sequence,
            graph,
            path_wavefronts,
            score
        ) 
    }        
}

fn expand<T>(
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts,
    mem_set: &mut HashSet<(usize, isize)>,
    score: usize,
    m: usize,
    ins: usize,
    del: usize,
    wavefront_impl: WavefrontImpl
) 
where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {

    let max_penalty = max(m, ins, del);

    let queue = &mut path_wavefronts.diagonals_queue;
    let mut mem_queue: Queue<(usize, isize)> = Queue::new();

    let mut next_wavefront;

    for _add in 1..=(score + max_penalty - path_wavefronts.wavefronts.get_max_score()) {
        path_wavefronts.wavefronts.add_wavefront::<T>(
            graph.get_min_diagonal(path_wavefronts.path),
            sequence.len() - 1,
            wavefront_impl,
            (*queue).size(),
            max_penalty
        );
    }

    while (*queue).size() > 0 {                
        let (penalty, k) = queue.remove().unwrap();

        if penalty == score {
            let wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

            let i = (*wavefront).get_diagonal_offset(k).unwrap();
            let v = (i as isize - k) as usize;

            if i < sequence.len() - 1 {
                next_wavefront = path_wavefronts.wavefronts.get_wavefront(score + ins).unwrap(); 

                if next_wavefront.get_diagonal_offset(k + 1).unwrap_or(0) < i + 1 {
                    next_wavefront.set_diagonal_offset(k + 1, i + 1);
                    next_wavefront.set_pred_diagonal(k + 1, k);

                    if !mem_set.contains(&(score + ins, k + 1)) {
                        mem_queue.add((score + ins, k + 1)).unwrap();
                        mem_set.insert((score + ins, k + 1));
                    }
                }
            }

            if i < sequence.len() - 1 && v < graph.get_path_length(path_wavefronts.path) - 1 {
                next_wavefront = path_wavefronts.wavefronts.get_wavefront(score + m).unwrap();

                if next_wavefront.get_diagonal_offset(k).unwrap_or(0) < i + 1 {
                    next_wavefront.set_diagonal_offset(k, i + 1);
                    next_wavefront.set_pred_diagonal(k, k);

                    if !mem_set.contains(&(score + m, k)) {
                        mem_queue.add((score + m, k)).unwrap(); 
                        mem_set.insert((score + m, k));
                    } 
                }
            }

            if v < graph.get_path_length(path_wavefronts.path) - 1 {
                next_wavefront = path_wavefronts.wavefronts.get_wavefront(score + del).unwrap();
                let other_offset = next_wavefront.get_diagonal_offset(k - 1).unwrap_or(0);

                if other_offset < i || (other_offset == 0 && i == 0) {
                    next_wavefront.set_diagonal_offset(k - 1, i);
                    next_wavefront.set_pred_diagonal(k - 1, k);

                    if !mem_set.contains(&(score + del, k - 1)) {
                        mem_queue.add((score + del, k - 1)).unwrap();
                        mem_set.insert((score + del, k - 1));
                    }
                }
            }
        }
        else if !mem_set.contains(&(penalty, k)) {
            mem_queue.add((penalty, k)).unwrap();
            mem_set.insert((penalty, k));
        }
    }
    *queue = mem_queue;
    mem_set.clear();    
}

fn traceback( 
    sequence: &[char], 
    graph: &PathStrings,
    path_wavefronts: &mut PathWavefronts, 
    final_diagonal: isize,
    alignment_penalty: usize,
    m: usize,
    ins: usize,
    del: usize,
    diagonals_pruned: usize
) -> Alignment {

    let mut i = sequence.len() - 1;
    let final_v = (i as isize - final_diagonal) as usize;
    let mut v = final_v;
    let mut k = i as isize - v as isize;

    let mut alignment = Alignment::new(
        path_wavefronts.path,
        graph.paths_mapping[path_wavefronts.path].len(),
        0,
        graph.get_path_node_abs_mapping(path_wavefronts.path, v),
        alignment_penalty,
        diagonals_pruned
    );

    let mut operations = Vec::with_capacity(graph.paths_mapping[path_wavefronts.path].len());
    let mut path_nodes = Vec::with_capacity(graph.paths_mapping[path_wavefronts.path].len());

    let mut score = alignment_penalty;

    let mut wavefront;

    while score > 0 {
        
        while v > 0 && i > 0 && graph.get_path_node_label(path_wavefronts.path, v) == sequence[i] {
            operations.push(AlignmentOperation::Match(sequence[i]));
            path_nodes.push(graph.get_path_node_abs_mapping(path_wavefronts.path, v));

            i -= 1;
            v -= 1;
        }

        wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

        if v > 0 || i > 0 {
            if let Some(prev_k) = wavefront.get_pred_diagonal(k) {
                path_nodes.push(graph.get_path_node_abs_mapping(path_wavefronts.path, v));

                if prev_k == k - 1 {
                    operations.push(AlignmentOperation::Insertion(sequence[i]));

                    i -= 1;
                    score -= ins;
                }
                else if prev_k == k {    
                    operations.push(AlignmentOperation::Mismatch(
                        graph.get_path_node_label(path_wavefronts.path, v), sequence[i]));

                    v -= 1;
                    i -= 1;
                    score -= m;
                }
                else if prev_k == k + 1 {
                    operations.push(AlignmentOperation::Deletion(
                        graph.get_path_node_label(path_wavefronts.path, v)));
                    
                    v -= 1;
                    score -= del;
                }

                k = prev_k;
            }
        }           
    }

    while v > 0 && i > 0 && graph.get_path_node_label(path_wavefronts.path, v) == sequence[i] {
        operations.push(AlignmentOperation::Match(sequence[i]));
        path_nodes.push(graph.get_path_node_abs_mapping(path_wavefronts.path, v));

        v -= 1; 
        i -= 1;
    }

    for op in operations.iter().rev() {
        alignment.operations.push(*op);
    }

    for v in path_nodes.iter().rev() {
        alignment.path_nodes.push(*v);
    }

    alignment.path_length = final_v - v + 1;
    alignment.path_start = graph.get_path_node_abs_mapping(path_wavefronts.path, v);

    alignment
}

/// Sperimental for pruning heuristics. Peform the pruning in semiglobal mode.
fn prune_semiglobal(
    path_wavefronts: &mut PathWavefronts,
    score: usize,
    max_delta_allowed: usize
) -> usize {
    let mut max_offset = 0;
    let mut diagonals_pruned = 0;

    let queue = &mut path_wavefronts.diagonals_queue;
    let mut mem_queue: Queue<(usize, isize)> = Queue::new();
    let mut tmp_queue: Queue<isize> = Queue::new();

    let wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

    while (*queue).size() > 0 {
        let (penalty, k) = queue.remove().unwrap();
        
        if penalty == score {
            let i = wavefront.get_diagonal_offset(k).unwrap();
 
            if i > max_offset {
                max_offset = i;
            }
            
            tmp_queue.add(k).unwrap();
        }
        else {
            mem_queue.add((penalty, k)).unwrap();
        }
    }

    while tmp_queue.size() > 0 {
        let k = tmp_queue.remove().unwrap();
        let i = wavefront.get_diagonal_offset(k).unwrap();

        if max_offset - i <= max_delta_allowed {
            mem_queue.add((score, k)).unwrap();
        }
        else {
            diagonals_pruned += 1;
        } 
    }

    *queue = mem_queue;
    diagonals_pruned
}

/// Sperimental for pruning heuristics (not working). Peform the pruning in global mode.
fn prune_global( 
    path_wavefronts: &mut PathWavefronts,
    score: usize,
    max_delta_allowed: usize
) -> usize {
    let mut max_offset = 0;
    let mut diagonals_pruned = 0;

    let queue = &mut path_wavefronts.diagonals_queue;
    let mut mem_queue: Queue<(usize, isize)> = Queue::new();
    let mut tmp_queue: Queue<isize> = Queue::new();

    let wavefront = path_wavefronts.wavefronts.get_wavefront(score).unwrap();

    while (*queue).size() > 0 {
        let (penalty, k) = queue.remove().unwrap();
        
        if penalty == score {

            let i = wavefront.get_diagonal_offset(k).unwrap();
            let v = (i as isize - k) as usize;
 
            if i + v > max_offset {
                max_offset = i + v;
            }
            tmp_queue.add(k).unwrap();
        }
        else {
            mem_queue.add((penalty, k)).unwrap();
        }
    }

    while tmp_queue.size() > 0 {
        let k = tmp_queue.remove().unwrap();
        let i = wavefront.get_diagonal_offset(k).unwrap();
        let v = (i as isize - k) as usize;

        if max_offset - (i + v) <= max_delta_allowed {
            mem_queue.add((score, k)).unwrap();
        }
        else {
            //mem_queue.add((score, k)).unwrap();
            //wavefront.remove_diagonal(k);
            diagonals_pruned += 1;
        } 
    }

    *queue = mem_queue;
    diagonals_pruned
}

/// Sperimental for pruning heuristics.
#[inline]
fn prune( 
    path_wavefronts: &mut PathWavefronts,
    score: usize,
    align_mode: AlignmentMod,
    max_delta_allowed: usize
) -> usize {
    match align_mode {
        AlignmentMod::Semiglobal | AlignmentMod::EndFreeGlobal => 
            prune_semiglobal(path_wavefronts, score, max_delta_allowed),
        AlignmentMod::Global | AlignmentMod::StartFreeGlobal => 
            prune_global(path_wavefronts, score, max_delta_allowed),
    }
}

/// Perform the alignment between a **single path** of the **graph** and the **sequence**.
fn wf_align_to_path<T>(
    sequence: &[char],
    graph: &PathStrings,
    path: usize,
    m: usize,
    ins: usize,
    del: usize,
    wavefront_impl: WavefrontImpl,
    modality: AlignmentMod,
    parallelize_match: bool,
    max_delta_allowed: usize,
    prune_condition: fn(&PathStrings, &[char], usize, usize, usize, usize) -> bool,
    find_all_alignments: bool,
    maybe_max_score_mutex: Arc<Mutex<Option<usize>>>,
) -> Alignment 
where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {

    let mut path_wavefronts = PathWavefronts::new::<T>(
        path, 
        graph.get_min_diagonal(path), 
        sequence.len() - 1, 
        wavefront_impl, 
        0, 
        max(m, ins, del)
    );

    let max_penalty = max(m, ins, del);

    let mut mem_set: HashSet<(usize, isize)> = HashSet::new();

    let final_diagonal : isize;
    let mut diagonals_pruned = 0;

    set_base_case(graph, &mut path_wavefronts, modality);

    let mut d: usize = 0;

    loop {

        extend(sequence, graph, &mut path_wavefronts, d, parallelize_match);
        
        if let Ok(ref mut mutex) = maybe_max_score_mutex.try_lock() {
            if let Some(max_score) = **mutex {
                if d > max_score {
                    return Alignment::new(path, 0, 0, 0, d, diagonals_pruned);
                }
            }
        }

        if let Some(final_diag) = check_termination(sequence, graph, &mut path_wavefronts, d, modality) {
            final_diagonal = final_diag;
            break;
        }

        //Sperimental for pruning heuristics
        
        if prune_condition(graph, sequence, path, path_wavefronts.diagonals_queue.size(), max_penalty, d) {
            diagonals_pruned += prune(&mut path_wavefronts, d, modality, max_delta_allowed);
        }        
        
        expand::<T>(sequence, graph, &mut path_wavefronts, &mut mem_set, d, m, ins, del, wavefront_impl);

        d += 1;
    }

    let mut retrieve_traceback = false;

    let mut maybe_max_score = maybe_max_score_mutex.lock().unwrap();
    if let Some(max_score) = *maybe_max_score {
        if d < max_score || (find_all_alignments && d == max_score) {
            *maybe_max_score = Some(d);
            retrieve_traceback = true;
        }
    }
    else {
        *maybe_max_score = Some(d);
        retrieve_traceback = true;
    }
    drop(maybe_max_score);      //Release the lock
    
    if retrieve_traceback {
        traceback(
            sequence, 
            graph, 
            &mut path_wavefronts, 
            final_diagonal, 
            d,
            m,
            ins,
            del,
            diagonals_pruned
        )
    }
    else {
        Alignment::new(path, 0, 0, 0, d, diagonals_pruned)
    }
}

/// Runs **WFA** to provide an **optimal alignment** between a <code>graph</code> and a <code>sequence</code>. 
/// # Arguments
/// - <code>graph</code>: **multiple string** rapresentation of the **graph** to align;
/// - <code>sequence</code>: **sequence** to align;
/// - <code>optimal_alignments</code>: reference to a <code>Vec</code> where the **optimal alignment(s)** 
/// will be stored;
/// - <code>m</code>: **mismatch penalty**;
/// - <code>ins</code>: **insertion penalty**;
/// - <code>del</code>: **deletion penalty**;
/// - <code>wavefront_impl</code>: **implementation** to store the wavefronts;
/// - <code>modality</code>: **alignment modality** to use;
/// - <code>parallelize_match</code>: if <code>true</code>, every diagonal will be extended by a different **thread**;
/// otherwise, all the diagonals will be extended **sequentially**;
/// - <code>max_delta_allowed</code>: **Sperimental for pruning heuristics**. Pruning maximum distance allowed;
/// - <code>prune_condition</code>: **Sperimental for pruning heuristics**. Function to decide wheter perform pruning or not;
/// - <code>find_all_alignments</code>: if <code>true</code>, all the best path alignments are returned; otherwise just one; 
/// - <code>max_threads</code>: max number of **threads running in parallel**.
/// # Return value
/// Returns the **value** of the **alignment**; the **sequence(s) of operations** performed are stored in 
/// <code>optimal_alignments</code>.
///
/// # Complexity
/// Let
/// - <code>n</code>: the length of the **path** which aligns better to the **sequence**,
/// - <code>p</code>: the numbers of path,
/// - <code>m</code>: the length of the **sequence**,
/// - <code>d</code>: the **best alignment score**
/// ## Time
/// The alignment to each path is performed by a different **thread** which stops when 
/// the best solution is found; consequentially, assuming every **thread running in parallel** (for a low
/// number of path), the **time complexity** is <code>O(max{n, m} d)</code>. Otherwise, considering an 
/// arbitrary number of paths, the **time complexity** becomes <code>O(p max{n, m} d)</code>.
/// ## Space
/// The **space complexity** depends on the wavefront implementation chosen; however, for each implementation,
/// in the **worst case** every diagonal is stored in the wavefront for **every score** lower than the 
/// **optimal alignment value**, so the **worst case space complexity** is <code>O(p (n + m) d)</code>.
pub fn wf_align<T>(
    graph: &PathStrings,
    sequence: &[char],
    optimal_alignments: &mut Vec<Alignment>,
    m: usize,
    ins: usize,
    del: usize,
    wavefront_impl: WavefrontImpl,
    modality: AlignmentMod,
    parallelize_match: bool,
    max_delta_allowed: usize,
    prune_condition: fn(&PathStrings, &[char], usize, usize, usize, usize) -> bool,
    find_all_alignments: bool,
    max_threads: usize,
) -> usize 
where T: num::NumCast + std::cmp::Eq + Hash + Copy + 'static {

    let mut threads_alignments: Vec<Alignment> = Vec::with_capacity(graph.paths_number);
    let maybe_max_score = None;
    let maybe_max_score_mutex = Arc::new(Mutex::new(maybe_max_score));

    thread::scope(|s| {
        let (tx, rx) = mpsc::channel();

        let mut path = 0;    
        let mut threads_running = 0;
        let mut alignment_terminated = 0;

        while alignment_terminated < graph.paths_number {
            while threads_running < max_threads && path < graph.paths_number {
                let tx_thread = tx.clone();
                let thread_maybe_max_score_mutex = Arc::clone(&maybe_max_score_mutex);

                s.spawn(move || {
                    let alignment = wf_align_to_path::<T>(
                        sequence,
                        graph,
                        path,
                        m,
                        ins,
                        del,
                        wavefront_impl,
                        modality,
                        parallelize_match,
                        max_delta_allowed,
                        prune_condition,
                        find_all_alignments,
                        thread_maybe_max_score_mutex,
                    );
                    tx_thread.send(alignment).unwrap();
                });

                threads_running += 1;
                path += 1;
            }
                
            match rx.recv() {
                Ok(alignment) => {
                    if alignment.score == (*maybe_max_score_mutex.lock().unwrap()).unwrap() {
                        threads_alignments.push(alignment);
                    }
                },
                Err(_) => (),
            }

            alignment_terminated += 1;
            threads_running -= 1;
        }
    });

    let mut i = 0; 
    let max_score = (*maybe_max_score_mutex.lock().unwrap()).unwrap();

    if find_all_alignments {
        while i < threads_alignments.len() {
            if threads_alignments[i].score == max_score {
                optimal_alignments.push(threads_alignments.swap_remove(i));
            }
            else {
                i += 1;
            }
        }
    }
    else {
        while i < threads_alignments.len() && optimal_alignments.len() == 0 {
            if threads_alignments[i].score == max_score &&
                threads_alignments[i].path_length > 0  {
                optimal_alignments.push(threads_alignments.swap_remove(i));
            }
            else {
                i += 1;
            }
        }
    }

    max_score
}