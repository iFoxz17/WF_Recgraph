/// Wavefront interface definition.
mod wavefront;
/// Wavefront implementations.
mod wf_implementation;
/// Wavefront alignment implementation.
mod wf_alignment;

extern crate num;

use crate::pathwise_graph::PathGraph;
use crate::gaf_output::*;
use crate::utils;

use wf_alignment::*;
use wf_implementation::*;

fn path_graph_to_path_strings(pred_graph: &PathGraph) -> PathStrings {

    let mut paths_mapping = Vec::with_capacity(pred_graph.paths_number);
    
    for _path in 0..pred_graph.paths_number { 
        paths_mapping.push(Vec::with_capacity(pred_graph.lnz.len()));
    }

    for v in 0..pred_graph.paths_nodes.len() {
        for path in 0..pred_graph.paths_nodes[v].len() {
            if pred_graph.paths_nodes[v][path] {
                paths_mapping[path].push(v);
            }
        }
    }

    PathStrings::new(pred_graph.lnz.clone(), paths_mapping, pred_graph.paths_number)
}

fn alignment_to_gaf_struct(sequence: &[char], path_graph: &PathGraph, alignment: &Alignment) -> GAFStruct {
    let comments = format!(
        "{}, best path: {}, score: {}\t{}",
        alignment.get_cigar(),
        alignment.path,
        alignment.score,
        alignment.get_path_string(),
    );

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set

    let mut sequence_path_nodes = Vec::with_capacity(alignment.path_nodes.len());
    let mut last_sequence_node = 0;
    let mut first_iteration = true;

    for v in &alignment.path_nodes {
        if path_graph.nodes_id_pos[*v] != last_sequence_node || first_iteration {
            last_sequence_node = path_graph.nodes_id_pos[*v];
            sequence_path_nodes.push(last_sequence_node as usize);
            first_iteration = false;
        }
    }

    let mut final_offset = 0;

    if sequence_path_nodes[sequence_path_nodes.len() - 1] == 0 {
        sequence_path_nodes.remove(sequence_path_nodes.len() - 1);
        final_offset = 1;
    }
 
    let (path_len, path_start, path_end) = utils::get_path_len_start_end(
        &path_graph.nodes_id_pos,
        alignment.path_start + if alignment.path_start == 0 {0} else {1},
        alignment.path_end - final_offset,
        alignment.path_length - final_offset,
    );

    let gaf = GAFStruct::build_gaf_struct(
        String::from("Temp"),
        sequence.len() - 1,
        0,
        sequence.len() - 2,
        '+',
        sequence_path_nodes,
        path_len - 1,
        path_start,
        path_end - 1,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    gaf
}

fn choose_uint_size(graph: &PathGraph, sequence: &[char]) -> usize {

    let max = |x, y| {if x >= y {x} else {y}};
    let choose_usize = false;

    if choose_usize {
        0
    }
    else if max(graph.lnz.len(), sequence.len()) < usize::pow(2, 7) - 1 {
        8
    } 
    else if max(graph.lnz.len(), sequence.len()) < usize::pow(2, 15) - 1 {
        16
    }
    else if max(graph.lnz.len(), sequence.len()) < usize::pow(2, 31) - 1 {
        32
    }
    else if max(graph.lnz.len(), sequence.len()) < usize::pow(2, 63) - 1 {
        64
    }
    else {
        128
    }
}

/// Runs **multiple threads WFA** to provide a **global optimal alignment** between a 
/// <code>canonical variation graph</code> and a <code>sequence</code>. 
/// # Arguments
/// - <code>sequence</code>: **sequence** to align;
/// - <code>graph</code>: **multiple string** rapresentation of the **graph** to align;
/// - <code>m</code>: **mismatch penalty**;
/// - <code>ins</code>: **insertion penalty**;
/// - <code>del</code>: **deletion penalty**;
/// # Return value
/// Returns a <code>GAFStruct</code> which contains all the alignment informations.
pub fn wf_pathwise_alignment_global(
    sequence: &[char],
    path_graph: &PathGraph,
    m: usize,
    ins: usize,
    del: usize,
) -> GAFStruct {

    let graph = path_graph_to_path_strings(path_graph);
    let mut optimal_alignments = Vec::new();
    let wf_implementation = WavefrontImpl::WavefrontVec;
    //let wf_implementation = WavefrontImpl::WavefrontHash;
    //let wf_implementation = WavefrontImpl::WavefrontMixed(0.7);
    let parallelize_match = false;
    let uint_type = choose_uint_size(path_graph, sequence);

    let _d = match uint_type {
        8 => wf_align::<u8>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match),
        16 => wf_align::<u16>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match),
        32 => wf_align::<u32>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match),
        64 => wf_align::<u64>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match),
        128 => wf_align::<u128>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match), 
        _ => wf_align::<usize>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Global, parallelize_match),
    };

    alignment_to_gaf_struct(sequence, path_graph, &optimal_alignments[0])
}

/// Runs **multiple threads WFA** to provide a **semiglobal optimal alignment** between a 
/// <code>canonical variation graph</code> and a <code>sequence</code>. 
/// # Arguments
/// - <code>sequence</code>: **sequence** to align;
/// - <code>graph</code>: **multiple string** rapresentation of the **graph** to align;
/// - <code>m</code>: **mismatch penalty**;
/// - <code>ins</code>: **insertion penalty**;
/// - <code>del</code>: **deletion penalty**;
/// # Return value
/// Returns a <code>GAFStruct</code> which contains all the alignment informations.
pub fn wf_pathwise_alignment_semiglobal(
    sequence: &[char],
    path_graph: &PathGraph,
    m: usize,
    ins: usize,
    del: usize,
) -> GAFStruct {

    let graph = path_graph_to_path_strings(path_graph);
    let mut optimal_alignments = Vec::new();
    //let wf_implementation = WavefrontImpl::WavefrontVec;
    //let wf_implementation = WavefrontImpl::WavefrontHash;
    let wf_implementation = WavefrontImpl::WavefrontMixed(0.3);
    let parallelize_match = false;
    let uint_type = choose_uint_size(path_graph, sequence);

    let _d = match uint_type {
        8 => wf_align::<u8>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match),
        16 => wf_align::<u16>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match),
        32 => wf_align::<u32>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match),
        64 => wf_align::<u64>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match),
        128 => wf_align::<u128>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match), 
        _ => wf_align::<usize>(&graph, &sequence, &mut optimal_alignments, m, ins, del, 
            wf_implementation, AlignmentMod::Semiglobal, parallelize_match),
    };

    alignment_to_gaf_struct(sequence, path_graph, &optimal_alignments[0])
}

