/// Wavefront interface definition.
mod wavefront;
/// Wavefront implementations.
mod wf_implementation;
/// Wavefront alignment implementation.
mod wf_alignment;

use crate::pathwise_graph::PathGraph;
use crate::gaf_output::*;

use wf_alignment::*;
use wf_implementation::*;

//Paths_nodes associa ad ogni nodo un bitvec dei paths a cui appartiene
//nodes_is_pos mappa ogni nodo della linearizzazione nel nodo a cui apparteneva nel file gaf iniziale
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

fn alignment_to_gaf_struct(sequence: &[char], alignment: &Alignment) -> GAFStruct {
    let comments = format!(
        "{}, best path: {}, score: {}\t{}",
        alignment.get_cigar(),
        alignment.path,
        alignment.score,
        alignment.get_path_string(),
    );

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set

    let gaf = GAFStruct::build_gaf_struct(
        String::from("Temp"),
        sequence.len(),
        0,
        sequence.len() - 1,
        '+',
        alignment.path_nodes.clone(),
        alignment.path_length,
        alignment.path_start,
        alignment.path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    gaf
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

    let _d = wf_align(
        &graph, 
        sequence, 
        &mut optimal_alignments, 
        m, 
        ins, 
        del, 
        wf_implementation, 
        AlignmentMod::Global, 
        parallelize_match
    );

    alignment_to_gaf_struct(sequence, &optimal_alignments[0])
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

    let _d = wf_align(
        &graph, 
        sequence, 
        &mut optimal_alignments, 
        m, 
        ins, 
        del, 
        wf_implementation, 
        AlignmentMod::Semiglobal, 
        parallelize_match
    );

    alignment_to_gaf_struct(sequence, &optimal_alignments[0])
}

