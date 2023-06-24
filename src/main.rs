use recgraph::args_parser;
use recgraph::gap_global_abpoa;
use recgraph::gap_local_poa;
use recgraph::global_abpoa;
use recgraph::graph;
use recgraph::local_poa;
use recgraph::pathwise_alignment;
use recgraph::pathwise_alignment_gap;
use recgraph::pathwise_alignment_gap_semi;
use recgraph::pathwise_alignment_recombination;
use recgraph::pathwise_alignment_semiglobal;
use recgraph::pathwise_graph;
use recgraph::pathwise_graph::nodes_displacement_matrix;
use recgraph::score_matrix;
use recgraph::sequences;
use recgraph::utils;

use recgraph::wfa;

use std::collections::HashMap;

use std::time::SystemTime;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

fn main() {
    let now = SystemTime::now();

    // get sequence
    let (sequences, seq_names) = sequences::get_sequences(args_parser::get_sequence_path());

    //get graph
    let graph_path = args_parser::get_graph_path();
    let graph_struct = graph::read_graph(&graph_path, false);

    //get score matrix
    let score_matrix = score_matrix::create_score_matrix();
    let scores_f32 = score_matrix::create_f32_scores_matrix();

    //get alignment option
    let align_mode = args_parser::get_align_mode();
    let amb_strand = args_parser::get_amb_strand_mode();
    let (b, f) = args_parser::get_b_f();

    //get handle position for output
    let hofp_forward = utils::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, false);
    let mut hofp_reverse = HashMap::new();

    match align_mode {
        //global alignment
        0 => {
            let r_values = utils::set_r_values(
                &graph_struct.nwp,
                &graph_struct.pred_hash,
                graph_struct.lnz.len(),
            );
            for (i, seq) in sequences.iter().enumerate() {
                let bases_to_add = (b + f * seq.len() as f32) as usize;
                let alignment = if is_x86_feature_detected!("avx2") {
                    unsafe {
                        global_abpoa::exec_simd(
                            seq,
                            (&seq_names[i], i + 1),
                            &graph_struct,
                            &scores_f32,
                            bases_to_add,
                            false,
                            &hofp_forward,
                            &r_values,
                        )
                    }
                } else {
                    global_abpoa::exec(
                        seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        bases_to_add,
                        false,
                        &hofp_forward,
                    )
                };
                if amb_strand && alignment.0 < 0 {
                    if hofp_reverse.is_empty() {
                        hofp_reverse =
                            utils::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, true);
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    let rev_alignment = global_abpoa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        bases_to_add,
                        true,
                        &hofp_reverse,
                    );
                    if rev_alignment.0 > alignment.0 {
                        utils::write_gaf(&rev_alignment.1.unwrap().to_string(), i + 1);
                    } else {
                        utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                    }
                } else {
                    utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                }
            }
        }
        //local alignment
        1 => {
            for (i, seq) in sequences.iter().enumerate() {
                let alignment = if is_x86_feature_detected!("avx2") {
                    unsafe {
                        let temp_score = local_poa::exec_simd(
                            seq,
                            (&seq_names[i], i + 1),
                            &graph_struct,
                            &scores_f32,
                            false,
                            &hofp_forward,
                        );
                        (temp_score.0 as i32, temp_score.1)
                    }
                } else {
                    local_poa::exec(
                        seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        false,
                        &hofp_forward,
                    )
                };
                if amb_strand {
                    if hofp_reverse.is_empty() {
                        hofp_reverse =
                            utils::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, true);
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    let alignment_rev = if is_x86_feature_detected!("avx2") {
                        unsafe {
                            let temp_alignment = local_poa::exec_simd(
                                &rev_seq,
                                (&seq_names[i], i + 1),
                                &graph_struct,
                                &scores_f32,
                                true,
                                &hofp_reverse,
                            );
                            (temp_alignment.0 as i32, temp_alignment.1)
                        }
                    } else {
                        local_poa::exec(
                            &rev_seq,
                            (&seq_names[i], i + 1),
                            &graph_struct,
                            &score_matrix,
                            true,
                            &hofp_reverse,
                        )
                    };
                    if alignment.0 < alignment_rev.0 {
                        utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                    } else {
                        utils::write_gaf(&alignment_rev.1.unwrap().to_string(), i + 1);
                    }
                } else {
                    utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1)
                }
            }
        }
        //affine gap global alignment
        2 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();

            for (i, seq) in sequences.iter().enumerate() {
                let bases_to_add = (b + f * seq.len() as f32) as usize;
                let alignment = gap_global_abpoa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    bases_to_add,
                    false,
                    &hofp_forward,
                );

                if amb_strand && alignment.0 < 0 {
                    if hofp_reverse.is_empty() {
                        hofp_reverse =
                            utils::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, true);
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    let rev_alignment = gap_global_abpoa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        bases_to_add,
                        true,
                        &hofp_reverse,
                    );
                    if rev_alignment.0 > alignment.0 {
                        utils::write_gaf(&rev_alignment.1.unwrap().to_string(), i + 1);
                    } else {
                        utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                    }
                } else {
                    utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                }
            }
        }
        //affine gap local alignment
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let alignment = gap_local_poa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    false,
                    &hofp_forward,
                );
                if amb_strand {
                    if hofp_reverse.is_empty() {
                        hofp_reverse =
                            utils::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, true);
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    let rev_alignment = gap_local_poa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        false,
                        &hofp_reverse,
                    );
                    if rev_alignment.0 > alignment.0 {
                        utils::write_gaf(&rev_alignment.1.unwrap().to_string(), i + 1);
                    } else {
                        utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                    }
                } else {
                    utils::write_gaf(&alignment.1.unwrap().to_string(), i + 1);
                }
            }
        }
        4 => {
            let graph = pathwise_graph::read_graph_w_path(&graph_path, false);

            let recgraph_mod = false;
            let print_status = true;
            let mut part_elapsed = 0;
            let mut last_partial;
            let (m, ins, del) = (1, 3, 2);

            let mut times = vec![0];

            if print_status {
                eprintln!("Starting alignment with options: ");
                eprint!("Algorithm used: ");
                if recgraph_mod {
                    eprintln!("recgraph pathwise alignment");
                }
                else {
                    eprintln!("wavefront multithread pathwise alignment: m={m}, ins={ins}, del={del}");
                }	
                eprintln!("Graph linearization length: {}", graph.lnz.len());
                eprintln!("Graph number of paths: {}", graph.paths_number);
                eprintln!("Number of reads: {}", sequences.len());
                eprintln!("Modality: Global");
                eprintln!("Print_status = {}", print_status);
            }

            if recgraph_mod {
                for (i, seq) in sequences.iter().enumerate() {
                    if print_status {
                        eprint!("Processing {}/{} ({:.1}%) -> ", 
                                i + 1, sequences.len(),
                                (i + 1) as f64 / sequences.len() as f64 * 100.0 
                        );
                    }
                    let mut gaf = pathwise_alignment::exec(seq, &graph, &score_matrix);
                    gaf.query_name = seq_names[i].clone();
                    utils::write_gaf(&gaf.to_string(), i);
                    if print_status {
                        last_partial = part_elapsed;
                        part_elapsed = now.elapsed().unwrap().as_millis();
                        times.push(part_elapsed - last_partial);
                        
                        eprintln!("Done: {:.1} s from start, expecting {:.1} s, {:.1} s remaining, mean for alignment {:.1} s, delta {:.1}", 
                                part_elapsed as f64 / 1000.0,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * sequences.len() as f64,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * (sequences.len() - i - 1) as f64,
				                part_elapsed as f64 / 1000.0 / (i + 1) as f64,
                                times[times.len() - 1] as f64 / 1000.0
                        );
                    }
                }
            }
            else {
                for (i, seq) in sequences.iter().enumerate() {
                    if print_status {
                        eprint!("Processing {}/{} ({:.1}%) -> ", 
                                i + 1, sequences.len(),
                                (i + 1) as f64 / sequences.len() as f64 * 100.0 
                        );
                    }
                    let mut gaf = wfa::wf_pathwise_alignment_global(
                        seq, 
                        &graph,
                        &wfa::path_graph_to_path_strings(&graph), 
                        m, 
                        ins, 
                        del
                    );
                    gaf.query_name = seq_names[i].clone();
                    utils::write_gaf(&gaf.to_string(), i);
                    if print_status {
                        last_partial = part_elapsed;
                        part_elapsed = now.elapsed().unwrap().as_millis();
                        times.push(part_elapsed - last_partial);
                        
                        eprintln!("Done: {:.1} s from start, expecting {:.1} s, {:.1} s remaining, mean {:.1} s, delta {:.1} s", 
                                part_elapsed as f64 / 1000.0,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * sequences.len() as f64,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * (sequences.len() - i - 1) as f64,
				                part_elapsed as f64 / 1000.0 / (i + 1) as f64,
                                times[times.len() - 1] as f64 / 1000.0,                                
                        );
                    }
                }
            }

            if print_status {
                eprintln!("----------------------------------------------------------------------------------");
                eprintln!("Stats:");
                eprintln!("----------------------------------------------------------------------------------");
                eprintln!("Number of alignments: {}", times.len());
                eprintln!("Mean time for alignment: {:.1}", 
                statistical::mean(
                    &(times.iter().map(|x| *x as f64 / 1000.0).collect::<Vec<f64>>()), 
                ));
                eprintln!("SD time for alignment: {:.1}", 
                statistical::standard_deviation(
                    &(times.iter().map(|x| *x as f64 / 1000.0).collect::<Vec<f64>>()), 
                    None
                ));
            }
        }
        5 => {
            let graph = pathwise_graph::read_graph_w_path(&graph_path, false);

            let recgraph_mod = false;
            let print_status = true;
            let mut part_elapsed = 0;
            let mut last_partial;
	        let (m, ins, del) = (1, 1, 1);
            let mut times = vec![];

            if print_status {
                eprintln!("Starting alignment with options: ");
                eprint!("Algorithm used: ");
                if recgraph_mod {
                    eprintln!("recgraph pathwise alignment");
                }
                else {
                    eprintln!("wavefront multithread pathwise alignment: m={m}, ins={ins}, del={del}");
                }
                eprintln!("Graph linearization length: {}", graph.lnz.len());
                eprintln!("Graph number of paths: {}", graph.paths_number);
                eprintln!("Number of reads: {}", sequences.len());
                eprintln!("Modality: Semiglobal");
                eprintln!("Print_status = {}", print_status);
            }

            if recgraph_mod {
                for (i, seq) in sequences.iter().enumerate() {
                    if print_status {
                        eprint!("Processing {}/{} ({:.1}%) -> ", 
                                i + 1, sequences.len(),
                                (i + 1) as f64 / sequences.len() as f64 * 100.0 
                        );
                    }
                    let mut gaf = pathwise_alignment_semiglobal::exec(seq, &graph, &score_matrix);
                    gaf.query_name = seq_names[i].clone();
                    utils::write_gaf(&gaf.to_string(), i);
                    if print_status {
                        last_partial = part_elapsed;
                        part_elapsed = now.elapsed().unwrap().as_millis(); 
                        times.push(part_elapsed - last_partial);
                        
                        eprintln!("Done: {:.1} s from start, expecting {:.1} s, {:.1} s remaining, mean for alignment {:.1} s, delta {:.1}", 
                                part_elapsed as f64 / 1000.0,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * sequences.len() as f64,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * (sequences.len() - i - 1) as f64,
				                part_elapsed as f64 / 1000.0 / (i + 1) as f64,
                                times[times.len() - 1] as f64 / 1000.0
                        );
                    }

                }
            }
            else {
                for (i, seq) in sequences.iter().enumerate() {
                    if print_status { 
                        eprint!("Processing {}/{} ({:.1}%) -> ", 
                                i + 1, sequences.len(),
                                (i + 1) as f64 / sequences.len() as f64 * 100.0 
                        );
                    }
                    let mut gaf = wfa::wf_pathwise_alignment_semiglobal(
                        seq, 
                        &graph,
                        &wfa::path_graph_to_path_strings(&graph), 
                        m, 
                        ins, 
                        del
                    );
                    gaf.query_name = seq_names[i].clone();
                    utils::write_gaf(&gaf.to_string(), i);
                    if print_status {
                        last_partial = part_elapsed;
                        part_elapsed = now.elapsed().unwrap().as_millis(); 
                        times.push(part_elapsed - last_partial);
                        
                        eprintln!("Done: {:.1} s from start, expecting {:.1} s, {:.1} s remaining, mean for alignment {:.1} s, delta {:.1}", 
                                part_elapsed as f64 / 1000.0,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * sequences.len() as f64,
                                part_elapsed as f64 / 1000.0 / (i + 1) as f64 * (sequences.len() - i - 1) as f64,
				                part_elapsed as f64 / 1000.0 / (i + 1) as f64,
                                times[times.len() - 1] as f64 / 1000.0
                        );
                    }
                }
            }
            if print_status {
                eprintln!("----------------------------------------------------------------------------------");
                eprintln!("Stats:");
                eprintln!("----------------------------------------------------------------------------------");
                eprintln!("Number of alignments: {}", times.len() - 1);
                eprintln!("Mean time for alignment: {:.1}", 
                statistical::mean(
                    &(times.iter().map(|x| *x as f64 / 1000.0).collect::<Vec<f64>>()), 
                ));
                eprintln!("SD time for alignment: {:.1}", 
                statistical::standard_deviation(
                    &(times.iter().map(|x| *x as f64 / 1000.0).collect::<Vec<f64>>()), 
                    None
                ));
            }
        }
        6 => {
            let graph = pathwise_graph::read_graph_w_path(&graph_path, false);
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let best_path =
                    pathwise_alignment_gap::exec(seq, &graph, &score_matrix, g_open, g_ext);
                println!("Best path sequence {i}: {best_path}");
            }
        }
        7 => {
            let graph = pathwise_graph::read_graph_w_path(&graph_path, false);
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let best_path =
                    pathwise_alignment_gap_semi::exec(seq, &graph, &score_matrix, g_open, g_ext);
                println!("Best path sequence {i}: {best_path}");
            }
        }
        8 | 9 => {
            let graph = pathwise_graph::read_graph_w_path(&graph_path, false);
            let rev_graph = pathwise_graph::create_reverse_path_graph(&graph);
            let displ_matrix = nodes_displacement_matrix(&graph, &rev_graph);

            let (base_rec_cost, multi_rec_cost) = args_parser::get_base_multi_recombination_cost();
            let rbw = args_parser::get_recombination_band_width();

            for (i, seq) in sequences.iter().enumerate() {
                let mut gaf = pathwise_alignment_recombination::exec(
                    align_mode,
                    seq,
                    &graph,
                    &rev_graph,
                    &score_matrix,
                    base_rec_cost,
                    multi_rec_cost,
                    &displ_matrix,
                    rbw,
                );
                gaf.query_name = seq_names[i].clone();

                utils::write_gaf(&gaf.to_string(), i);
            }
        }

        _ => {
            panic!("Alignment mode must be in [0..9]");
        }
    }
    match now.elapsed() {
        Ok(elapsed) => {
            // it prints '2'
            eprintln!("Done in {} seconds.", elapsed.as_millis() as f64 / 1000.0);
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {e:?}");
        }
    }
}

#[cfg(test)]
mod tests {

    use std::time::*;
    use crate::*;

    #[test]
    fn test_pathwise_wfa_vs_recgraph_global() {

        // get sequence
        let (sequences, seq_names) = sequences::get_sequences("example\\reads.fa".to_string());

        //get graph
        let graph_struct = graph::read_graph("example\\graph.gfa", false);

        //get score matrix
        let score_matrix = score_matrix::create_score_matrix();

        let graph = pathwise_graph::read_graph_w_path("example\\graph.gfa", false);

        let mut recgraph_time = 0;
        let mut wf_time = 0;
        let mut now;

        now = Instant::now();
        for (i, seq) in sequences.iter().enumerate() {
            let mut gaf = wfa::wf_pathwise_alignment_global(seq, &graph, 1, 2, 3);
            gaf.query_name = seq_names[i].clone();
            utils::write_gaf(&gaf.to_string(), i);
        }
        wf_time += now.elapsed().as_millis();

        println!("Wavefront time: {}", wf_time as f64 / 1000.0);

        println!("Starting recgraph alignment");

        now = Instant::now();
        for (i, seq) in sequences.iter().enumerate() {
            let mut gaf = pathwise_alignment::exec(seq, &graph, &score_matrix);
            gaf.query_name = seq_names[i].clone();
            utils::write_gaf(&gaf.to_string(), i);
        }
        recgraph_time += now.elapsed().as_millis();
        
        println!("Recgraph time: {}", recgraph_time as f64 / 1000.0);
    }
}
