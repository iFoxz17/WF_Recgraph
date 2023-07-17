use std::collections::HashMap;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use recgraph::{global_abpoa, graph, local_poa, score_matrix, utils};

fn bench_local_simd_no_simd(c: &mut Criterion) {
    let seq = "TGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATAGGCAGGCCACTCTGGGATCCCTGGGATGCAAGCCCAGGGACAGCAGAGTCCCCAGGTGGGAAATCTACACACACACCCCAGGGATGTCCCAGAGACTTCTACCCTAAGAGGAGATCCTGGGCAGGATGTGAGAAATCTGAGCATCCTCTGTTTGGATGGCCGAAGCTGCTGGCATCAAACTCTGGTCTGGAAGAATCAGTCTGGGGGAGAGACAGGGATGGAGGAAAGGCATCAGGGGATCCATCCTCCTCCTCCTTCTCCTCCTCCTCCTCCCCCACAAAGGCCTTGCTCGCCCTGCCTGCACCACACCCTGCAGAAGTTGATCTCTCCTTGTTCCCAAATCATCTCCAAGCACCCTTCCTACAGCACCCCATGATTCCTTTTTTCACTCAAAGCAATTCTTGTGACCCATAACTGTGTGTGTGTAACTGGGTCCCCAACTGGGAAGATGTGCCCCCATGGTGCTGGATACAGGCCCCCACACCCAAGGGCCTGAGGATCGCTATATGTCCCCCCATGCCACAAAATAATCCTGACACATGCACGCATGCACCACTGTATCTGGCTCCCACAGGCTCACCCGCCCCCTCCAGATGACATACCACCTGAGCAAGGCTTCCGGAAGTAGATGATGAGAACAATGCCCACGATGATGCCCAGCACACCCAGGCCAAAGGCCACGCCACACAGCACATTCTCCAGCAGATCTGAGGGCAGTGCGTTCCGGGGTACTGGAGGAAATGAGTGGCTCAGCCTGGGGACCTAGTTAGGGAGCCTCCCACCCAGGGAAATGACGTGGGTGTCTGGGATGACATGGGAGACTGGGATGGGCTTAGGGTAGGAATGGACTAAACAAGGTACCAGTGGAGAAAGAAGCCTCCTCCCATGGATCTATCCCTTTTTGCCCCCAAAAGGACCAGAATTCCAGGGAGAAAGCCTCACCCCAATAGGCAATTGCTGTGTAGCGGTCAATTTCGTGAGTCACAATGCAGGAGAAAATGTCAGAAGGTTCTGGTGTGAAGTTTAAGTAAGAAAAGGCCTGGAAGCTGAGTCCATCGACAGCTGAGACAAAAGTAGGCCCAAATCCTTCCACAGGGACGGAATGATGCTGCCAGTTCACTGTCAGCATGGGTGGGAAGAGATTACTGACAAAACAGACCAAAGTGTTGGGCTTGCCAAACTCCAGGGGCTTCAGCGTGAACACTTCAGCGATAGGAAACCCTGGTGGGGGGATTGAAGTGTAGGGGGAAAAAGAGACTAGTTTAGATGGTATCTCTGTGTTTGGAGGGGCCATGGCATATGGAGGGGAGGGCAGAGAAGAACACAGTGGGTCAGGCTTTGGGAGACAGAGATGAGCGAGGAGCTGGGCTCTGAAGGGAGGTCTTCTTCCAGGCAAGGACTGCAGCTAGACATAGAAGCAGAGCCAGATCCAGGCTACTCTGGACCCCTCCACCATGACTTCCTTCAGCACTTCCTGTCTAGAGCTCACATTGATGTCTAACCATGCACTGTCTTCTCACTAAGACATAGTCACGTCATCAGATATTTCCACTCTTCCCATCCATCTTGCTGGGCATAGTAGCACAAGTGTTAATATTCAGTAGGTATCAGTTGGTACCTGTTGAATTCATCACATTCAATACATAGTTCTGAATGCCTACTACATGCTAGGTACTTCGGCCCACCAAAAGAACACAGGGTGCAGACCAAGGCTGGTGGAAAAATTAAGGTGATGAAGAGAACCAGAAAGTATTTGAGATGGGGAGCTGGTATCAAGGGGAATTATTCAGTGTACAGATCAATGAGGTTAATGCAGCCCTCCTCCCTTCACTCCCCAGAAAACTCCTGACCTCTGGACACCGGGATTTTCCCATCAAGTTTTGGCCCTATTTGCTGGATCATCCACTCGCAGAACTCTTTGTCAAATAAAATGGCAGGAGCATCTCCCTGTTCCTGAGCCCAGTCAGCAAATTCGGGCAGGCGAGGCACCCGAGTGTTCTGGGAAAAGTCGAAGAAGAAAAGCTGGTCCTCGTCGTAGGCCTCAGAGAGTCCCACACTGGGACTCCCATCCTGGCAGTACACTGTGTGCAGGAATGTGTGGTTTTGCAGGTCATCTGGCCACATTGGAGTAGGAGCTGCAAAGGACACAGGGTGAGGTTCAGGGAGGTGGGAGCCTTCTCCTCCAACTTAAAAAACAGCAAGGTGGGGCTAGGCGCAGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGTGGATCATGAGGTCAGGAGTTTGAGACCAGCCTGGCCAGCATGGTGAAACTCCGTCTCTACTAAAAATACAAAAAAGTAGCTGGGCATGTTGGCATGCGCCTGTAGCTACTCGGGAGGCTGAGGGAGGAGAATTGCTTGAACCAGGGAGGCAGAGGTTGCCGGGAGCTAAGATTAAGCCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCTGTCTCAAAACAAAACAACAAAAACAAGCAAGGCCTGCTTAAGGAGCGTGGGCTGAGGTGAGACCCTTTCCTGTGTCTGTTATTTAGACTCCCCCTCCCAAAGGGGGTGAAGAACAAATTATGGCATCTCTCCAAGCTTCCCCTGCCTATAAAAAGGCCAGTTGGCAAAAGTAAAGAGTTCTACTTTCTAAAGTGACAGATTCAGGCCAGGCATGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAGGCAGATTGCTTGAGCCCAGGAGTTCAAGACCAACCTGGGCAACACAGCGAGACCCTGTCTCTACAAAAAATACAAAAACTTAGCCAGGTGTGGTGGCAAACACCTGTGGTCTCAGCTACTCTGGAGGCTGAGGCAGGAGGATTGCTTGTGCCTAGGAAGTTGGGGCTGCAGTGAGCCATGATTGTGCCACTGGACTCCAGCCCAGGTGACAGAATGAGCCCGTCTCAAAAAATATATATATAAAGGCCGGGCGCGGTGGCTCAAGCTTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCAAACATGATGAAACCCCATCTCTACTAAAAATACAAAAATCAGCTGGGTGTGGTGGCATGCGCCTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAGTCTCTTGAACCCCAGAGGCAGGGGTTGCAGGGAGCCGAGATCACGTCACTGCACTCTACCCTGGGTGACAGAGCGAGATGCCGTGTCAAAAAAAATAAATTAAATCAAATAAAAAATTTAAAAATGTATATATATAAAATAAAGTGACAGATTCAGAGTCACTGTTCATTGTGTGTTTGGGGGCTGCACAAAGACACCTAGCCAAAGAAGCAAGTGAAAGCCTGCATTCTGCTCACCATGCCATACATCCTGGCATAGGGCTGTATCCTCCCAAAGGGGATTCCTTTGTCTAATTCATACCAGGCCACTGTATTGACTAGAGAAGGCCATGGATGGGTTTCTCACTCTTAGAAGGGAAAGAGGAGGAATGGCTACAGCCTCCCCAAGCCATAGATGGGACTGCCTCCCACTATCCCCAGACACAAATGGTAAATTGGAAAACCTGTATCCAGACATTTCTTCAGCCACTTCATTGGCACCAAGCGTCTCTCAAAATGTCTTCTGTTCCTTAACCTACCAGGCCTCCCAAAGACAGCAATGGGAGAAGTGACCCCATAACTGCATAAAATAATCCCTCTTCTTTGAAGCTCTTGGCAGGAATCGCTCAGCCAGCAGGAAACCTTTAACCCAATACCCAGAAAAACAGACATTTGGAGGAAGAGGGATCTTCCAGATTATTCTTCCATTCTGCCCCATCCTCTACAGAGAAGGAAACTAAGACACTTTTCAAGAATCACAAGATAAGTTAATGATAGAAAGCAGAGTAGAATCTTGAGTGGAGGAGTGAAAATAACATTCACTTTGTTCAAATCCCAGCTCTACCACTTTCCAATGGTGTGAACTTGCACAAATAACTCTGAGTCTCATTTTCTTCATTTGTAAAATGGAGAGAACAATCTCCGCTTCAAGAGATTGTCTTAAATGGAACATGCAAAGCATCACTGATATCGTTTACCAACCACACATAGCAGCTGTCTTTCCCCACTCCCCTGTTGTTTCCACTGCCTCATAAGACTTCCCACCACTCACAAAGCACAGCGCTTTTCCTCACAAAGCTGAGTGGGCTCCCTAGGTTCAGGATGGAAGTAAATAGGAGTACCATCTTACCTTCAGGGACGGCCCAGGAGTGGGGTAGCAGCCACAGAAGTGGTAACATCTGTAGCAGCGCAGCTCCTTGGTTCTGTTCATGACCCATACCTTCTTGCCACACAGTAGGTAGGAGCTACCAACCCAGCCAACCCAGCTTCCCCAACTCCCTCCCCGAGAGGGTGGCCTTAGAT";
    let mut sequence = seq.chars().collect::<Vec<char>>();
    sequence.insert(0, '$');

    let graph_struct = graph::read_graph(
        &"tests/DMA-3108.fa.353ea42.34ee7b1.1576367.smooth.fix.gfa",
        false,
    );

    let score_matrix = score_matrix::create_score_matrix_match_mis(2, -4);
    let hofp = HashMap::new();
    let mut scores_f32 = HashMap::new();
    for (k, v) in score_matrix.iter() {
        scores_f32.insert(*k, *v as f32);
    }
    let mut group = c.benchmark_group("Local POA");
    for i in [1, 2, 3] {
        unsafe {
            group.bench_with_input(BenchmarkId::new(" Simd", i), &i, |b, _i| {
                b.iter(|| {
                    local_poa::exec_simd(
                        black_box(&sequence),
                        black_box(("name", 0)),
                        black_box(&graph_struct),
                        black_box(&scores_f32),
                        false,
                        &hofp,
                    )
                })
            });

            group.bench_with_input(BenchmarkId::new("No Simd", i), &i, |b, _i| {
                b.iter(|| {
                    local_poa::exec(
                        black_box(&sequence),
                        black_box(("name", 0)),
                        black_box(&graph_struct),
                        black_box(&score_matrix),
                        false,
                        &hofp,
                    )
                })
            });
        }
    }
}

fn bench_global_simd_no_simd(c: &mut Criterion) {
    let seq = "TGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATAGGCAGGCCACTCTGGGATCCCTGGGATGCAAGCCCAGGGACAGCAGAGTCCCCAGGTGGGAAATCTACACACACACCCCAGGGATGTCCCAGAGACTTCTACCCTAAGAGGAGATCCTGGGCAGGATGTGAGAAATCTGAGCATCCTCTGTTTGGATGGCCGAAGCTGCTGGCATCAAACTCTGGTCTGGAAGAATCAGTCTGGGGGAGAGACAGGGATGGAGGAAAGGCATCAGGGGATCCATCCTCCTCCTCCTTCTCCTCCTCCTCCTCCCCCACAAAGGCCTTGCTCGCCCTGCCTGCACCACACCCTGCAGAAGTTGATCTCTCCTTGTTCCCAAATCATCTCCAAGCACCCTTCCTACAGCACCCCATGATTCCTTTTTTCACTCAAAGCAATTCTTGTGACCCATAACTGTGTGTGTGTAACTGGGTCCCCAACTGGGAAGATGTGCCCCCATGGTGCTGGATACAGGCCCCCACACCCAAGGGCCTGAGGATCGCTATATGTCCCCCCATGCCACAAAATAATCCTGACACATGCACGCATGCACCACTGTATCTGGCTCCCACAGGCTCACCCGCCCCCTCCAGATGACATACCACCTGAGCAAGGCTTCCGGAAGTAGATGATGAGAACAATGCCCACGATGATGCCCAGCACACCCAGGCCAAAGGCCACGCCACACAGCACATTCTCCAGCAGATCTGAGGGCAGTGCGTTCCGGGGTACTGGAGGAAATGAGTGGCTCAGCCTGGGGACCTAGTTAGGGAGCCTCCCACCCAGGGAAATGACGTGGGTGTCTGGGATGACATGGGAGACTGGGATGGGCTTAGGGTAGGAATGGACTAAACAAGGTACCAGTGGAGAAAGAAGCCTCCTCCCATGGATCTATCCCTTTTTGCCCCCAAAAGGACCAGAATTCCAGGGAGAAAGCCTCACCCCAATAGGCAATTGCTGTGTAGCGGTCAATTTCGTGAGTCACAATGCAGGAGAAAATGTCAGAAGGTTCTGGTGTGAAGTTTAAGTAAGAAAAGGCCTGGAAGCTGAGTCCATCGACAGCTGAGACAAAAGTAGGCCCAAATCCTTCCACAGGGACGGAATGATGCTGCCAGTTCACTGTCAGCATGGGTGGGAAGAGATTACTGACAAAACAGACCAAAGTGTTGGGCTTGCCAAACTCCAGGGGCTTCAGCGTGAACACTTCAGCGATAGGAAACCCTGGTGGGGGGATTGAAGTGTAGGGGGAAAAAGAGACTAGTTTAGATGGTATCTCTGTGTTTGGAGGGGCCATGGCATATGGAGGGGAGGGCAGAGAAGAACACAGTGGGTCAGGCTTTGGGAGACAGAGATGAGCGAGGAGCTGGGCTCTGAAGGGAGGTCTTCTTCCAGGCAAGGACTGCAGCTAGACATAGAAGCAGAGCCAGATCCAGGCTACTCTGGACCCCTCCACCATGACTTCCTTCAGCACTTCCTGTCTAGAGCTCACATTGATGTCTAACCATGCACTGTCTTCTCACTAAGACATAGTCACGTCATCAGATATTTCCACTCTTCCCATCCATCTTGCTGGGCATAGTAGCACAAGTGTTAATATTCAGTAGGTATCAGTTGGTACCTGTTGAATTCATCACATTCAATACATAGTTCTGAATGCCTACTACATGCTAGGTACTTCGGCCCACCAAAAGAACACAGGGTGCAGACCAAGGCTGGTGGAAAAATTAAGGTGATGAAGAGAACCAGAAAGTATTTGAGATGGGGAGCTGGTATCAAGGGGAATTATTCAGTGTACAGATCAATGAGGTTAATGCAGCCCTCCTCCCTTCACTCCCCAGAAAACTCCTGACCTCTGGACACCGGGATTTTCCCATCAAGTTTTGGCCCTATTTGCTGGATCATCCACTCGCAGAACTCTTTGTCAAATAAAATGGCAGGAGCATCTCCCTGTTCCTGAGCCCAGTCAGCAAATTCGGGCAGGCGAGGCACCCGAGTGTTCTGGGAAAAGTCGAAGAAGAAAAGCTGGTCCTCGTCGTAGGCCTCAGAGAGTCCCACACTGGGACTCCCATCCTGGCAGTACACTGTGTGCAGGAATGTGTGGTTTTGCAGGTCATCTGGCCACATTGGAGTAGGAGCTGCAAAGGACACAGGGTGAGGTTCAGGGAGGTGGGAGCCTTCTCCTCCAACTTAAAAAACAGCAAGGTGGGGCTAGGCGCAGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGTGGATCATGAGGTCAGGAGTTTGAGACCAGCCTGGCCAGCATGGTGAAACTCCGTCTCTACTAAAAATACAAAAAAGTAGCTGGGCATGTTGGCATGCGCCTGTAGCTACTCGGGAGGCTGAGGGAGGAGAATTGCTTGAACCAGGGAGGCAGAGGTTGCCGGGAGCTAAGATTAAGCCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCTGTCTCAAAACAAAACAACAAAAACAAGCAAGGCCTGCTTAAGGAGCGTGGGCTGAGGTGAGACCCTTTCCTGTGTCTGTTATTTAGACTCCCCCTCCCAAAGGGGGTGAAGAACAAATTATGGCATCTCTCCAAGCTTCCCCTGCCTATAAAAAGGCCAGTTGGCAAAAGTAAAGAGTTCTACTTTCTAAAGTGACAGATTCAGGCCAGGCATGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAGGCAGATTGCTTGAGCCCAGGAGTTCAAGACCAACCTGGGCAACACAGCGAGACCCTGTCTCTACAAAAAATACAAAAACTTAGCCAGGTGTGGTGGCAAACACCTGTGGTCTCAGCTACTCTGGAGGCTGAGGCAGGAGGATTGCTTGTGCCTAGGAAGTTGGGGCTGCAGTGAGCCATGATTGTGCCACTGGACTCCAGCCCAGGTGACAGAATGAGCCCGTCTCAAAAAATATATATATAAAGGCCGGGCGCGGTGGCTCAAGCTTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCAAACATGATGAAACCCCATCTCTACTAAAAATACAAAAATCAGCTGGGTGTGGTGGCATGCGCCTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAGTCTCTTGAACCCCAGAGGCAGGGGTTGCAGGGAGCCGAGATCACGTCACTGCACTCTACCCTGGGTGACAGAGCGAGATGCCGTGTCAAAAAAAATAAATTAAATCAAATAAAAAATTTAAAAATGTATATATATAAAATAAAGTGACAGATTCAGAGTCACTGTTCATTGTGTGTTTGGGGGCTGCACAAAGACACCTAGCCAAAGAAGCAAGTGAAAGCCTGCATTCTGCTCACCATGCCATACATCCTGGCATAGGGCTGTATCCTCCCAAAGGGGATTCCTTTGTCTAATTCATACCAGGCCACTGTATTGACTAGAGAAGGCCATGGATGGGTTTCTCACTCTTAGAAGGGAAAGAGGAGGAATGGCTACAGCCTCCCCAAGCCATAGATGGGACTGCCTCCCACTATCCCCAGACACAAATGGTAAATTGGAAAACCTGTATCCAGACATTTCTTCAGCCACTTCATTGGCACCAAGCGTCTCTCAAAATGTCTTCTGTTCCTTAACCTACCAGGCCTCCCAAAGACAGCAATGGGAGAAGTGACCCCATAACTGCATAAAATAATCCCTCTTCTTTGAAGCTCTTGGCAGGAATCGCTCAGCCAGCAGGAAACCTTTAACCCAATACCCAGAAAAACAGACATTTGGAGGAAGAGGGATCTTCCAGATTATTCTTCCATTCTGCCCCATCCTCTACAGAGAAGGAAACTAAGACACTTTTCAAGAATCACAAGATAAGTTAATGATAGAAAGCAGAGTAGAATCTTGAGTGGAGGAGTGAAAATAACATTCACTTTGTTCAAATCCCAGCTCTACCACTTTCCAATGGTGTGAACTTGCACAAATAACTCTGAGTCTCATTTTCTTCATTTGTAAAATGGAGAGAACAATCTCCGCTTCAAGAGATTGTCTTAAATGGAACATGCAAAGCATCACTGATATCGTTTACCAACCACACATAGCAGCTGTCTTTCCCCACTCCCCTGTTGTTTCCACTGCCTCATAAGACTTCCCACCACTCACAAAGCACAGCGCTTTTCCTCACAAAGCTGAGTGGGCTCCCTAGGTTCAGGATGGAAGTAAATAGGAGTACCATCTTACCTTCAGGGACGGCCCAGGAGTGGGGTAGCAGCCACAGAAGTGGTAACATCTGTAGCAGCGCAGCTCCTTGGTTCTGTTCATGACCCATACCTTCTTGCCACACAGTAGGTAGGAGCTACCAACCCAGCCAACCCAGCTTCCCCAACTCCCTCCCCGAGAGGGTGGCCTTAGAT";
    let mut sequence = seq.chars().collect::<Vec<char>>();
    sequence.insert(0, '$');

    let graph_struct = graph::read_graph(
        &"tests/DMA-3108.fa.353ea42.34ee7b1.1576367.smooth.fix.gfa",
        false,
    );

    let score_matrix = score_matrix::create_score_matrix_match_mis(2, -4);
    let hofp = HashMap::new();
    let mut scores_f32 = HashMap::new();
    for (k, v) in score_matrix.iter() {
        scores_f32.insert(*k, *v as f32);
    }
    let r_values = utils::set_r_values(
        &graph_struct.nwp,
        &graph_struct.pred_hash,
        graph_struct.lnz.len(),
    );
    let mut group = c.benchmark_group("Global abPOA");
    for i in [1, 2, 3] {
        unsafe {
            group.bench_with_input(BenchmarkId::new("Simd", i), &i, |b, _i| {
                b.iter(|| {
                    global_abpoa::exec_simd(
                        black_box(&sequence),
                        black_box(("name", 0)),
                        black_box(&graph_struct),
                        black_box(&scores_f32),
                        300,
                        false,
                        &hofp,
                        &r_values,
                    )
                })
            });

            group.bench_with_input(BenchmarkId::new("No Simd", i), &i, |b, _i| {
                b.iter(|| {
                    global_abpoa::exec(
                        black_box(&sequence),
                        black_box(("name", 0)),
                        black_box(&graph_struct),
                        black_box(&score_matrix),
                        300,
                        false,
                        &hofp,
                    )
                })
            });
        }
    }
}

criterion_group!(benches, bench_global_simd_no_simd, bench_local_simd_no_simd,);
criterion_main!(benches);
