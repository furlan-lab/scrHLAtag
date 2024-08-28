#![allow(non_snake_case)]
// #![allow(unused_imports)]
#[macro_use] 

extern crate arrayref;

use simple_log::LogConfigBuilder;
use simple_log::info;
// use std::path::Path;

pub mod vc;
pub mod scrhlatag;
use crate::scrhlatag::*;
fn main() {

    let params = load_params().unwrap();
    let in_params = load_params().unwrap();
    let checked_params = check_params(in_params).unwrap();
    if checked_params.verbose {
        eprintln!("Writing log file: '{}'", &checked_params.output_path.join("scrHLAtag.log").to_str().unwrap());
    }

    let config = LogConfigBuilder::builder()
        .path(checked_params.output_path.join("scrHLAtag.log").to_str().unwrap())
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("info")
        .output_file()
        .build();
    let _ = simple_log::new(config);

    let runs = create_runs(&checked_params);

    let alleles_query = read_allelesfile(&checked_params);

    
    info!("\t\t\t\tRunning with {} thread(s)!", &checked_params.threads);
    if checked_params.verbose {
        eprintln!("Running with {} thread(s)!", &checked_params.threads);
    }

    if checked_params.verbose {
        eprintln!("Checking command line programs!");
    }
    info!("\t\t\t\tChecking command line programs!");
    let _prog_test_res = test_progs("minimap2".to_string());
    let _prog_test_res = test_progs("samtools".to_string());


    if checked_params.verbose {
        eprintln!("\nMaking fastq from BAM: '{}'", &checked_params.bam);
    }
    info!("\t\t\t\tMaking fastq from BAM: '{}'", &checked_params.bam);
    let _fq = make_fastq(&checked_params);

            
    for run in runs {
            if run.params.verbose {
                eprintln!("Matching HLA alleles with known {} reference!", &run.level.descriptor);
            }
            info!("\t\t\t\tMatching HLA alleles with known references!");

            if run.params.verbose {
                eprintln!("Making partial reference: '{}'", &run.level.mini_fasta.to_str().unwrap());
            }
            info!("\t\t\t\tMaking partial reference: '{}'", &run.level.mini_fasta.to_str().unwrap());

            let _fasta_names = make_partial_reference(&alleles_query, &run);
            // align sort and count have their own loggers
            let _al = align(&run);
            let _so = sort(&run);
            let (counts, molecule_info) = count(&run);

            let _wc = write_counts(counts, &run);
            let _wm = write_molecules(molecule_info, &run);

            if run.params.verbose {
                eprintln!("Cleaning up!\n");
            }  
            info!("\t\t\t\tCleaning up!");
            let _cu = cleanup(&run.level.out_unsorted_sam, true);
            let _cu = cleanup(&run.level.out_sorted_sam, true);

    }
    let _cu = cleanup(&checked_params.output_path.join("fastq.fq.gz"), true);
    if params.verbose {
        eprintln!("Done!!!");
    }
    info!("\t\t\t\tDone!!!");

}
