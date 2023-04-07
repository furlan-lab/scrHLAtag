#![allow(non_snake_case)]

pub mod scrhlatag;
use crate::scrhlatag::*;
fn main() {
    // let config = LogConfigBuilder::builder()
    //     .path("./scrHLAtag.log")
    //     .size(1 * 100)
    //     .roll_count(10)
    //     .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
    //     .level("debug")
    //     .output_file()
    //     .build();
    // let _ = simple_log::new(config);
    // info!("starting!");

    let params = load_params().unwrap();
    if params.verbose {
        eprintln!("Parsing Parameters!");
        eprintln!("Running with {} thread(s)!", &params.threads);
    }

    if params.verbose {
        eprintln!("Checking command line programs!");
    }
    let _prog_test_res = test_progs("minimap2".to_string());
    let _prog_test_res = test_progs("samtools".to_string());
    let alleles_query = read_allelesfile(&params);
    
    if params.verbose {
        eprintln!("Matching HLA alleles with known references!");
    }

    let align_fasta_path = &params.output.join("align.fa");
    let align_fasta = align_fasta_path.to_str().unwrap();
    if params.verbose {
        eprintln!("Making partial reference: '{}'", &align_fasta);
    }
    let _fasta_names = make_partial_reference(alleles_query, &params);
    if params.verbose {
        eprintln!("\nMaking fastq from BAM: '{}'", &params.bam);
    }
    let _fq = make_fastq(&params);
    let _al = align(&params);
    let _so = sort(&params);
    let counts = count(&params);
    // let _ar = align_and_count(&params, fasta_names.unwrap());  
    if params.verbose {
        eprintln!("Cleaning up!\n");
    }  
    let _cu = cleanup(&params.output.join("Aligned_mm2.sam"));
    let _cu = cleanup(&params.output.join("Aligned_mm2_sorted.sam"));
    let _cu = cleanup(&params.output.join("fastq.fq.gz"));
    let _wc = write_counts(counts, &params);
    // let _so = sort(&params);
    if params.verbose {
        eprintln!("Done!!!");
    }

}
