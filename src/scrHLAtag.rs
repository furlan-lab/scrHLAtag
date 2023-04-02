#![allow(non_snake_case)]


/**

target/release/scrHLAtag -b data/AML_403_34_HLA.dedup.bam -a data/testhla.tsv

**/
// scrHLA typing, alignment, : single cell rna-based HLA typing and alignment
// this program takes a pacbio bam and performs alignment and prediction
// single cell hla typing, alignment, and grouping - scrHLAtag


// #[allow(non_snake_case)]
// #[allow(unused_must_use)]
// #[allow(dead_code)]

extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
// extern crate itertools;
// extern crate fasta;
extern crate seq_io;
extern crate minimap2;


// use std::io
// use std::fs;
use clap::{App, load_yaml};
use std::str;
use std::error::Error;
use serde::Deserialize;
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{BufReader};
// use itertools::Itertools;
// use flate2::GzBuilder;
// use flate2::Compression;
use simple_log::info;
use std::path::Path;
// use fastq::parse_path;
// use fastq::each_zipped;
use simple_log::LogConfigBuilder;
// use fastq::RefRecord;
// use crate::fastq::Record;
// use std::ffi::OsStr;
// use flate2::{read};
use std::process::{Command, Stdio};
// use std::ffi::OsStr;
// use fasta::read::FastaReader;
use std::collections::HashMap;
use seq_io::fasta::{Reader,Record};
 use minimap2::Aligner;



// #[derive(Clone)]

#[allow(dead_code)]
#[derive(Deserialize)]
struct HLAalleles {
    allele: String,
    index: Option<usize>,
}

// // Implement `Display` for `Variant`.
// impl fmt::Display for Variant {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         // Use `self.number` to refer to each positional data point.
//         write!(f, "seq: {} start: {} ref_nt: {} query_nt: {} name: {}", self.seq, self.start, self.ref_nt, self.query_nt, self.name)
//     }
// }

#[allow(dead_code)]
#[derive(Debug)]
struct Params {
    bam: String,
    genome: String,
    threads: usize,
    alleles_file: String,
    output: String,
    verbose: bool,
}


fn load_params() -> Params {
    let yaml = load_yaml!("params_scrHLAtag.yml");
    let params = App::from_yaml(yaml).get_matches();
        let bam = params.value_of("bam").unwrap();
        let alleles_file = params.value_of("alleles").unwrap();
        let output = params.value_of("output").unwrap_or("out.bam");
        let genome = params.value_of("genome").unwrap_or("genome.fasta");
        let threads = params.value_of("threads").unwrap_or("1");
        let threads = threads.to_string().parse::<usize>().unwrap();
        let mut verbose = true;
        if params.is_present("verbose") {
                verbose = false
        };

        Params{
                bam: bam.to_string(),
                genome: genome.to_string(),
                threads: threads as usize,
                alleles_file: alleles_file.to_string(),
                output: output.to_string(),
                verbose: verbose,
        }
}




fn read_allelesfile(params: &Params) -> Vec<HLAalleles> {
    eprintln!("Opening alleles file: {}\n", &params.alleles_file.to_string());
    let file = File::open(&params.alleles_file.to_string()).unwrap();
    let reader = BufReader::new(file);
    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(reader);
    let mut csvdata: Vec<String> = Vec::new();
    for result in rdr.deserialize() {
        csvdata.push(result.unwrap());
    }
    csvdata.sort();
    csvdata.dedup();
    let mut allele_vec: Vec<HLAalleles>  = Vec::new();
    for data in csvdata {
        let record: HLAalleles = HLAalleles{allele: data, index: None};
        allele_vec.push(record);
    }
    allele_vec
}

fn parse_fasta_header (input: String)-> String {
    let chunks: Vec<_> = input.split(" ").collect();
    let substring1 = chunks[0].replace(">", "").replace("|", "*");
    return substring1
}

#[allow(unused_must_use)]
fn make_partial_reference (ref_file: String, alleles_query: Vec<HLAalleles>) -> Result<(), Box<dyn Error>>{
    let mut header: HashMap<String, usize> = HashMap::new();
    let infile = Path::new(&ref_file);
    let mut i: usize = 0;
    let mut reader = Reader::from_path(&infile).unwrap();
    while let Some(record) = reader.next() {
        i = i+1;
        let record = record.unwrap();
        header.insert(parse_fasta_header(record.id().unwrap().to_string()), i);
        // println!("{}", record.id().unwrap());
    }
    let mut alleles_confirmed: Vec<HLAalleles> = Vec::new();
    let mut simpleindices = Vec::new();
    for data in alleles_query {
        let matched = header.get(&data.allele);
        if matched.is_some(){
            eprintln!("Found: {}", &data.allele);
            let index = *matched.unwrap();
            alleles_confirmed.push(HLAalleles{allele: data.allele, index: Some(index)});
            simpleindices.push(index)
        } else{ 
            eprintln!("Could not find: {}", &data.allele)
        };
    }
    i = 0;
    let align_fasta = "data/align.fasta".to_string();
    let align_fasta = Path::new(&align_fasta);
    let mut reader = Reader::from_path(&infile).unwrap();
    // let mut output = io::stdout();
    let mut output = File::create(align_fasta)?;
    while let Some(record) = reader.next() {
        i = i+1;
        let record = record.unwrap();
        if simpleindices.contains(&i){
            record.write_wrap(&mut output, 80);
        }
    }
    Ok(())
}



fn main() {
    let config = LogConfigBuilder::builder()
        .path("./scrHLAtag.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();
    let _ = simple_log::new(config);
    info!("starting!");

    let params = load_params();
        if params.verbose {
        eprintln!("\n\n\n\nParsing Parameters!\n");
        eprintln!("\n\n\n\nRunning with {} threads!\n", &params.threads);
    }
    if params.verbose {
        eprintln!("\n\n\n\nChecking programs and parsing alleles_file!\n");
    }
    // let _prog_test_res = test_progs();
    let alleles_query = read_allelesfile(&params);
    
    if params.verbose {
        eprintln!("\n\n\n\nMatching HLA alleles with known references!\n");
    }

    let ref_file = "data/hla_mRNA.fasta".to_string();
    if params.verbose {
        eprintln!("\n\n\n\nChecking programs, parsing reference file, and making partial reference: {}\n", &ref_file);
    }
    let _mpr = make_partial_reference(ref_file, alleles_query);
    
    let _ar = align(&params);

}

// minimap2 --MD -a $fa -t 8 mutcaller_R1.fq.gz -o Aligned.mm2.sam
// samtools sort -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sam
// samtools view -b -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sorted.bam
// samtools index -@ 8 Aligned.mm2.ssorted.bam

// fn test_progs () -> Result<(), Box<dyn Error>>{
//     let _output = Command::new("minimap2")
//                     .arg("-h")
//                     .stderr(Stdio::piped())
//                     .stdout(Stdio::piped())
//                      .output()
//                      .expect("\n\n*******Failed to execute minimap2*******\n\n");
//     let _output = Command::new("samtools")
//                     .arg("-h")
//                     .stderr(Stdio::piped())
//                     .stdout(Stdio::piped())
//                      .output()
//                      .expect("\n\n*******Failed to execute samtools*******\n\n");
//     // eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     Ok(())
// }


fn align (params: &Params)-> Result<(), Box<dyn Error>> {
    let aligner = Aligner::builder()
        .map_hifi()
        .with_threads(1)
        .with_cigar()
        .with_index("data/align.fasta", None)
        .expect("Unable to build index");

    //  test with baked seq
    // let seq: Vec<u8> = b"TTTCTTATATGGGGAGAATCTCCTCAGACGCCGAGATGCGGGTCACGGCACCCCGAACCGTCCTCCTGCTGCTCTGGGGGGCAGTGGCCCTGACCGAGACCTGGGCCGGCTCCCACTCCATGAGGTATTTCTACACCGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGGGGAGACACGGAACATGAAGGCCTCCGCGCAGACTTACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGGTCTCACATCATCCAGGTGATGTATGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATAACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGTGGCGGAGCAGCGGAGAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCAAAGACACATGTGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACACCGAGCTTGTGGAGACCAGACCAGCAGGAGATAGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACATGCCATGTACAGCATGAGGGGCTGCCAAAGCCCCTCACCCTGAGATGGGAGCCATCTTCCCAATCCACCGTCCCCATCGTGGGCATTGTTGCTGGCCTGGCTGTCCTAGCAGTTGTGGTCATCGGAGCTGTGGTCGCTGCTGTGATGTGTAGGAGGAAGAGCTCAGGTGGAAAAGGAGGGAGCTACTCTCAGGCTGCGTGCAGCGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGAAAAGCCTGAGACAGCTGTCTTGTGAGGGACTGAGATGCAGGATTTCTTCACGCCTCCCCTTTGTGACTTCAAGAGCCTCTGGCATCTCTTTCTGCAAAGGCACCTGAATGTGTCTGCGTCCCTGTTAGCCTAATGTGAGGAGGTGGAGAGACAGCCCAACCTTGTGTCCACTGTGACCCCTGTTCCCATGCTGACCTGTGTTTCCTCCCCAGTCATCTTTCTTGTTCCAGAGAGGTGGGGCTGGATGTCTCCATCTCTGTCTCAACTTTATGTGCACTGAGCTGCAACTTCTTACTTCCCTGCTGAAAATAAGAATCTGAATATCAATTTGTTTTCTCAAATATTTGCTATGAGAGGTTGATGGATTAATTAAATAAGTCAATTCCTGGAATTTGAGAGAGCAAATAAAGACCTGAGAACCTTCCAG".to_vec();
    // eprintln!("{:?}", &seq);
    // let alignment = aligner
    //     .map(&seq, false, false, None, None)
    //     .expect("Unable to align");
    //     eprintln!("{:?}", alignment);

    let reader = bam::BamReader::from_path(&params.bam, 0).unwrap();

    for record in reader {
        let seq = record.unwrap().sequence().to_vec();
        let alignment = aligner
            .map(&seq, false, false, None, None)
            .expect("Unable to align");
        eprintln!("{:?}", alignment);
    }
    Ok(())
}

// fn align (params: &Params)-> Result<(), Box<dyn Error>> {
//     let keep = true;
//     // let mm_cmd = "/Users/sfurlan/.local/bin/minimap2";
//     // let mm_args = "-h";
//     // let mm_args_pre = format!("-a {} -t {} mutcaller_R1.fq.gz | samtools sort -@ {} | samtools view -@ {} -o Aligned.mm2.bam", params.genome.to_string(), params.threads.to_string(), params.threads.to_string(), params.threads.to_string());
//     // let mm_args_vec = mm_args_pre.split(" ").collect::<Vec<&str>>();
//     // let mm_args = OsString::new();
//     // let mm_args = OsString::from(mm_args_pre);
//     // let mm_cmd = "/Users/sfurlan/.local/bin/minimap2";
//     // let mm_cmd = "ls";
//     // let st_cmd = format!("sammtools index -@ {} Aligned.mm2.bam", params.threads.to_string());
//     // let output = Command::new("echo")
//     //                  .arg("Hello world")
//     //                  .output()
//     //                  .expect("Failed to execute command");

//     // eprintln!("{:?}", &mm_args);
//     // let output = Command::new(mm_cmd)
//     //                 .arg(mm_args)
//     //                  .output()
//     //                  .expect("Failed to execute minimap2");
//     // let output = Command::new("/Users/sfurlan/.local/bin/minimap2")
//     //                 .arg("-h")
//     //                  .output()
//     //                  .expect("Failed to execute minimap2");
//     eprintln!("{}", "Aligning reads using minimap2");
//     let output = Command::new("minimap2")
//                     .arg("--MD")
//                     .arg("-a")
//                     .arg("data/align.fasta".to_string())
//                     .arg("-t")
//                     .arg(params.threads.to_string())
//                     .arg(params.bam.to_string())
//                     .arg("-o")
//                     .arg(params.output.to_string())
//                     .stderr(Stdio::piped())
//                     .stdout(Stdio::piped())
//                      .output()
//                      .expect("\n\n*******Failed to execute minimap2*******\n\n");
//     eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     eprintln!("{}", "Minimap2 complete; Running samtools sort");
//     let output = Command::new("samtools")
//                     .arg("sort")
//                     .arg("-@")
//                     .arg(params.threads.to_string())
//                     .arg("-o")
//                     .arg("Aligned.mm2.sorted.sam")
//                     .arg(params.output.to_string())
//                     .stderr(Stdio::piped())
//                     .stdout(Stdio::piped())
//                     .output()
//                      .expect("\n\n*******Failed to execute samtools view*******\n\n");
//     eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     eprintln!("{}", "Samtools sort complete; Running samtools view");
//     let output = Command::new("samtools")
//                     .arg("view")
//                     .arg("-b")
//                     .arg("-@")
//                     .arg(params.threads.to_string())
//                     .arg("-o")
//                     .arg("Aligned.mm2.sorted.bam")
//                     .arg("Aligned.mm2.sorted.sam")
//                     .stderr(Stdio::piped())
//                     .stdout(Stdio::piped())
//                     .output()
//                      .expect("\n\n*******Failed to execute samtools sort*******\n\n");
//     eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     eprintln!("{}", "Samtools view complete; Running samtools index");
//     let output = Command::new("samtools")
//                     .arg("index")
//                     .arg("-@")
//                     .arg(params.threads.to_string())
//                     .arg("Aligned.mm2.sorted.bam")
//                     .stdout(Stdio::piped())
//                     .stderr(Stdio::piped())
//                     .output()
//                      .expect("\n\n*******Failed to execute samtools index*******\n\n");
//     eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     if keep {
//         fs::remove_file("Aligned.mm2.sorted.sam")?;
//         fs::remove_file("Aligned.mm2.sam")?;
//         fs::remove_file("mutcaller_R1.fq.gz")?;
//     }
//     Ok(())
// }



// fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, fname: String) -> Result<(), Box<dyn Error>> {
//         let f = File::create(fname)?;
//         let mut gz = GzBuilder::new()
//                         .filename("counts_mm.txt.gz")
//                         .write(f, Compression::default());
//         for result in count_vec {
//             for line in result {
//                 gz.write_all(&line)?;
//             }
//         }
//         gz.finish()?;
//         Ok(())
// }




// fn remove_whitespace(s: &mut String) {
//     s.retain(|c| !c.is_whitespace());
// }



// fn fastq(params: &Params) -> Result<(), Box<dyn Error>>{
//     let split = "|BARCODE=".to_string();
//     let outfastq = "mutcaller_R1.fq.gz".to_string();
//     let mut cbvec = lines_from_file(&params.bcs);
//     cbvec.sort_unstable();
//     let _zip = true;
//     let mut total_count: usize = 0;
//     let mut nfound_count: usize = 0;
//     let mut mmcb_count: usize = 0;
//     let split_at = &params.umi_len + &params.cb_len;
//     // let sep: Vec::<u8> = params.name_sep.as_bytes().to_vec();

//     let fastq1 = &params.fastq1;
//     let fastq2 = &params.fastq2;
//     let _counts = (0u64, 0u64);
//     let path = Path::new(&outfastq);
//     let _file = match File::create(&path) {
//         Err(_why) => panic!("couldn't open {}", path.display()),
//         Ok(file) => file,
//     };
//     // let mut writer = io::stdout();
//     let f = File::create(&outfastq)?;
//     let mut writer = GzBuilder::new()
//                             .filename(outfastq)
//                             .write(f, Compression::default());
//     parse_path(Some(fastq1), |parser1| {
//         parse_path(Some(fastq2), |parser2| {
//             each_zipped(parser1, parser2, |rec1, rec2| {
//                 if rec1.is_some() & rec2.is_some(){
//                     let r1 = &rec1.unwrap();
//                     let r2 = &rec2.unwrap();
//                     if r1.seq().contains(&b"N"[0]) | r2.seq().contains(&b"N"[0]){
//                         nfound_count += 1;
//                         total_count +=1;
//                     }else{
//                         total_count +=1;
//                         let (barcode, _seq) = &r1.seq().split_at(split_at.into());
//                         let (cb, _seq) = barcode.split_at(params.cb_len as usize);
//                         match cbvec.binary_search(&std::str::from_utf8(cb).unwrap().to_string()) {
//                             Ok(_u) => {
//                                 let mut readout = RefRecord::to_owned_record(&r2);
//                                 let _some_x = vec![b" "];
//                                 let mut new_header = std::str::from_utf8(&readout.head()).unwrap().to_string();
//                                 remove_whitespace(&mut new_header);
//                                 let _ = new_header.push_str(&split);
//                                 let _ = new_header.push_str(&std::str::from_utf8(&barcode).unwrap().to_string());
//                                 readout.head = new_header.as_bytes().to_vec();

//                                 let _ = readout.write(&mut writer);
//                             }
//                             Err(_e) => {
//                                 mmcb_count +=1;
//                             }
//                         }
//                     }
//                 }
//                 (true, true)
//             })
//             .expect("Invalid record.");
//         })
//         .expect("Unknown format for file 2.");
//     })
//     .expect("Unknown format for file 1.");
//     eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist\n", total_count, nfound_count, mmcb_count);
//     Ok(())
// }



// fn process_variant(ref_id: u32, start: u32)->bam::Region{
//     let region = bam::Region::new(ref_id,start - 1,start - 1);
//     return region;
// }


// fn count_variants_mm2(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
//     eprintln!("Processing using cb and umi in header");
//     let split = "|BARCODE=".to_string();
//     let ibam = "Aligned.mm2.sorted.bam";
//     let mut total: usize = 0;
//     let mut err: usize = 0;
//     let seqname = variant.seq;
//     let start = variant.start.parse::<u32>().unwrap();
//     let vname = variant.name;
//     let mut reader = bam::IndexedReader::build()
//         .additional_threads(*&params.threads as u16)
//         .from_path(ibam).unwrap();
//     let mut seqnames = Vec::new();
//     let mut _result = "";
//     let query_nt = variant.query_nt as char;
//     let header = reader.header().clone();
//     let hdata = header.reference_names();
//     for seq in hdata {
//         seqnames.push(seq)
//     }
//     let mut data = Vec::new();
//     let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
//     let region = process_variant(ref_id as u32, start);
//     for record in reader.fetch_by(&&region, |record| record.mapq() >= 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
//         total+=1;
//         let readheader = match str::from_utf8(record.as_ref().unwrap().name()) {
//             Ok(v) => v,
//             Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
//         };
//         let cbumi = match readheader.split(&split).nth(1){
//             Some(v) => v.to_string(),
//             None => {
//                 err+=1;
//                 continue
//             },
//         };
//         if cbumi.len() <= params.cb_len+1 {
//             err+=1;
//             continue
//         }
//         let (cb, umi) = cbumi.split_at((params.cb_len+1).into());
//         for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
//             if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
//                 if region.start() == ref_pos {
//                     if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
//                         if ref_nt as char == record_nt as char {
//                             _result = "ref";
//                         } else if record_nt as char == query_nt{
//                             _result = "query";
//                         } else {
//                             _result = "other";
//                         }
//                             data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, _result))
//                         }
//                     } else {
//                         continue
//                     }
//             } else {
//                 continue
//             }        }
//     }
//     eprintln!("Found {} reads spanning this variant!\n\tNumbers of errors: {}", total, err);
//     data.sort();
//     let mut out_vec = Vec::new();
//     let cdata = data.into_iter().dedup_with_count();
//     for (count, record) in cdata {
//        let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
//         out_vec.push(count_str.as_bytes().to_owned());
//     }
//     return out_vec;
// }


// fn count_variants_mm(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
//     eprintln!("Processing using cb and umi in BAM tags");
//     // let split = "|BARCODE=".to_string();
//     let ibam = "Aligned.mm2.bam";
//     let mut total: usize = 0;
//     let seqname = variant.seq;
//     let start = variant.start.parse::<u32>().unwrap();
//     let vname = variant.name;
//     let mut reader = bam::IndexedReader::build()
//         .additional_threads(*&params.threads as u16)
//         .from_path(ibam).unwrap();
//     let mut seqnames = Vec::new();
//     let mut cb;
//     let mut umi;
//     let mut result = "null";
//     let query_nt = variant.query_nt as char;
//     let header = reader.header().clone();
//     let hdata = header.reference_names();
//     for seq in hdata {
//         seqnames.push(seq)
//     }
//     let mut data = Vec::new();
//     let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
//     let region = process_variant(ref_id as u32, start);
//     for record in reader.fetch_by(&&region, |record| record.mapq() > 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
//         total+=1;
//         match record.as_ref().unwrap().tags().get(b"CB") {
//             Some( bam::record::tags::TagValue::String(cba, _)) => {
//                 cb = str::from_utf8(&cba).unwrap().to_string();
//             },
//             _ => panic!("Unexpected type"),
//         }
//         match record.as_ref().unwrap().tags().get(b"UB") {
//             Some( bam::record::tags::TagValue::String(uba, _)) => {
//                 umi = str::from_utf8(&uba).unwrap().to_string();
//             },
//             _ => panic!("Unexpected type"),
//         }
//         for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
//             if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
//                 if region.start() == ref_pos {

//                     if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
//                         if ref_nt as char == record_nt as char {
//                             result = "ref";
//                         } else if record_nt as char == query_nt{
//                             result = "query";
//                         } else {
//                             result = "other";
//                         }
//                             data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, result))
//                         }
//                     } else {
//                         continue
//                     }
//             } else {
//                 continue
//             }        }
//     }
//     eprintln!("Found {} reads spanning this variant!", total);
//     data.sort();
//     let mut out_vec = Vec::new();
//     let cdata = data.into_iter().dedup_with_count();
//     for (count, record) in cdata {
//         let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
//         // let count_str = record+&" ".to_owned()+&(count.to_string())+&"\n".to_owned();
//         out_vec.push(count_str.as_bytes().to_owned());
//     }
//     return out_vec;
// }


// fn count_star(params: &Params) {
//     let ibam = "Aligned.mm2.bam";
//     let split = "|BARCODE=".to_string();
//     let joiner = "_".to_string();
//     eprintln!("Counting star reads");
//     let mut total: usize = 0;
//     let mut goodreadcount: usize = 0;
//     let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
//         (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
//     } else {
//         (0 as usize, 0 as usize)
//     };

//     let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
//     let _output = std::io::BufWriter::new(io::stdout());
//     let header = reader.header().clone();
//     let data = header.reference_names();
//     let mut seqnames = Vec::new();
//     for seq in data {
//         seqnames.push(seq)
//     }
//     for record in reader {
//         total += 1;
//         let newrecord = record.unwrap();
//         let seqname = match str::from_utf8(&newrecord.name()) {
//             Ok(v) => v,
//             Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
//         };
//         let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
//         let _modified_name = seqname.replace(&split, &joiner);
//         let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
//         let mut good_read = false;
//         let cigarmatch = format!("{}M", *&params.read_len);
//         let cigar = newrecord.cigar().to_string();
//         if cigar == cigarmatch{
//             good_read = true
//         }
//         if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
//             goodreadcount += 1;
//             println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
//         }
//     }
//     eprintln!("Completed; {} total reads processed!", &total);
//     eprintln!("{} good reads counted!", &goodreadcount);
// }



// fn count_kallisto(params: &Params) {
//     eprintln!("Counting kallisto reads");
//     let mut total: usize = 0;
//     let ibam = "Aligned.mm2.bam";
//     let split = "|BARCODE=".to_string();
//     let joiner = "_".to_string();
//     let mut goodreadcount: usize = 0;
//     let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
//         (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
//     } else {
//         (0 as usize, 0 as usize)
//     };
//     let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
//     let _output = io::BufWriter::new(io::stdout());
//     let header = reader.header().clone();
//     let data = header.reference_names();
//     let mut seqnames = Vec::new();
//     for seq in data {
//         seqnames.push(seq)
//     }
//     for record in reader {
//         total += 1;
//         let newrecord = record.unwrap();
//         let seqname = match str::from_utf8(&newrecord.name()) {
//             Ok(v) => v,
//             Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
//         };
//         let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
//         let _modified_name = seqname.replace(&split, &joiner);
//         let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
//         let mut good_read = false;
//         let cigarmatch = format!("{}M", *&params.read_len);
//         let cigar = newrecord.cigar().to_string();
//         if cigar == cigarmatch{
//             good_read = true
//         }
//         if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
//             goodreadcount += 1;
//             println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
//         }
//     }
//     eprintln!("Completed; {} total alignments processed!", &total);
//     eprintln!("{} good alignments counted!", &goodreadcount);
// }

// fn lines_from_file(filename: &str) -> Vec<String> {
//     let path = Path::new(filename);
//     let file = match File::open(&path) {
//         Err(_why) => panic!("\n\n*******couldn't open {}*******\n\n", path.display()),
//         Ok(file) => file,
//     };
//     if path.extension() == Some(OsStr::new("gz")){
//         let buf = BufReader::new(read::GzDecoder::new(file));
//         buf.lines()
//             .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
//             .collect()
//     }else{
//         let buf = BufReader::new(file);
//         buf.lines()
//             .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
//             .collect()
//     }
// }

