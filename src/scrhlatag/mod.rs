/**

~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/test.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out -l transcriptome -s
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/full_dedup.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/full_corrected_sorted.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/mini_corrected_sorted.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out

**/
// scrHLA typing, alignment, : single cell rna-based HLA typing and alignment
// this program takes a pacbio bam and performs alignment and prediction
// single cell hla typing, alignment, and grouping - scrHLAtag


extern crate csv;
extern crate clap;
extern crate serde;
extern crate bam;
extern crate fastq;
extern crate kseq;


use std::{env, io::Write, str, fs, fs::File, error::Error, path::{Path, PathBuf}, collections::HashMap, process::{Command, Stdio }};
use std::io::{ BufReader, BufWriter};
use clap::{App, load_yaml};
use serde::Deserialize;
use csv::ReaderBuilder;
use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
use bam::{BamReader, record::tags::TagValue};
use flate2::{Compression, GzBuilder};
use fastq::{OwnedRecord, Record};
use itertools::Itertools;
use kseq::parse_path;
use simple_log::LogConfigBuilder;
use simple_log::{info, warn, error};
use toml::from_str;





#[derive(Deserialize)]
pub struct HLAalleles {
    allele: String
}

#[derive(Clone)]
pub struct InputParams {
    pub bam: String,
    pub threads: usize,
    pub alleles_file: String,
    pub output_path: Box<Path>,
    pub verbose: bool,
    pub return_sequence: bool, 
    pub cb_tag: String,
    pub umi_tag: String,
    pub level: String,
}

pub struct AlignmentLevel {
    pub file_tag: String,
    pub hla_ref: PathBuf,
    pub descriptor: String,
    pub mini_fasta: PathBuf,
    pub out_unsorted_sam: PathBuf,
    pub out_sorted_sam: PathBuf,
    pub out_sorted_bam: PathBuf,
    pub counts: PathBuf,
    pub molecule: PathBuf,
}

pub struct Run {
    pub level: AlignmentLevel,
    pub params: InputParams,
}

pub fn create_runs(params: &InputParams) -> Vec<Run> {
    let package_dir_path = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut runs = Vec::new();
    if params.level == "transcriptome" {
        runs.push(Run{
            level: AlignmentLevel{
                file_tag: "mRNA".to_string(),
                hla_ref: package_dir_path.join("data/HLA_DB_3field_mRNA.fa.gz"),
                descriptor: "transcriptome".to_string(),
                mini_fasta: params.output_path.join("align_mRNA.fa"),
                out_unsorted_sam: params.output_path.join("Aligned_mm2_mRNA.sam"),
                out_sorted_sam: params.output_path.join("Aligned_mm2_sorted_mRNA.sam"),
                out_sorted_bam: params.output_path.join("Aligned_mm2_sorted_mRNA.bam"),
                counts: params.output_path.join("counts_mRNA.txt.gz"),
                molecule: params.output_path.join("molecule_info_mRNA.txt.gz"),
            },
            params: params.clone(),
        });
        return runs;
    } else {
        runs.push(Run{
            level: AlignmentLevel{
                file_tag: "mRNA".to_string(),
                hla_ref: package_dir_path.join("data/HLA_DB_3field_gene.fa.gz"),
                descriptor: "genome".to_string(),
                mini_fasta: params.output_path.join("align_gene.fa"),
                out_unsorted_sam: params.output_path.join("Aligned_mm2_gene.sam"),
                out_sorted_sam: params.output_path.join("Aligned_mm2_sorted_gene.sam"),
                out_sorted_bam: params.output_path.join("Aligned_mm2_sorted_gene.bam"),
                counts: params.output_path.join("counts_gene.txt.gz"),
                molecule: params.output_path.join("molecule_info_gene.txt.gz"),
            },
            params: params.clone(),
        });
        if params.level == "both" {
            runs.push(Run{
                level: AlignmentLevel{
                    file_tag: "mRNA".to_string(),
                    hla_ref: package_dir_path.join("data/HLA_DB_3field_mRNA.fa.gz"),
                    descriptor: "transcriptome".to_string(),
                    mini_fasta: params.output_path.join("align_mRNA.fa"),
                    out_unsorted_sam: params.output_path.join("Aligned_mm2_mRNA.sam"),
                    out_sorted_sam: params.output_path.join("Aligned_mm2_sorted_mRNA.sam"),
                    out_sorted_bam: params.output_path.join("Aligned_mm2_sorted_mRNA.bam"),
                    counts: params.output_path.join("counts_mRNA.txt.gz"),
                    molecule: params.output_path.join("molecule_info_mRNA.txt.gz"),
                },
                params: params.clone(),
            });
        }
        return runs
    }
}


pub fn load_params() -> Result<InputParams, Box<dyn Error>> {
    let yaml = load_yaml!("../cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    // eprintln!("{:?}", params);
    let bam = params.value_of("bam").unwrap().to_string();
    let alleles_file = params.value_of("alleles_file").unwrap().to_string();
    let level = params.value_of("align_level").unwrap_or("both").to_string();
    let cb = params.value_of("cb").unwrap_or("CB").to_string();
    let umi = params.value_of("umi").unwrap_or("XM").to_string();
    let output = params.value_of("output_folder").unwrap_or("out").to_string();
    let threads = params.value_of("threads").unwrap_or("1").to_string().parse::<usize>().unwrap();
    let mut return_sequence = false;
    if params.is_present("return_sequence") {
        return_sequence = true;
    }
    let mut verbose = false;
    if params.is_present("verbose") {
        verbose = true;
    }
    let outpath = Path::new(&output);
    Ok(InputParams{
            bam: bam,
            threads: threads as usize,
            alleles_file: alleles_file,
            output_path: outpath.into(),
            verbose: verbose,
            return_sequence: return_sequence,
            // hla_ref: hla_ref.to_string(),
            cb_tag: cb,
            umi_tag: umi,
            level: level.to_string(),
    })
}


pub fn check_params(params: InputParams) -> Result<InputParams, Box<dyn Error>>{
    let _cu = cleanup(&params.output_path.join("scrHLAtag.log"), false);
    let config = LogConfigBuilder::builder()
        .path(params.output_path.join("scrHLAtag.log").to_str().unwrap())
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("info")
        .output_file()
        .build();
    let _ = simple_log::new(config);
    let allowables = vec!["genome".to_string(), "transcrptome".into(), "both".into()];
    if !allowables.contains(&params.level){
        error!("\t\tInput 'level' parameter is not valid.  User supplied {}", params.level);
        panic!("Input 'level' parameter is not valid.  User supplied {}", params.level);
    }
    info!("\t\t\tStarting!");
    let wdpb= get_current_working_dir().unwrap();
    let wdir = wdpb.to_str().unwrap();
    info!("\t\t\tCurrent working directory: '{}'", wdir);
    if params.verbose {
        eprintln!("\n\nCurrent working directory: '{}'", wdir);
    }
    if params.output_path.is_relative(){
        let a1 = Path::new(wdir).join(&params.output_path);
        let abs_outpath = a1.to_str().unwrap();
        if params.output_path.exists() {
                info!("\t\t\tFound existing output directory: '{}'", &abs_outpath);
                warn!("\t\t\t{}", "Existing data in this folder could be lost!!!");
               if params.verbose {
                    eprintln!("Found existing output directory: '{}'", &abs_outpath);
                    eprintln!("\t{}", "Existing data in this folder could be lost!!!");
                }
        } else {
            info!("\t\tCreating output directory: '{}'", &abs_outpath);
            if params.verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(&params.output_path)?;
        }
    } else {
        let abs_outpath = &params.output_path.to_str().unwrap();
        if params.output_path.exists() {
            info!("\t\tFound existing output directory: '{}'", &abs_outpath);
            warn!("\t\t{}", "Existing data in this folder could be lost!!!");
            if params.verbose {
                eprintln!("Found existing output directory: '{}'", &abs_outpath);
                eprintln!("\t{}", "Existing data in this folder could be lost!!!");
            }
        } else {
            info!("\t\tCreating output directory: '{}'", &abs_outpath);
            if params.verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(&params.output_path)?;
        }
    }
    Ok(InputParams{
            bam: params.bam,
            threads: params.threads,
            alleles_file: params.alleles_file,
            output_path: params.output_path,
            verbose: params.verbose,
            return_sequence: params.return_sequence,
            // hla_ref: hla_ref.to_string(),
            cb_tag: params.cb_tag,
            umi_tag: params.umi_tag,
            level: params.level,
    })
}


pub fn read_allelesfile(params: &InputParams) -> Vec<HLAalleles> {
    info!("\t\tOpening alleles file: '{}'", &params.alleles_file.to_string());
    if params.verbose {
        eprintln!("Opening alleles file: '{}'", &params.alleles_file.to_string());
    }
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
        let record: HLAalleles = HLAalleles{allele: data};
        allele_vec.push(record);
    }
    allele_vec
}


/** "make_partial_reference" and its companion function "parse_fasta_header" take a fasta file 
 * and keep only those entries that match a list of HLA alleles passed as an allele_query. In doing so, 
 * it also parses the name to remove fasta friendly characters, replacing them with 
 * characters standard in HLA designations to match the characters in alleles_query.
 * 
 * It then writes a new fasta containing only the desired alleles in alleles_query.  It returns the names
 * of the entries for later writing into a bam header.
**/
pub fn make_partial_reference (alleles_query: &Vec<HLAalleles>, run: &Run) -> Result<Vec<String>, Box<dyn Error>>{
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut i: usize = 0;
    info!("\t\tRead HLA reference file: '{}'", &run.level.hla_ref.to_str().unwrap());
    if run.params.verbose {
        eprintln!("Read HLA reference file: '{}'", &run.level.hla_ref.to_str().unwrap());
    }

    let mut records = parse_path(&run.level.hla_ref).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        i = i+1;
        names.insert(parse_fasta_header(record.head().to_string()), i);
    }

    let mut simpleindices = Vec::new();
    for data in alleles_query {
        let matched = names.get(&data.allele);
        if matched.is_some(){
            info!("\t\t\tFound: {}", &data.allele);
            if run.params.verbose {
                eprintln!("\tFound: {}", &data.allele);
            }
            simpleindices.push(*matched.unwrap())
        } else{ 
            warn!("\t\tCould not find: {}", &data.allele);
            if run.params.verbose {
                eprintln!("Could not find: {}", &data.allele)
            }
        };
    }
    i = 0;
    let align_fasta_path = &run.level.mini_fasta;
    let align_fasta = align_fasta_path.to_str().unwrap();
    let f = File::create(align_fasta)?;

    // write align.fasta
    let mut fasta_writer = fasta::writer::Builder::default().build_with_writer(BufWriter::new(f));
    let mut names_out = Vec::new();
    let mut records = parse_path(&run.level.hla_ref).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        i = i+1;
        if simpleindices.contains(&i){
            names_out.push(record.head().to_string());
            let definition = Definition::new(record.head().to_string(), None);
            let seq = record.seq().to_string().into_bytes();
            let sequence = Sequence::from(seq);
            let frecord = noodles_fasta::record::Record::new(definition, sequence);
            let _fw = fasta_writer.write_record(&frecord);
        }
    }
    Ok(names_out.to_owned())
}


fn parse_fasta_header (input: String)-> String {
    let chunks: Vec<_> = input.split(" ").collect();
    let substring1 = chunks[0].replace(">", "").replace("|", "*");
    return substring1
}


pub fn test_progs (software: String) -> Result<(), Box<dyn Error>>{
    let _output = Command::new(software.clone())
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect(&format!("\n\n*******Failed to execute {}*******\n\n", software));
    Ok(())
}

fn calc_threads(params: &InputParams) -> u16 {
    if params.threads>2{
        (params.threads-1).try_into().unwrap()
    } else {
        0
    }

}

#[allow(unused_assignments)]
pub fn make_fastq (params: &InputParams)-> Result<(), Box<dyn Error>> {

    // declare local variables
    let split = "|BARCODE=".to_string();
    let fastq_path = &params.output_path.join("fastq.fq.gz");
    let bam_fn = Path::new(&params.bam);
    // let cb_b: [u8; 2] = clone_into_array(&params.cb_tag.as_bytes().to_vec()[0..1]);
    let binding = params.cb_tag.as_bytes().to_vec();
    let cb_b = pop2(&binding);
    let binding = params.umi_tag.as_bytes().to_vec();
    let umi_b = pop2(&binding);
    let binding = "nb".as_bytes().to_vec();
    let nb_b = pop2(&binding);
    let mut nb_present = false;
    let mut nb = 0;
    let mut new_readname: String = Default::default();

    // set counters
    let mut total_count: usize = 0;
    let mut err_count: usize = 0;

    // make bam reader
    let bam_reader = BamReader::from_path(bam_fn, calc_threads(&params)).unwrap();
    
    // make fastq writer;
    let _file = match File::create(&fastq_path) {
        Err(_why) => panic!("couldn't open {}", fastq_path.display()),
        Ok(file) => file,
    };
    let f = File::create(&fastq_path)?;
    let mut writer = GzBuilder::new()
                            .filename(fastq_path.to_str().unwrap())
                            .write(f, Compression::default());

    // read bam
    for record in bam_reader{
        let mut cb: String = Default::default();
        let mut umi: String = Default::default();
        total_count+=1;
        let rec = record.as_ref().unwrap();

        // get name
        let old_readname = match str::from_utf8(rec.name()) {
            Ok(value) => {
                remove_whitespace(value)
            },
            Err(_e) => {
                err_count+=1;
                continue
            }
        };

        // get cb and umi
        match rec.tags().get(&cb_b) {
            Some( bam::record::tags::TagValue::String(cba, _)) => {
                cb = str::from_utf8(&cba).unwrap().to_string();
                // eprintln!("{:?}", &cb);
            },
            _ => {
                err_count+=1;
                continue
            },
        }
        match rec.tags().get(&umi_b) {
            Some( bam::record::tags::TagValue::String(uba, _)) => {
                umi = str::from_utf8(&uba).unwrap().to_string();
                // eprintln!("{:?}", &umi);
            },
            _ => {
                err_count+=1;
                continue
            },
        }
        match rec.tags().get(&nb_b) {
            Some( bam::record::tags::TagValue::Int(nba, _)) => {
                nb_present = true;
                nb = nba;
            },
            _ => {
                nb_present = false;
            },
        }
        if nb_present {
            new_readname = format!("{}{}{}_{}_nb_{}",&old_readname, &split, cb, umi, nb);
        } else {
            new_readname = format!("{}{}{}_{}",&old_readname, &split, cb, umi );
        }
        
        
        // write to fastq.gz
        let new_record: OwnedRecord = OwnedRecord{head: new_readname.as_bytes().to_vec(),
                                    seq: rec.sequence().to_vec(),
                                    sep: None,
                                    qual: vec![255; rec.sequence().len()]};

        let _nr = new_record.write(&mut writer);
    }
    info!("\t\t\tTotal reads processed: {}\tReads with errors: {}", total_count, err_count);
    if params.verbose {
        eprintln!("\tTotal reads processed: {}\n\tReads with errors: {}\n", total_count, err_count);
    }
    Ok(())
}


pub fn align (run: &Run)-> Result<(), Box<dyn Error>> {
    let align_fasta_path = &run.level.mini_fasta;
    let align_fasta = align_fasta_path.to_str().unwrap();
    if !align_fasta_path.exists() {
        error!("\t\tAlignment fasta not found at: {}", align_fasta);
        panic!("Alignment fasta not found at: {}", align_fasta);
    }
    let fastq_path = run.params.output_path.join("fastq.fq.gz");
    let fastq_file = fastq_path.to_str().unwrap();
    if !fastq_path.exists() {
        error!("\t\tAlignment fastq not found at: {}", fastq_file);
        panic!("Alignment fastq not found at: {}", fastq_file);
    }
    let sam_path = &run.level.out_unsorted_sam;
    let sam_file = sam_path.to_str().unwrap();
    info!("\t\t{}", "Aligning reads using minimap2 - Output below:");
    if run.params.verbose {
        eprintln!("{}", "Aligning reads using minimap2 - Output below:\n");
    }
    let output = Command::new("minimap2")
                    .arg("--cs=long")
                    .arg("--secondary=no")
                    .arg("-x")
                    .arg("map-hifi")
                    .arg("-Q")  // TODO: this ignores base qual.  fix this later
                    .arg("--MD")
                    .arg("-a")
                    .arg("-t")
                    .arg(run.params.threads.to_string())
                    .arg(align_fasta)
                    .arg(fastq_file)
                    .arg("-o")
                    .arg(sam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute minimap2*******\n\n");
    info!("{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    Ok(())
}



pub fn sort (run: &Run)-> Result<(), Box<dyn Error>> {
    let sam_path = &run.level.out_unsorted_sam;
    let sam_file = sam_path.to_str().unwrap();
    if !sam_path.exists() {
        error!("\t\tSAM not found at: {}", sam_file);
        panic!("SAM not found at: {}", sam_file);
    }
    let ssam_path = &run.level.out_sorted_sam;
    let ssam_file = ssam_path.to_str().unwrap();
    let bam_path = &run.level.out_sorted_bam;
    let bam_file = bam_path.to_str().unwrap();
    info!("\t\t{}", "Minimap2 complete; Running samtools sort");
    if run.params.verbose {
        eprintln!("{}", "Minimap2 complete; Running samtools sort");
    }
    let output = Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg(run.params.threads.to_string())
                    .arg("-o")
                    .arg(ssam_file)
                    .arg(sam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools view*******\n\n");
    info!("\t\t{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
        
    }
    info!("\t\t{}", "Samtools sort complete; Running samtools view");
    if run.params.verbose {
        eprintln!("{}", "Samtools sort complete; Running samtools view");
    }
    
    if !ssam_path.exists() {
        error!("\t\tSorted SAM not found at: {}", ssam_file);
        panic!("Sorted SAM not found at: {}", ssam_file);
    }
    let output = Command::new("samtools")
                    .arg("view")
                    .arg("-b")
                    .arg("-@")
                    .arg(run.params.threads.to_string())
                    .arg("-o")
                    .arg(bam_file)
                    .arg(ssam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools sort*******\n\n");
    info!("\t\t{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    info!("\t\t{}", "Samtools view complete; Running samtools index");
    if run.params.verbose {
        eprintln!("{}", "Samtools view complete; Running samtools index");
    }
    if !bam_path.exists() {
        error!("Sorted SAM not found at: {}", bam_file);
        panic!("\t\tSorted SAM not found at: {}", bam_file);
    }
    let output = Command::new("samtools")
                    .arg("index")
                    .arg("-@")
                    .arg(run.params.threads.to_string())
                    .arg(bam_file)
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools index*******\n\n");
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    Ok(())
}

pub fn count(run: &Run) -> (Vec<Vec<u8>>, Vec<String>){
    info!("\t\tCounting reads:");
    if run.params.verbose {
        eprintln!("Counting reads:");
    }
    // declare some stuff
    let split = "|BARCODE=".to_string();
    let bam_path = &run.level.out_sorted_bam;
    let bam_file = bam_path.to_str().unwrap();
    let mut total_count: usize = 0;
    let mut err_count: usize = 0;
    let mut unmapped_count: usize = 0;
    let mut mapped_count: usize = 0;
    let mut bc_perfect_count = 0;
    let mut bc_imperfect_count = 0;
    let mut bc_notfound_count = 0;

    // bam reader and get seqnames
    let bam_reader = bam::BamReader::from_path(bam_file, calc_threads(&run.params)).unwrap();
    let mut seqnames = Vec::new();
    let mut _result = "";
    let header = bam_reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }

    // count
    let mut data = Vec::new();
    let mut molecule_data: Vec<String> = Vec::new();
    for record in bam_reader{
        total_count+=1;
        let readname = match str::from_utf8(record.as_ref().unwrap().name()) {
            Ok(v) => v,
            Err(_e) => {
                err_count+=1;
                continue
            },
        };
        let mut cbumi = readname.split(&split).nth(1).unwrap();
        let mut nb_present = false;
        let mut nb: i32 = 0;
        // eprintln!("{}", &cbumi);
        if cbumi.contains("_nb_") {
            nb_present = true;
            nb = cbumi.split("_nb_").nth(1).unwrap().parse::<i32>().unwrap();
            // nb = from_str::<i32>(cbumi.split("_nb_").nth(1).unwrap()).unwrap();
            cbumi  = cbumi.split("_nb_").nth(0).unwrap();
        }
        let cb = cbumi.split("_").nth(0);
        let umi = cbumi.split("_").nth(1);
        if record.as_ref().unwrap().ref_id() < 0 {
                unmapped_count+=1;
                continue;
        } else if record.as_ref().unwrap().ref_id() > seqnames.len().try_into().unwrap() {
                err_count+=1;
                continue;
        } else {
                mapped_count+=1;
                let rec = record.as_ref().unwrap();
                let index = rec.ref_id() as usize;
                let name = rec.name();
                // eprintln!("{} {} {}", &cb.unwrap(), &umi.unwrap(), seqnames[index]);
                data.push(format!("{} {} {}", &cb.unwrap(), &umi.unwrap(), seqnames[index]));

                let nm_tag = match rec.tags().get(b"NM") {
                    Some(TagValue::Int(value, _)) => value,
                    _ => {
                            warn!("NM tag not returned correctly for read {:?}", str::from_utf8(rec.name()).unwrap());
                            -1
                        },
                    };
                let as_tag = match rec.tags().get(b"AS") {
                    Some(TagValue::Int(value, _)) => value,
                    _ => {
                            warn!("AS tag not returned correctly for read {:?}", rec.name());
                            0
                        },
                };
                let s1_tag = match rec.tags().get(b"s1") {
                    Some(TagValue::Int(value, _)) => value,
                    _ => {
                            warn!("s1 tag not returned correctly for read {:?}", rec.name());
                            0
                        },
                };
                let de_tag = match rec.tags().get(b"de") {
                    Some(TagValue::Float(value)) => value,
                    _ => {
                            warn!("de tag not returned correctly for read {:?}", rec.name());
                            0.0
                        },
                };

                let _de_tag = match rec.tags().get(b"nb") {
                    Some(TagValue::Int(value, _)) =>  value,  
                    _ => {  
                            // let e = rec.tags().get(b"nb");
                            // eprintln!("{:?}", e);
                            warn!("nb tag not returned correctly for read {:?}", rec.name());
                            -1
                        },
                };

                if nb_present {
                    if nb == 0 {
                        bc_perfect_count+=1
                    } else {
                        if nb < 0 {
                            bc_notfound_count+=1
                        }
                        if nb > 0 {
                            bc_imperfect_count+=1
                        }
                    }
                    // columns cb, nb, umi, seqname, query_len, start, mapq, cigar, NM, AS, chaining_score, de (per base sequence divergence), seq (optional)
                    if run.params.return_sequence {
                        molecule_data.push(format!("{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n", String::from_utf8_lossy(name), &cb.unwrap(), nb, &umi.unwrap(), seqnames[index], rec.query_len(), rec.start(), rec.mapq(), rec.cigar(), nm_tag, as_tag, s1_tag, de_tag, String::from_utf8_lossy(&rec.sequence().to_vec())));
                    } else {
                        molecule_data.push(format!("{} {} {} {} {} {} {} {} {} {} {} {}\n", String::from_utf8_lossy(name), &cb.unwrap(), nb, &umi.unwrap(), seqnames[index], rec.start(), rec.mapq(), rec.cigar(), nm_tag, as_tag, s1_tag, de_tag));
                    }
                } else{
                    // columns cb, umi, seqname, query_len, start, mapq, cigar, NM, AS, chaining_score, de (per base sequence divergence), seq (optional)
                    if run.params.return_sequence {
                        molecule_data.push(format!("{} {} {} {} {} {} {} {} {} {} {} {} {}\n", String::from_utf8_lossy(name), &cb.unwrap(), &umi.unwrap(), seqnames[index], rec.query_len(), rec.start(), rec.mapq(), rec.cigar(), nm_tag, as_tag, s1_tag, de_tag, String::from_utf8_lossy(&rec.sequence().to_vec())));
                    } else {
                        molecule_data.push(format!("{} {} {} {} {} {} {} {} {} {} {}\n", String::from_utf8_lossy(name), &cb.unwrap(), &umi.unwrap(), seqnames[index], rec.start(), rec.mapq(), rec.cigar(), nm_tag, as_tag, s1_tag, de_tag));
                    }
                }
                
                
        }
    }
    info!("\t\t\tTotal reads processed: {}\tReads with errors: {}", total_count, err_count);
    info!("\t\t\tUnmapped reads: {}\tMapped reads: {}", unmapped_count, mapped_count);
    if bc_perfect_count + bc_notfound_count + bc_imperfect_count > 0 {
        info!("\tPerfect Barcodes: {}\n\tCorrected Barcodes: {}\n\tUnmatchable Barcodes: {}\n", bc_perfect_count, bc_imperfect_count, bc_notfound_count);
    }
    if run.params.verbose {
        eprintln!("\tTotal reads processed: {}\n\tReads with errors: {}\n\tUnmapped reads: {}\n\tMapped reads: {}\n", total_count, err_count, unmapped_count, mapped_count);
        if bc_perfect_count + bc_notfound_count + bc_imperfect_count > 0 {
            eprintln!("\tPerfect Barcodes: {}\n\tCorrected Barcodes: {}\n\tUnmatchable Barcodes: {}\n", bc_perfect_count, bc_imperfect_count, bc_notfound_count);
        }
    }
    data.sort();
    let mut out_vec = Vec::new();
    let cdata = data.into_iter().dedup_with_count();
    for (count, record) in cdata {
       let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
        out_vec.push(count_str.as_bytes().to_owned());
    }
    molecule_data.sort();

    return (out_vec, molecule_data);
}

pub fn write_counts (count_vec: Vec<Vec<u8>>, run: &Run) -> Result<(), Box<dyn Error>> {
        let counts_path = &run.level.counts;
        let counts_file = counts_path.to_str().unwrap();
        info!("\t\tWriting counts to : '{}'", counts_file);
        if run.params.verbose{
            eprintln!("Writing counts to : '{}'\n", counts_file);
        }
        let f = File::create(counts_file)?;
        let mut gz = GzBuilder::new()
                        .filename(counts_file)
                        .write(f, Compression::default());
        for result in count_vec {
                gz.write_all(&result)?;
        }
        gz.finish()?;
        Ok(())
}

pub fn write_molecules (molecule_vec: Vec<String>, run: &Run) -> Result<(), Box<dyn Error>> {
        let molecule_path = &run.level.molecule;
        let molecule_file = molecule_path.to_str().unwrap();
        info!("\t\tWriting molecule info to : '{}'", molecule_file);
        if run.params.verbose{
            eprintln!("Writing molecule info to : '{}'\n", molecule_file);
        }
        let f = File::create(molecule_file)?;
        let mut gz = GzBuilder::new()
                        .filename(molecule_file)
                        .write(f, Compression::default());
        for result in molecule_vec {
                gz.write_all(&result.as_bytes())?;
        }
        gz.finish()?;
        Ok(())
}



pub fn cleanup(filename: &Path, warn: bool) -> std::io::Result<()> {
    if Path::new(filename).exists(){
        fs::remove_file(filename.to_str().unwrap())?;
        Ok(())
    }else {
        if warn{
            warn!("\t\tFile does not exist: '{:?}'", filename);
        }
        Ok(())
    }
}



fn pop2(barry: &[u8]) -> &[u8; 2] {
    array_ref!(barry, 0, 2)
}



fn remove_whitespace( s: &str) ->  String{
    s.to_string().retain(|c| !c.is_whitespace());
    return s.to_string()
}

fn get_current_working_dir() -> std::io::Result<PathBuf> {
    env::current_dir()
}

