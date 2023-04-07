/**

~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/test.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/full_dedup.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out


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





#[derive(Deserialize)]
pub struct HLAalleles {
    allele: String
}


pub struct Params {
    pub bam: String,
    pub threads: usize,
    pub alleles_file: String,
    pub output: Box<Path>,
    pub verbose: bool,
    pub hla_ref: String,
    pub cb_tag: String,
    pub umi_tag: String,
}


pub fn load_params() -> Result<Params, Box<dyn Error>> {
    let wdpb= get_current_working_dir().unwrap();
    let wdir = wdpb.to_str().unwrap();
    let yaml = load_yaml!("params_scrHLAtag.yml");
    let params = App::from_yaml(yaml).get_matches();
    // eprintln!("{:?}", params);
    let bam = params.value_of("bam").unwrap().to_string();
    let alleles_file = params.value_of("alleles_file").unwrap().to_string();
    let cb = params.value_of("c").unwrap_or("CB").to_string();
    let umi = params.value_of("u").unwrap_or("XM").to_string();
    let output = params.value_of("output").unwrap_or("out").to_string();
    let threads = params.value_of("threads").unwrap_or("1").to_string().parse::<usize>().unwrap();
    let mut verbose = false;
    if params.is_present("verbose") {
        verbose = true;
    }
    let package_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let hla_path = package_dir.join("data/hla_mRNA.fasta.gz");
    let hla_ref = hla_path.to_str().unwrap();
    let _cu = cleanup(&Path::new(&output).join("scrHLAtag.log"), false);
    let config = LogConfigBuilder::builder()
        .path(Path::new(&output).join("scrHLAtag.log").to_string_lossy())
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("info")
        .output_file()
        .build();
    let _ = simple_log::new(config);
    info!("\t\t\tStarting!");
    info!("\t\t\tCurrent working directory: '{}'", wdir);
    if verbose {
        eprintln!("\n\nCurrent working directory: '{}'", wdir);
    }
    let outpath = Path::new(&output);
    if outpath.is_relative(){
        let a1 = Path::new(wdir).join(&output);
        let abs_outpath = a1.to_str().unwrap();
        if outpath.exists() {
                info!("\t\t\tFound existing output directory: '{}'", &abs_outpath);
                warn!("\t\t\t{}", "Existing data in this folder could be lost!!!");
               if verbose {
                    eprintln!("Found existing output directory: '{}'", &abs_outpath);
                    eprintln!("\t{}", "Existing data in this folder could be lost!!!");
                }
        } else {
            info!("\t\tCreating output directory: '{}'", &abs_outpath);
            if verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(outpath)?;
        }
    } else {
        let abs_outpath = outpath.to_str().unwrap();
        if outpath.exists() {
            info!("\t\tFound existing output directory: '{}'", &abs_outpath);
            warn!("\t\t{}", "Existing data in this folder could be lost!!!");
            if verbose {
                eprintln!("Found existing output directory: '{}'", &abs_outpath);
                eprintln!("\t{}", "Existing data in this folder could be lost!!!");
            }
        } else {
            info!("\t\tCreating output directory: '{}'", &abs_outpath);
            if verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(outpath)?;
        }
    }
    Ok(Params{
            bam: bam,
            threads: threads as usize,
            alleles_file: alleles_file,
            output: outpath.into(),
            verbose: verbose,
            hla_ref: hla_ref.to_string(),
            cb_tag: cb,
            umi_tag: umi,
    })
}


pub fn read_allelesfile(params: &Params) -> Vec<HLAalleles> {
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
pub fn make_partial_reference (alleles_query: Vec<HLAalleles>, params: &Params) -> Result<Vec<String>, Box<dyn Error>>{
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut i: usize = 0;
    info!("\t\tRead HLA reference file: '{}'", &params.hla_ref);
    if params.verbose {
        eprintln!("Read HLA reference file: '{}'", &params.hla_ref);
    }

    let mut records = parse_path(&params.hla_ref).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        i = i+1;
        names.insert(parse_fasta_header(record.head().to_string()), i);
    }

    let mut simpleindices = Vec::new();
    for data in alleles_query {
        let matched = names.get(&data.allele);
        if matched.is_some(){
            info!("\t\t\tFound: {}", &data.allele);
            if params.verbose {
                eprintln!("\tFound: {}", &data.allele);
            }
            simpleindices.push(*matched.unwrap())
        } else{ 
            warn!("\t\tCould not find: {}", &data.allele);
            if params.verbose {
                eprintln!("Could not find: {}", &data.allele)
            }
        };
    }
    i = 0;
    let align_fasta_path = params.output.join("align.fa");
    let align_fasta = align_fasta_path.to_str().unwrap();
    let f = File::create(align_fasta)?;

    // write align.fasta
    let mut fasta_writer = fasta::writer::Builder::default().build_with_writer(BufWriter::new(f));
    let mut names_out = Vec::new();
    let mut records = parse_path(&params.hla_ref).unwrap();
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

fn calc_threads(params: &Params) -> u16 {
    if params.threads>2{
        (params.threads-1).try_into().unwrap()
    } else {
        0
    }

}

#[allow(unused_assignments)]
pub fn make_fastq (params: &Params)-> Result<(), Box<dyn Error>> {

    // declare local variables
    let split = "|BARCODE=".to_string();
    let fastq_path = &params.output.join("fastq.fq.gz");
    let bam_fn = Path::new(&params.bam);
    // let cb_b: [u8; 2] = clone_into_array(&params.cb_tag.as_bytes().to_vec()[0..1]);
    let binding = params.cb_tag.as_bytes().to_vec();
    let cb_b = pop2(&binding);
    let binding = params.umi_tag.as_bytes().to_vec();
    let umi_b = pop2(&binding);

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
        let new_readname = format!("{}{}{}_{}",&old_readname, &split, cb, umi );
        
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

pub fn align (params: &Params)-> Result<(), Box<dyn Error>> {
    let align_fasta_path = params.output.join("align.fa");
    let align_fasta = align_fasta_path.to_str().unwrap();
    if !align_fasta_path.exists() {
        error!("\t\tAlignment fasta not found at: {}", align_fasta);
        panic!("Alignment fasta not found at: {}", align_fasta);
    }
    let fastq_path = params.output.join("fastq.fq.gz");
    let fastq_file = fastq_path.to_str().unwrap();
    if !fastq_path.exists() {
        error!("\t\tAlignment fastq not found at: {}", fastq_file);
        panic!("Alignment fastq not found at: {}", fastq_file);
    }
    let sam_path = params.output.join("Aligned_mm2.sam");
    let sam_file = sam_path.to_str().unwrap();
    info!("\t\t{}", "Aligning reads using minimap2 - Output below:");
    if params.verbose {
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
                    .arg(params.threads.to_string())
                    .arg(align_fasta)
                    .arg(fastq_file)
                    .arg("-o")
                    .arg(sam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute minimap2*******\n\n");
    info!("{}", String::from_utf8_lossy(&output.stderr));
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    Ok(())
}


pub fn sort (params: &Params)-> Result<(), Box<dyn Error>> {
    let sam_path = params.output.join("Aligned_mm2.sam");
    let sam_file = sam_path.to_str().unwrap();
    if !sam_path.exists() {
        error!("\t\tSAM not found at: {}", sam_file);
        panic!("SAM not found at: {}", sam_file);
    }
    let ssam_path = params.output.join("Aligned_mm2_sorted.sam");
    let ssam_file = ssam_path.to_str().unwrap();
    let bam_path = params.output.join("Aligned_mm2_sorted.bam");
    let bam_file = bam_path.to_str().unwrap();
    info!("\t\t{}", "Minimap2 complete; Running samtools sort");
    if params.verbose {
        eprintln!("{}", "Minimap2 complete; Running samtools sort");
    }
    let output = Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg(ssam_file)
                    .arg(sam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools view*******\n\n");
    info!("\t\t{}", String::from_utf8_lossy(&output.stderr));
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
        
    }
    info!("\t\t{}", "Samtools sort complete; Running samtools view");
    if params.verbose {
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
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg(bam_file)
                    .arg(ssam_file)
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools sort*******\n\n");
    info!("\t\t{}", String::from_utf8_lossy(&output.stderr));
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    info!("\t\t{}", "Samtools view complete; Running samtools index");
    if params.verbose {
        eprintln!("{}", "Samtools view complete; Running samtools index");
    }
    if !bam_path.exists() {
        error!("Sorted SAM not found at: {}", bam_file);
        panic!("\t\tSorted SAM not found at: {}", bam_file);
    }
    let output = Command::new("samtools")
                    .arg("index")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg(bam_file)
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools index*******\n\n");
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    Ok(())
}

pub fn count(params: &Params) -> (Vec<Vec<u8>>, Vec<String>){
    info!("\t\tCounting reads:");
    if params.verbose {
        eprintln!("Counting reads:");
    }
    // declare some stuff
    let split = "|BARCODE=".to_string();
    let bam_path = params.output.join("Aligned_mm2_sorted.bam");
    let bam_file = bam_path.to_str().unwrap();
    let mut total_count: usize = 0;
    let mut err_count: usize = 0;
    let mut unmapped_count: usize = 0;
    let mut mapped_count: usize = 0;

    // bam reader and get seqnames
    let bam_reader = bam::BamReader::from_path(bam_file, calc_threads(params)).unwrap();
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
        let cbumi = readname.split(&split).nth(1).unwrap();
        // eprintln!("{}", &cbumi);
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
                            warn!("AS tag not returned correctly for read {:?}", rec.name());
                            0
                        },
                };
                let de_tag = match rec.tags().get(b"de") {
                    Some(TagValue::Float(value)) => value,
                    _ => {
                            warn!("AS tag not returned correctly for read {:?}", rec.name());
                            0.0
                        },
                };
                // columns cb, umi, seqname, star, mapq, cigar, NM, AS, chaining_score, de (per base sequence divergence)
                molecule_data.push(format!("{} {} {} {} {} {} {} {} {} {}\n", &cb.unwrap(), &umi.unwrap(), seqnames[index], rec.start(), rec.mapq(), rec.cigar(), nm_tag, as_tag, s1_tag, de_tag));
                
        }
    }
    info!("\t\t\tTotal reads processed: {}\tReads with errors: {}", total_count, err_count);
    info!("\t\t\tUnmapped reads: {}\tMapped reads: {}", unmapped_count, mapped_count);
    if params.verbose {
        eprintln!("\tTotal reads processed: {}\n\tReads with errors: {}\n\tUnmapped reads: {}\n\tMapped reads: {}\n", total_count, err_count, unmapped_count, mapped_count);
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

pub fn write_counts (count_vec: Vec<Vec<u8>>, params: &Params) -> Result<(), Box<dyn Error>> {
        let counts_path = params.output.join("counts.txt.gz");
        let counts_file = counts_path.to_str().unwrap();
        info!("\t\tWriting counts to : '{}'", counts_file);
        if params.verbose{
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

pub fn write_molecules (molecule_vec: Vec<String>, params: &Params) -> Result<(), Box<dyn Error>> {
        let molecule_path = params.output.join("molecule_txt.gz");
        let molecule_file = molecule_path.to_str().unwrap();
        info!("\t\tWriting molecule info to : '{}'", molecule_file);
        if params.verbose{
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

