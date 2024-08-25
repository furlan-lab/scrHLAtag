/**
cd ~/develop/scrHLAtag
cargo build --release && cp target/release/scrHLAtag ~/.local/bin
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/test.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out -l transcriptome -s
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/full_dedup.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/full_corrected_sorted.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out
~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/mini_corrected_sorted.bam -a ~/develop/scrHLAtag/data/testhla.tsv -o out

~/develop/scrHLAtag/target/release/scrHLAtag -v -b ~/develop/scrHLAtag/data/test.bam -o out -l transcriptome -s


#makes all alleles file

zcat < HLA_DB_3field_gene.fa.gz | grep -e ">" - > all_alleles.tmp
zcat < HLA_DB_3field_mRNA.fa.gz | grep -e ">" - >> all_alleles.tmp
R
library(tidyr)
d<-read.table("all_alleles.tmp")[,1] %>% gsub(">", "", .) %>% gsub("\\|", "*", .)
d<-unique(d)
write.table(d, "all_alleles.tsv", row.names=F, quote=F, col.names=F)
q()
n
rm all_alleles.tmp
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

use std::{
    collections::HashMap,
    env,
    error::Error,
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    str,
};

use bam::{BamReader, record::tags::TagValue, record::tags::StringType};
use clap::{App, load_yaml};
use csv::ReaderBuilder;
use fastq::{OwnedRecord, Record};
use flate2::{Compression, GzBuilder};
use itertools::Itertools;
use kseq::parse_path;
use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
use serde::Deserialize;
use simple_log::{info, warn, error, LogConfigBuilder};

#[derive(Deserialize, PartialEq, Ord, PartialOrd, Eq)]
pub struct HLAalleles {
    allele: String,
}

#[derive(Clone)]
pub struct InputParams {
    pub bam: String,
    pub threads: usize,
    pub alleles_file: String,
    pub output_path: Box<Path>,
    pub hla_sep: String,
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

    let mut add_run = |level: &str, hla_ref: &str| {
        runs.push(Run {
            level: AlignmentLevel {
                file_tag: "mRNA".to_string(),
                hla_ref: package_dir_path.join(hla_ref),
                descriptor: level.to_string(),
                mini_fasta: params.output_path.join(format!("align_{}.fa", level)),
                out_unsorted_sam: params.output_path.join(format!("Aligned_mm2_{}.sam", level)),
                out_sorted_sam: params.output_path.join(format!("Aligned_mm2_sorted_{}.sam", level)),
                out_sorted_bam: params.output_path.join(format!("Aligned_mm2_sorted_{}.bam", level)),
                counts: params.output_path.join(format!("counts_{}.txt.gz", level)),
                molecule: params.output_path.join(format!("molecule_info_{}.txt.gz", level)),
            },
            params: params.clone(),
        });
    };

    match params.level.as_str() {
        "transcriptome" => add_run("transcriptome", "data/HLA_DB_3field_mRNA.fa.gz"),
        "genome" => add_run("genome", "data/HLA_DB_3field_gene.fa.gz"),
        "both" => {
            add_run("genome", "data/HLA_DB_3field_gene.fa.gz");
            add_run("transcriptome", "data/HLA_DB_3field_mRNA.fa.gz");
        }
        _ => panic!("Invalid alignment level: {}", params.level),
    }

    runs
}

pub fn load_params() -> Result<InputParams, Box<dyn Error>> {
    let yaml = load_yaml!("../cli.yml");
    let params = App::from_yaml(yaml).get_matches();

    Ok(InputParams {
        bam: params.value_of("bam").unwrap().to_string(),
        threads: params.value_of("threads").unwrap_or("1").parse::<usize>()?,
        alleles_file: params.value_of("alleles_file").unwrap_or("none_provided").to_string(),
        output_path: Path::new(params.value_of("output_folder").unwrap_or("out")).into(),
        hla_sep: params.value_of("hsep").unwrap_or("*").to_string(),
        verbose: params.is_present("verbose"),
        return_sequence: params.is_present("return_sequence"),
        cb_tag: params.value_of("cb").unwrap_or("CB").to_string(),
        umi_tag: params.value_of("umi").unwrap_or("XM").to_string(),
        level: params.value_of("align_level").unwrap_or("both").to_string(),
    })
}

pub fn check_params(params: InputParams) -> Result<InputParams, Box<dyn Error>> {
    cleanup(&params.output_path.join("scrHLAtag.log"), false)?;
    
    let log_config = LogConfigBuilder::builder()
        .path(params.output_path.join("scrHLAtag.log").to_str().unwrap())
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f")
        .level("info")
        .output_file()
        .build();
    simple_log::new(log_config)?;

    let allowables = vec!["genome", "transcriptome", "both"];
    if !allowables.contains(&params.level.as_str()) {
        error!("Input 'level' parameter is not valid: {}", params.level);
        panic!("Invalid 'level' parameter: {}", params.level);
    }

    info!("Starting!");
    let working_dir = get_current_working_dir()?.to_str().unwrap().to_string();
    info!("Current working directory: '{}'", working_dir);

    if params.output_path.is_relative() {
        let abs_outpath = Path::new(&working_dir).join(&params.output_path).to_str().unwrap().to_string();
        if params.output_path.exists() {
            info!("Found existing output directory: '{}'", &abs_outpath);
            warn!("Existing data in this folder could be lost!");
        } else {
            info!("Creating output directory: '{}'", &abs_outpath);
            fs::create_dir(&params.output_path)?;
        }
    } else if params.output_path.exists() {
        info!("Found existing output directory: '{}'", params.output_path.to_str().unwrap());
        warn!("Existing data in this folder could be lost!");
    } else {
        info!("Creating output directory: '{}'", params.output_path.to_str().unwrap());
        fs::create_dir(&params.output_path)?;
    }

    if params.alleles_file == "none_provided" {
        info!("No alleles file provided, running against all alleles!");
        let all_alleles_file = Path::new(env!("CARGO_MANIFEST_DIR")).join("data/all_alleles.tsv");
        return Ok(InputParams { alleles_file: all_alleles_file.to_str().unwrap().to_string(), ..params });
    }

    Ok(params)
}

pub fn read_allelesfile(params: &InputParams) -> Vec<HLAalleles> {
    info!("Opening alleles file: '{}'", params.alleles_file);

    let file = File::open(&params.alleles_file).unwrap();
    let reader = BufReader::new(file);
    let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(reader);

    let mut alleles: Vec<HLAalleles> = rdr.deserialize().filter_map(Result::ok).collect();
    alleles.sort();
    alleles.dedup();

    alleles
}

fn parse_fasta_header (input: String)-> String {
    let chunks: Vec<_> = input.split(" ").collect();
    let substring1 = chunks[0].replace(">", "").replace("|", "*");
    return substring1
}


pub fn make_partial_reference(alleles_query: &[HLAalleles], run: &Run) -> Result<Vec<String>, Box<dyn Error>> {
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut i: usize = 0;

    info!("\t\tRead HLA reference file: '{}'", run.level.hla_ref.to_str().unwrap());
    if run.params.verbose {
        eprintln!("Read HLA reference file: '{}'", run.level.hla_ref.to_str().unwrap());
    }

    let mut records = parse_path(&run.level.hla_ref).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        i = i+1;
        names.insert(parse_fasta_header(record.head().to_string()), i);
    }


    let mut simpleindices: Vec<usize> = Vec::new();
    for data in alleles_query {
        if let Some(&matched) = names.get(&data.allele) {
            info!("\t\t\tFound: {}", &data.allele);
            if run.params.verbose {
                eprintln!("\tFound: {}", &data.allele);
            }
            simpleindices.push(matched);
        } else {
            warn!("\t\tCould not find: {}", &data.allele);
            if run.params.verbose {
                eprintln!("Could not find: {}", &data.allele);
            }
        }
    }

    let align_fasta = run.level.mini_fasta.to_str().unwrap();
    let file = File::create(align_fasta)?;
    let mut fasta_writer = fasta::writer::Builder::default().build_with_writer(BufWriter::new(file));

    let mut names_out: Vec<String> = Vec::new();
    
    let mut records = parse_path(&run.level.hla_ref)?;
    i = 0;
    while let Some(record) = records.iter_record()? {
        i += 1;
        if simpleindices.contains(&i) {
            let header = record.head().to_string();
            names_out.push(header.clone());

            let definition = Definition::new(header, None);
            let sequence = Sequence::from(record.seq().to_string().into_bytes());
            let fasta_record = noodles_fasta::record::Record::new(definition, sequence);

            fasta_writer.write_record(&fasta_record)?;
        }
    }


    Ok(names_out)
}
pub fn test_progs(software: String) -> Result<(), Box<dyn Error>> {
    Command::new(software.clone())
        .arg("-h")
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect(&format!("\n\n*******Failed to execute {}*******\n\n", software));
    Ok(())
}

fn calc_threads(params: &InputParams) -> u16 {
    if params.threads > 2 { params.threads as u16 - 1 } else { 0 }
}

fn convert_tag_value(tag: TagValue) -> String {
    match tag {
        TagValue::Char(v) => v.to_string(),
        TagValue::Int(v, _) => v.to_string(),
        TagValue::Float(v) => v.to_string(),
        TagValue::String(v, _) => str::from_utf8(v).unwrap().to_string(),
        _ => 0.to_string()
    }
}

pub fn make_fastq(params: &InputParams) -> Result<(), Box<dyn Error>> {
    let split = "|BARCODE=".to_string();
    let fastq_path = params.output_path.join("fastq.fq.gz");
    let bam_fn = Path::new(&params.bam);
    let cb_b = pop2(params.cb_tag.as_bytes());
    let umi_b = pop2(params.umi_tag.as_bytes());
    let nb_b = pop2("nb".as_bytes());

    let mut total_count = 0;
    let mut err_count = 0;

    let bam_reader = BamReader::from_path(bam_fn, calc_threads(params))?;
    let f = File::create(&fastq_path)?;
    let mut writer = GzBuilder::new().filename(fastq_path.to_str().unwrap()).write(f, Compression::default());

    for record in bam_reader {
        total_count += 1;
        let rec = record.as_ref().unwrap();

        let old_readname = match str::from_utf8(rec.name()) {
            Ok(value) => remove_whitespace(value),
            Err(_) => {
                err_count += 1;
                continue;
            }
        };

        let cb = match rec.tags().get(&cb_b) {
            Some(TagValue::String(cba, _)) => str::from_utf8(cba).unwrap().to_string(),
            _ => {
                err_count += 1;
                continue;
            }
        };

        let umi = match rec.tags().get(&umi_b) {
            Some(TagValue::String(uba, _)) => str::from_utf8(uba).unwrap().to_string(),
            _ => {
                err_count += 1;
                continue;
            }
        };

        let nb_present = rec.tags().get(&nb_b).is_some();
        let nb = if nb_present { rec.tags().get(&nb_b).unwrap() } else { TagValue::Char(0)  };

        let new_readname = if nb_present {
            format!("{}{}{}_{}_nb_{}", old_readname, split, cb, umi, convert_tag_value(nb))
        } else {
            format!("{}{}{}_{}", old_readname, split, cb, umi)
        };

        let new_record = OwnedRecord {
            head: new_readname.as_bytes().to_vec(),
            seq: rec.sequence().to_vec(),
            sep: None,
            qual: vec![255; rec.sequence().len()],
        };

        new_record.write(&mut writer)?;
    }

    info!("Total reads processed: {}\tReads with errors: {}", total_count, err_count);
    if params.verbose {
        eprintln!("Total reads processed: {}\nReads with errors: {}", total_count, err_count);
    }

    Ok(())
}

pub fn align(run: &Run) -> Result<(), Box<dyn Error>> {
    if !run.level.mini_fasta.exists() {
        error!("Alignment fasta not found at: {}", run.level.mini_fasta.to_str().unwrap());
        panic!("Alignment fasta not found at: {}", run.level.mini_fasta.to_str().unwrap());
    }

    let fastq_path = run.params.output_path.join("fastq.fq.gz");
    if !fastq_path.exists() {
        error!("Alignment fastq not found at: {}", fastq_path.to_str().unwrap());
        panic!("Alignment fastq not found at: {}", fastq_path.to_str().unwrap());
    }

    info!("Aligning reads using minimap2:");
    if run.params.verbose {
        eprintln!("Aligning reads using minimap2:");
    }

    let output = Command::new("minimap2")
        .arg("--cs=long")
        .arg("--secondary=no")
        .arg("-x")
        .arg("map-hifi")
        .arg("-Q")
        .arg("--MD")
        .arg("-a")
        .arg("-t")
        .arg(run.params.threads.to_string())
        .arg(run.level.mini_fasta.to_str().unwrap())
        .arg(fastq_path.to_str().unwrap())
        .arg("-o")
        .arg(run.level.out_unsorted_sam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("Failed to execute minimap2");

    info!("{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }

    Ok(())
}

pub fn sort(run: &Run) -> Result<(), Box<dyn Error>> {
    if !run.level.out_unsorted_sam.exists() {
        error!("SAM not found at: {}", run.level.out_unsorted_sam.to_str().unwrap());
        panic!("SAM not found at: {}", run.level.out_unsorted_sam.to_str().unwrap());
    }

    info!("Minimap2 complete; Running samtools sort");
    if run.params.verbose {
        eprintln!("Minimap2 complete; Running samtools sort");
    }

    let output = Command::new("samtools")
        .arg("sort")
        .arg("-@")
        .arg(run.params.threads.to_string())
        .arg("-o")
        .arg(run.level.out_sorted_sam.to_str().unwrap())
        .arg(run.level.out_unsorted_sam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("Failed to execute samtools sort");

    info!("{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }

    info!("Samtools sort complete; Running samtools view");
    if run.params.verbose {
        eprintln!("Samtools sort complete; Running samtools view");
    }

    if !run.level.out_sorted_sam.exists() {
        error!("Sorted SAM not found at: {}", run.level.out_sorted_sam.to_str().unwrap());
        panic!("Sorted SAM not found at: {}", run.level.out_sorted_sam.to_str().unwrap());
    }

    let output = Command::new("samtools")
        .arg("view")
        .arg("-b")
        .arg("-@")
        .arg(run.params.threads.to_string())
        .arg("-o")
        .arg(run.level.out_sorted_bam.to_str().unwrap())
        .arg(run.level.out_sorted_sam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("Failed to execute samtools view");

    info!("{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }

    info!("Samtools view complete; Running samtools index");
    if run.params.verbose {
        eprintln!("Samtools view complete; Running samtools index");
    }

    if !run.level.out_sorted_bam.exists() {
        error!("Sorted SAM not found at: {}", run.level.out_sorted_bam.to_str().unwrap());
        panic!("Sorted SAM not found at: {}", run.level.out_sorted_bam.to_str().unwrap());
    }

    let output = Command::new("samtools")
        .arg("index")
        .arg("-@")
        .arg(run.params.threads.to_string())
        .arg(run.level.out_sorted_bam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("Failed to execute samtools index");

    info!("{}", String::from_utf8_lossy(&output.stderr));
    if run.params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }

    Ok(())
}

pub fn count(run: &Run) -> (Vec<Vec<u8>>, Vec<String>) {
    info!("Counting reads:");
    if run.params.verbose {
        eprintln!("Counting reads:");
    }

    let split = "|BARCODE=".to_string();
    let bam_reader = BamReader::from_path(&run.level.out_sorted_bam, calc_threads(&run.params)).unwrap();
    let seqnames: Vec<String> = bam_reader
        .header()
        .reference_names()
        .iter()
        .map(|seq| seq.replace("|", &run.params.hla_sep))
        .collect();

    let mut total_count = 0;
    let mut err_count = 0;
    let mut unmapped_count = 0;
    let mut mapped_count = 0;
    let mut bc_perfect_count = 0;
    let mut bc_imperfect_count = 0;
    let mut bc_notfound_count = 0;

    let mut data = Vec::new();
    let mut molecule_data = Vec::new();

    for record in bam_reader {
        total_count += 1;
        let rec = match record.as_ref() {
            Ok(rec) => rec,
            Err(_) => {
                err_count += 1;
                continue;
            }
        };

        let readname = match str::from_utf8(rec.name()) {
            Ok(v) => v,
            Err(_) => {
                err_count += 1;
                continue;
            }
        };

        let mut cbumi = readname.split(&split).nth(1).unwrap_or_default().to_string();
        let mut nb_present = false;
        let mut nb = 0;

        if cbumi.contains("_nb_") {
            nb_present = true;
            nb = cbumi.split("_nb_").nth(1).unwrap().parse::<i32>().unwrap();
            cbumi = cbumi.split("_nb_").nth(0).unwrap().to_string();
        }

        let cb = cbumi.split("_").nth(0);
        let umi = cbumi.split("_").nth(1);

        if rec.ref_id() < 0 {
            unmapped_count += 1;
            continue;
        } else if rec.ref_id() as usize >= seqnames.len() {
            err_count += 1;
            continue;
        } else {
            mapped_count += 1;
            let index = rec.ref_id() as usize;
            data.push(format!("{} {} {}", cb.unwrap_or_default(), umi.unwrap_or_default(), seqnames[index]));

            if nb_present {
                if nb == 0 {
                    bc_perfect_count += 1;
                } else if nb < 0 {
                    bc_notfound_count += 1;
                } else {
                    bc_imperfect_count += 1;
                }
                molecule_data.push(format_molecule_data(
                    rec, &cb.unwrap_or_default(), nb, &umi.unwrap_or_default(), &seqnames[index], run.params.return_sequence,
                ));
            } else {
                molecule_data.push(format_molecule_data(
                    rec, &cb.unwrap_or_default(), 0, &umi.unwrap_or_default(), &seqnames[index], run.params.return_sequence,
                ));
            }
        }
    }

    info!("Total reads processed: {}\tReads with errors: {}", total_count, err_count);
    info!("Unmapped reads: {}\tMapped reads: {}", unmapped_count, mapped_count);

    if bc_perfect_count + bc_notfound_count + bc_imperfect_count > 0 {
        info!(
            "Perfect Barcodes: {}\nCorrected Barcodes: {}\nUnmatchable Barcodes: {}\n",
            bc_perfect_count, bc_imperfect_count, bc_notfound_count
        );
    }

    if run.params.verbose {
        eprintln!(
            "Total reads processed: {}\nReads with errors: {}\nUnmapped reads: {}\nMapped reads: {}\n",
            total_count, err_count, unmapped_count, mapped_count
        );
        if bc_perfect_count + bc_notfound_count + bc_imperfect_count > 0 {
            eprintln!(
                "Perfect Barcodes: {}\nCorrected Barcodes: {}\nUnmatchable Barcodes: {}\n",
                bc_perfect_count, bc_imperfect_count, bc_notfound_count
            );
        }
    }

    let count_vec = data.into_iter().dedup_with_count()
        .map(|(count, record)| format!("{} {}\n", record, count).into_bytes())
        .collect();

    molecule_data.sort();

    (count_vec, molecule_data)
}

fn format_molecule_data(
    rec: &bam::Record, cb: &str, nb: i32, umi: &str, seqname: &str, return_sequence: bool,
) -> String {
    let nm_tag = rec.tags().get(b"NM").unwrap_or(TagValue::String("-1".as_bytes(), StringType::String));
    let as_tag = rec.tags().get(b"AS").unwrap_or(TagValue::Char(0));
    let s1_tag = rec.tags().get(b"s1").unwrap_or(TagValue::Char(0));
    let de_tag = rec.tags().get(b"de").unwrap_or(TagValue::String("0.0".as_bytes(), StringType::String));

    if return_sequence {
        format!(
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n",
            str::from_utf8(rec.name()).unwrap_or_default(),
            cb,
            nb,
            umi,
            seqname,
            rec.query_len(),
            rec.start(),
            rec.mapq(),
            rec.cigar(),
            convert_tag_value(nm_tag),
            convert_tag_value(as_tag),
            convert_tag_value(s1_tag),
            convert_tag_value(de_tag),
            str::from_utf8(&rec.sequence().to_vec()).unwrap_or_default()
        )
    } else {
        format!(
            "{} {} {} {} {} {} {} {} {} {} {} {} {}\n",
            str::from_utf8(rec.name()).unwrap_or_default(),
            cb,
            nb,
            umi,
            seqname,
            rec.query_len(),
            rec.start(),
            rec.mapq(),
            rec.cigar(),
            convert_tag_value(nm_tag),
            convert_tag_value(as_tag),
            convert_tag_value(s1_tag),
            convert_tag_value(de_tag)
        )
    }
}

pub fn write_counts(count_vec: Vec<Vec<u8>>, run: &Run) -> Result<(), Box<dyn Error>> {
    info!("Writing counts to : '{}'", run.level.counts.to_str().unwrap());

    let f = File::create(&run.level.counts)?;
    let mut gz = GzBuilder::new().filename(run.level.counts.to_str().unwrap()).write(f, Compression::default());

    for result in count_vec {
        gz.write_all(&result)?;
    }

    gz.finish()?;
    Ok(())
}

pub fn write_molecules(molecule_vec: Vec<String>, run: &Run) -> Result<(), Box<dyn Error>> {
    info!("Writing molecule info to : '{}'", run.level.molecule.to_str().unwrap());

    let f = File::create(&run.level.molecule)?;
    let mut gz = GzBuilder::new().filename(run.level.molecule.to_str().unwrap()).write(f, Compression::default());

    for result in molecule_vec {
        gz.write_all(result.as_bytes())?;
    }

    gz.finish()?;
    Ok(())
}

pub fn cleanup(filename: &Path, warn: bool) -> std::io::Result<()> {
    if filename.exists() {
        fs::remove_file(filename)?;
    } else if warn {
        warn!("File does not exist: '{:?}'", filename);
    }
    Ok(())
}

fn pop2(barry: &[u8]) -> &[u8; 2] {
    array_ref!(barry, 0, 2)
}

fn remove_whitespace(s: &str) -> String {
    s.chars().filter(|c| !c.is_whitespace()).collect()
}

fn get_current_working_dir() -> std::io::Result<PathBuf> {
    env::current_dir()
}
