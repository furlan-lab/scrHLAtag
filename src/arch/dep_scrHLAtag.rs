#![allow(non_snake_case)]
/**

target/release/scrHLAtag -b data/test.bam -a data/testhla.tsv

**/
// scrHLA typing, alignment, : single cell rna-based HLA typing and alignment
// this program takes a pacbio bam and performs alignment and prediction
// single cell hla typing, alignment, and grouping - scrHLAtag


extern crate csv;
extern crate clap;
extern crate serde;
extern crate minimap2;
extern crate noodles;

use std::{io, str, fs, fs::File, error::Error, path::Path, collections::HashMap, process::{Command, Stdio }};
use std::io::{BufReader, BufWriter};
use clap::{App, load_yaml};
use serde::Deserialize;
use csv::ReaderBuilder;
use simple_log::{info, LogConfigBuilder};
use rust_htslib::{bam, bam::{Read, Header, Record}, bam::header::HeaderRecord, bam::record::{Cigar, CigarString}, };
use minimap2::{Aligner, Mapping};
use noodles::fasta;


#[derive(Deserialize)]
struct HLAalleles {
    allele: String
}


struct Params {
    bam: String,
    // genome: String,
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
        // let genome = params.value_of("genome").unwrap_or("genome.fasta");
        let threads = params.value_of("threads").unwrap_or("1");
        let threads = threads.to_string().parse::<usize>().unwrap();
        let mut verbose = true;
        if params.is_present("verbose") {
                verbose = false
        };

        Params{
                bam: bam.to_string(),
                // genome: genome.to_string(),
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
fn make_partial_reference (ref_file: String, alleles_query: Vec<HLAalleles>) -> Result<Vec<String>, Box<dyn Error>>{
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut i: usize = 0;
    // let mut reader = fasta::indexed_reader::Builder::default().build_from_path(src)?;
    // TODO: use compressed fasta: https://github.com/zaeleus/noodles/issues/142
    let mut fasta_reader = File::open(&ref_file)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    for result in fasta_reader.records() {
        i = i+1;
        let record = result?;
        names.insert(parse_fasta_header(record.name().to_string()), i);
    }
    let mut simpleindices = Vec::new();
    for data in alleles_query {
        let matched = names.get(&data.allele);
        if matched.is_some(){
            eprintln!("Found: {}", &data.allele);
            simpleindices.push(*matched.unwrap())
        } else{ 
            eprintln!("Could not find: {}", &data.allele)
        };
    }
    i = 0;
    let align_fasta = "data/align.fasta".to_string();
    let mut fasta_reader = File::open(&ref_file)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;
    let f = File::create(align_fasta)?;
    let mut fasta_writer = fasta::writer::Builder::default().build_with_writer(BufWriter::new(f));
    // let mut output =  fasta::Writer::new(Vec::new());
    let mut names_out = Vec::new();
    for result in fasta_reader.records() {
        i = i+1;
        let record = result?;
        if simpleindices.contains(&i){
            names_out.push(record.name().to_string());
            let _fw = fasta_writer.write_record(&record);
        }
    }
    Ok(names_out.to_owned())
}


fn parse_fasta_header (input: String)-> String {
    let chunks: Vec<_> = input.split(" ").collect();
    let substring1 = chunks[0].replace(">", "").replace("|", "*");
    return substring1
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
        eprintln!("\nParsing Parameters!\n");
        eprintln!("\nRunning with {} threads!\n", &params.threads);
    }
    if params.verbose {
        eprintln!("\nChecking programs and parsing alleles_file!\n");
    }
    let _prog_test_res = test_progs();
    let alleles_query = read_allelesfile(&params);
    
    if params.verbose {
        eprintln!("\nMatching HLA alleles with known references!\n");
    }

    let ref_file = "data/hla_mRNA.fasta".to_string();
    if params.verbose {
        eprintln!("\nChecking programs, parsing reference file, and making partial reference: {}\n", &ref_file);
    }
    let fasta_names = make_partial_reference(ref_file, alleles_query);
    
    let _ar = align_and_count(&params, fasta_names.unwrap());
    // let _cu = cleanup("data/align.fasta".to_string());
    // let _so = sort(&params);

}

fn cleanup(filename: String) -> std::io::Result<()> {
    fs::remove_file(filename)?;
    // fs::remove_file("out.bam")?;
    Ok(())
}


// minimap2 --MD -a $fa -t 8 mutcaller_R1.fq.gz -o Aligned.mm2.sam
// samtools sort -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sam
// samtools view -b -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sorted.bam
// samtools index -@ 8 Aligned.mm2.ssorted.bam

fn test_progs () -> Result<(), Box<dyn Error>>{
    let _output = Command::new("samtools")
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute samtools*******\n\n");
    Ok(())
}




fn align_and_count (params: &Params, fasta_names: Vec<String>)-> Result<(), Box<dyn Error>> {
    let aligner = Aligner::builder()
        .map_hifi()
        .with_threads(1)
        .with_cigar()
        .with_index("data/align.fasta", None)
        .expect("Unable to build index");
    // let thresh: u8= 30;
    let mut bam = bam::Reader::from_path(&params.bam).unwrap();
    let mut header = Header::from_template(bam.header());
    for name in &fasta_names {
        let mut new_header_record = HeaderRecord::new(b"SQ");
        new_header_record.push_tag(b"SN", &name);
        new_header_record.push_tag(b"LN", &1000); // TODO: get length correct
        Header::push_record(& mut header, &new_header_record);
    }
   let mut writer = bam::Writer::from_stdout(& mut header, bam::Format::Sam).unwrap();
    // let mut writer = bam::Writer::from_path(Path::new(&params.output), & mut header, bam::Format::Sam).unwrap();
    let mut i = 0;
    // let f = File::create(align_fasta)?;
    // let mut counts_writer = BufWriter::new(f);
    let mut count_writer = io::stdout();
    for r in bam.records() {
        i = i +1;
        let record = r.unwrap();
        let mapping = aligner
            .map(&record.seq().as_bytes(), false, false, None, None)
            .expect("Unable to align");
        eprintln!("{:?}", &mapping);
        if mapping.len() > 0 {
            let first_mapping = mapping.first().unwrap();
            let mut aligned_record = record.clone();
            aligned_record = add_mapping_to_record(Some(&first_mapping), aligned_record, &fasta_names);
            let _w = writer.write(&aligned_record);
        } else {
            let _w = writer.write(&record);
        }
        // cb = aligned.record.aux(b"CB")
        // out_writer.write(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, _result))

    }
    Ok(())
}




pub fn add_mapping_to_record(
    mapping: Option<&Mapping>, rec: Record, tid_strings: &Vec<String>) -> Record {
    let mut newrec = Record::new();
    let qname = rec.qname();
    let qual =  vec![255; rec.seq().len()];
    let cigar: Option<CigarString> = mapping
        .and_then(|m| m.alignment.clone()) // FIXFIX: we probably don't need a clone here
        .and_then(|a| a.cigar)
        .map(|c| cigar_to_cigarstr(&c));
    newrec.set(qname, cigar.as_ref(), &rec.seq().as_bytes(), &qual[..]); 
    match mapping {
        Some(m) => {
            // println!("Strand {m:?}");
            // if m.strand == Strand::Reverse {
            //     println!("here");
            //     rec.set_reverse();
            // }
            // TODO: set secondary/supplementary flags
            let tid = tid_strings.iter().position(|r| r == &m.target_name.clone().unwrap()).unwrap() as i32;
            newrec.set_tid(tid);
            // /** 
            //  * query_len
            //  * query_start
            //  * query_end
            //  * Strand
            //  * target_len
            //  * target_end
            //  * match_len
            //  * block_len
            //  * md?
            //  * cs?
            //  * 
            //  **/
            newrec.set_pos(m.target_start as i64);
            newrec.set_mapq(m.mapq as u8);
            newrec.set_mpos(-1);
            newrec.set_mtid(-1);
            newrec.set_insert_size(0);
        }
        None => {
            newrec.set_unmapped();
            newrec.set_tid(-1);
            newrec.set_pos(-1);
            newrec.set_mapq(255);
            newrec.set_mpos(-1);
            newrec.set_mtid(-1);
            newrec.set_insert_size(-1);
        }
    };
    for aux_data in rec.aux_iter(){
        let tag = aux_data.unwrap();
        let _s = newrec.push_aux(tag.0, tag.1);
    }
    // TODO: set AUX flags for cs/md if available
    newrec
}


fn cigar_to_cigarstr(cigar: &Vec<(u32, u8)>) -> CigarString {
    let op_vec: Vec<Cigar> = cigar
        .to_owned()
        .iter()
        .map(|(len, op)| match op {
            0 => Cigar::Match(*len),
            1 => Cigar::Ins(*len),
            2 => Cigar::Del(*len),
            3 => Cigar::RefSkip(*len),
            4 => Cigar::SoftClip(*len),
            5 => Cigar::HardClip(*len),
            6 => Cigar::Pad(*len),
            7 => Cigar::Equal(*len),
            8 => Cigar::Diff(*len),
            _ => panic!("Unexpected cigar operation"),
        })
        .collect();
    CigarString(op_vec)
}


fn sort (params: &Params)-> Result<(), Box<dyn Error>> {
    let keep = false;
    let output = Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg("out.sorted.bam")
                    .arg(params.output.to_string())
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools view*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    eprintln!("{}", "Samtools sort complete; Running samtools index");
    let output = Command::new("samtools")
                    .arg("index")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("out.sorted.bam")
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools index*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    if keep {
        fs::remove_file("out.sorted.bam")?;
    }
    Ok(())
}



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

