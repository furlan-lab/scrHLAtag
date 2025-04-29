
extern crate clap;
extern crate bam;

use bam::RecordWriter;
use clap::{App, load_yaml};
use simple_log::LogConfigBuilder;
use simple_log::info;
// use bam::{BamReader, BamWriter};

pub mod scrhlatag;
// use crate::scrhlatag::remove_whitespace;
// use std::path::Path;

fn main() {
    let mut verbose = false;
    let yaml = load_yaml!("addback.yml");
    let params = App::from_yaml(yaml).get_matches();
    // eprintln!("{:?}", params);
    let bam = params.value_of("bam").unwrap().to_string();
    if params.is_present("verbose") {
        verbose = true;
    }
    if verbose {
        eprintln!("Writing log file: '{}'", "addback.log");
    }
    let out = "out.bam".to_string();

    let config = LogConfigBuilder::builder()
        .path("addback.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("info")
        .output_file()
        .build();
    let _ = simple_log::new(config);

    addback_bam(bam, out);
    if verbose {
        eprintln!("Done!!!");
    }
    info!("\t\t\t\tDone!!!");

}

pub fn addback_bam(bam: String, out: String) {
    let split = "|BARCODE=".to_string();
    let reader = bam::BamReader::from_path(bam, 0).unwrap();
    let mut writer = bam::BamWriter::from_path(out, reader.header().clone()).unwrap();
    for record in reader {
        let mut record = record.unwrap();
        // get name
        let _new_readname = match std::str::from_utf8(record.name()) {
            Ok(value) => {
                // remove_whitespace(value)
                
                let data = value.split(split.as_str()).collect::<Vec<&str>>()[1].to_string();
                let data_vec = data.split("_").collect::<Vec<&str>>();
                let cb = data_vec[0].to_string();
                let umi = data_vec[1].to_string();
                let nb = data_vec[3].parse::<i32>().unwrap_or(0);
                record.tags_mut().push_string(b"CB", cb.as_bytes());
                record.tags_mut().push_string(b"XM", umi.as_bytes());
                record.tags_mut().push_num(b"nb", nb);

            },
            Err(_e) => {
                continue
            }
        };
        let _write = writer.write(&record);
    }
}