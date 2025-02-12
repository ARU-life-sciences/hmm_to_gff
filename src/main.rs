// Format the output of hmmscan/cmscan to a GFF file
// Specifically we are interested in the output from ORFfinder
// fasta -> ORFfinder -> hmmscan -> GFF
// We use the Pfam HMM database to achieve this.
//
// OR:
//
// fasta -> cmscan -> GFF
//
// We use RFam for this.

use calm_io::*;
use hmm_tblout::Reader;

// program used is ORFfinder
const SOURCE: [&str; 2] = ["ORFfinder", "cmscan"];
// feature type is CDS
const FEATURE: [&str; 2] = ["CDS", "intron"];
// frames are not available
const FRAME: &str = ".";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // load in the command line arguments
    let args: Vec<String> = std::env::args().collect();
    // we only expect one arg, the path to the HMM file
    if args.len() != 3 {
        let _ = stderrln!("hmm_to_gff by Max Brown");
        let _ = stderrln!("Usage: {} <HMM tblout file> <hmmscan/cmscan>", args[0]);
        std::process::exit(1);
    }

    let file_path = &args[1];
    let hmmer_type = &args[2];

    if hmmer_type != "hmmscan" && hmmer_type != "cmscan" {
        let _ = stderrln!("Usage: {} <HMM tblout file> <**hmmscan/cmscan**>", args[0]);
        std::process::exit(1);
    }

    // now create a HMMER tblout reader
    let mut reader = Reader::from_path(file_path)?;

    let mut gff_records = Vec::new();

    match hmmer_type.as_str() {
        "hmmscan" => {
            // hmmscan output
            process_hmmscan_from_orffinder(&mut reader, &mut gff_records)?;
        }
        "cmscan" => {
            // cmscan output
            process_cmscan_from_infernal(&mut reader, &mut gff_records)?;
        }
        _ => {
            // this should never happen
            let _ = stderrln!("Invalid HMMER type");
            std::process::exit(1);
        }
    }

    Ok(())
}

fn process_hmmscan_from_orffinder<R: std::io::Read>(
    reader: &mut Reader<R>,
    gff_records: &mut Vec<(String, &str, &str, i32, i32, f32, char, &str, String)>,
) -> Result<(), Box<dyn std::error::Error>> {
    for record in reader.records() {
        let record = record?;

        let query = record.query_name();

        // format of the query name is:
        // lcl|ORF26_u14:7416:7751
        // get name of the ORF, start, end
        let left = query.split("|").last().unwrap().to_string();
        // remove everything up to and including the _
        let seqname = left.split(":").next().unwrap().to_string();
        let seqname = seqname.split("_").last().unwrap().to_string();

        let right = query.split("|").last().unwrap();
        let name: Vec<&str> = right.split(":").collect();
        let start = name[1].parse::<u64>()?;
        let end = name[2].parse::<u64>()?;

        let score = record.score_full().unwrap();

        // sort out starts and ends
        let final_start;
        let final_end;

        let strand = if start < end {
            final_start = start;
            final_end = end;
            '+'
        } else {
            final_start = end;
            final_end = start;
            '-'
        };

        // and we make the attribute here in the form:
        // query_name_full=;accession=;target_name=;description=;
        let attr = format!(
            "query_name={};pfam_accession={};target_name={};description={};e_value={:.2e}",
            record.query_name(),
            record.target_accession(),
            record.target_name(),
            record.description().replace(" ", "_"),
            record.e_value_full().unwrap(),
        );

        gff_records.push((
            seqname,
            SOURCE[0],
            FEATURE[0],
            final_start as i32,
            final_end as i32,
            score,
            strand,
            FRAME,
            attr,
        ));
    }

    // sort the GFF on sequence name and start position
    gff_records.sort_by(|a, b| {
        if a.0 == b.0 {
            a.3.cmp(&b.3)
        } else {
            a.0.cmp(&b.0)
        }
    });

    let _ = stdoutln!("##gff-version 3");

    for record in gff_records {
        let gff_line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.0,
            record.1,
            record.2,
            record.3,
            record.4,
            record.5,
            record.6,
            record.7,
            record.8
        );
        let _ = stdoutln!("{}", gff_line);
    }

    Ok(())
}

fn process_cmscan_from_infernal<R: std::io::Read>(
    reader: &mut Reader<R>,
    gff_records: &mut Vec<(String, &str, &str, i32, i32, f32, char, &str, String)>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut id = 1;
    for record in reader.records() {
        let record = record?;

        let seqid = record.query_name();
        let source = SOURCE[1];
        let feature_type = FEATURE[1];

        let start = record.seq_from().unwrap();
        let end = record.seq_to().unwrap();

        // account for the fact that sometimes the start is larger than the end
        let updated_start = if start < end { start } else { end };
        let updated_end = if start < end { end } else { start };

        let score = record.score().unwrap();

        let strand = record.strand().unwrap();

        // Rfam accession, description
        let attributes = format!(
            "ID={}_{};rfam_accession={};description={};e_value={:.2e}",
            seqid,
            id,
            record.target_accession(),
            record.description().replace(" ", "_"),
            record.e_value().unwrap()
        );

        gff_records.push((
            seqid,
            source,
            feature_type,
            updated_start,
            updated_end,
            score,
            // hacky af
            strand.to_string().chars().next().unwrap(),
            FRAME,
            attributes,
        ));

        id += 1;
    }

    // sort the GFF on sequence name and start position
    gff_records.sort_by(|a, b| {
        if a.0 == b.0 {
            a.3.cmp(&b.3)
        } else {
            a.0.cmp(&b.0)
        }
    });

    let _ = stdoutln!("##gff-version 3");

    for record in gff_records {
        let gff_line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.0,
            record.1,
            record.2,
            record.3,
            record.4,
            record.5,
            record.6,
            record.7,
            record.8
        );
        let _ = stdoutln!("{}", gff_line);
    }

    Ok(())
}
