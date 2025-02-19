// Format the output of hmmscan/cmscan/nhmmscan to a GFF file
// Specifically we are interested in the output from ORFfinder
// fasta -> ORFfinder -> hmmscan -> GFF
// We use the Pfam HMM database to achieve this.
//
// OR:
//
// fasta -> cmscan -> GFF
//
// We use RFam for this.
//
// This is kind of a hacky script which needs to be tidied up to be made more
// generic than my specific use case.

use std::{collections::HashMap, fmt::Display, str::FromStr};

use calm_io::*;
use hmm_tblout::Reader;

enum Source {
    Hmmscan,
    Cmscan,
    Nhmmscan,
}

impl Display for Source {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let source = match self {
            // specific to orffinder
            Source::Hmmscan => "ORFfinder",
            Source::Cmscan => "cmscan",
            // specific to oatkdb output
            Source::Nhmmscan => "oatkDB",
        };
        write!(f, "{}", source)
    }
}

impl FromStr for Source {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "hmmscan" => Ok(Source::Hmmscan),
            "cmscan" => Ok(Source::Cmscan),
            "nhmmscan" => Ok(Source::Nhmmscan),
            e => Err(From::from(format!(
                "Found {}, second arg must be <hmmscan/nhmmscan/cmscan>",
                e
            ))),
        }
    }
}

enum Features {
    CDS,
    Intron,
    Gene,
    Trna,
    Rrna,
}

impl Display for Features {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let feature = match self {
            Features::CDS => "CDS",
            Features::Intron => "intron",
            Features::Gene => "gene",
            Features::Trna => "tRNA",
            Features::Rrna => "rRNA",
        };
        write!(f, "{}", feature)
    }
}

// frames are not available
const FRAME: &str = ".";
// matches are partial if they are less than 80% the length of the gene/hit
const COVERAGE: f32 = 0.8;

// now include gene: length map
const PFAM_DATA: &str = include_str!("../data/Pfam-A.lengths");
const RFAM_INTRONS: &str = include_str!("../data/introns.length");
const VIRIDIPLANTAE_DATA: &str = include_str!("../data/viridiplantae_mito_v20250217.lengths");

fn pfam_lengths() -> HashMap<&'static str, usize> {
    PFAM_DATA
        .lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            let gene = parts.next()?;
            let length = parts.next()?.parse::<usize>().ok()? * 3; // *3 to account for aa -> nt. Note this isn't perfect (frame shifts, start/stop codons/non coding regions...)
            Some((gene, length))
        })
        .collect()
}

fn rfam_lengths() -> HashMap<&'static str, usize> {
    RFAM_INTRONS
        .lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            let gene = parts.next()?;
            let length = parts.next()?.parse::<usize>().ok()?;
            Some((gene, length))
        })
        .collect()
}

fn viridiplantae_lengths() -> HashMap<&'static str, usize> {
    VIRIDIPLANTAE_DATA
        .lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            let gene = parts.next()?;
            let length = parts.next()?.parse::<usize>().ok()?;
            Some((gene, length))
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // load in the command line arguments
    let args: Vec<String> = std::env::args().collect();
    // we only expect one arg, the path to the HMM file
    if args.len() != 3 {
        let _ = stderrln!("hmm_to_gff by Max Brown");
        let _ = stderrln!(
            "Usage: {} <HMM tblout file> <hmmscan/cmscan/nhmmscan>",
            args[0]
        );
        let _ = stderrln!("Note hmmscan output is specific for ORFfinder output.");
        std::process::exit(1);
    }

    let file_path = &args[1];
    let source = Source::from_str(&args[2])?;

    // now create a HMMER tblout reader
    let mut reader = Reader::from_path(file_path)?;

    let mut gff_records = Vec::new();

    match source {
        Source::Hmmscan => {
            // hmmscan output
            process_hmmscan_from_orffinder(&mut reader, &mut gff_records)?;
        }
        Source::Cmscan => {
            // cmscan output
            process_cmscan_from_infernal(&mut reader, &mut gff_records)?;
        }
        Source::Nhmmscan => {
            process_nhmmscan_from_oatkdb(&mut reader, &mut gff_records)?;
        }
    }

    Ok(())
}

fn process_hmmscan_from_orffinder<R: std::io::Read>(
    reader: &mut Reader<R>,
    gff_records: &mut Vec<(String, String, String, i32, i32, f32, char, &str, String)>,
) -> Result<(), Box<dyn std::error::Error>> {
    // get the lengths of each of the genes
    let pfam_lengths = pfam_lengths();

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

        // do the length analysis
        let aligned_length = (final_end - final_start) as f32;
        let gene_length = *pfam_lengths.get(record.target_name().as_str()).unwrap() as f32;

        // make the coverage a descriptive tag, with % coverage included
        // if partial
        let coverage = if aligned_length / gene_length < COVERAGE {
            format!("partial:{:.2}", (aligned_length / gene_length) * 100.0)
        } else {
            "full".to_string()
        };

        // and we make the attribute here in the form:
        // query_name_full=;accession=;target_name=;description=;
        let attr = format!(
            "query_name={};pfam_accession={};target_name={};description={};coverage={};e_value={:.2e}",
            record.query_name(),
            record.target_accession(),
            record.target_name(),
            record.description().replace(" ", "_"),
            coverage,
            record.e_value_full().unwrap(),
        );

        gff_records.push((
            seqname,
            Source::Hmmscan.to_string(),
            Features::CDS.to_string(),
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
    gff_records: &mut Vec<(String, String, String, i32, i32, f32, char, &str, String)>,
) -> Result<(), Box<dyn std::error::Error>> {
    let rfam_intron_lengths = rfam_lengths();
    let mut id = 1;
    for record in reader.records() {
        let record = record?;

        // skip records that don't meet the significance threshold
        // i.e. below 1e-05

        if record.e_value().unwrap() > 1e-05 {
            continue;
        }

        let seqid = record.query_name();
        let source = Source::Cmscan.to_string();
        let feature_type = Features::Intron.to_string();

        let start = record.seq_from().unwrap();
        let end = record.seq_to().unwrap();

        // account for the fact that sometimes the start is larger than the end
        let updated_start = if start < end { start } else { end };
        let updated_end = if start < end { end } else { start };

        let score = record.score().unwrap();

        let strand = record.strand().unwrap();

        // do the length analysis
        let aligned_length = (updated_end - updated_start) as f32;
        let gene_length = *rfam_intron_lengths
            .get(record.target_name().as_str())
            .unwrap() as f32;

        // make the coverage a descriptive tag, with % coverage included
        // if partial
        let coverage = if aligned_length / gene_length < COVERAGE {
            format!("partial:{:.2}", (aligned_length / gene_length) * 100.0)
        } else {
            "full".to_string()
        };

        // Rfam accession, description
        let attributes = format!(
            "ID={}_{};rfam_accession={};description={};coverage={};e_value={:.2e}",
            seqid,
            id,
            record.target_accession(),
            record.description().replace(" ", "_"),
            coverage,
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

fn process_nhmmscan_from_oatkdb<R: std::io::Read>(
    reader: &mut Reader<R>,
    gff_records: &mut Vec<(String, String, String, i32, i32, f32, char, &str, String)>,
) -> Result<(), Box<dyn std::error::Error>> {
    let viridiplantae_lengths = viridiplantae_lengths();

    let mut id = 1;
    for record in reader.records() {
        let record = record?;

        // skip records that don't meet the significance threshold
        // i.e. below 1e-05

        if record.e_value().unwrap() > 1e-05 {
            continue;
        }

        let seqid = record.query_name();
        let genid = record.target_name();
        let source = Source::Nhmmscan.to_string();

        let feature_type = match record.target_name().contains("trn") {
            true => Features::Trna.to_string(),
            false => match record.target_name().contains("rrn") {
                true => Features::Rrna.to_string(),
                false => Features::Gene.to_string(),
            },
        };

        let start = record.env_from().unwrap();
        let end = record.env_to().unwrap();

        // account for the fact that sometimes the start is larger than the end
        let updated_start = if start < end { start } else { end };
        let updated_end = if start < end { end } else { start };

        let score = record.score().unwrap();

        let strand = record.strand().unwrap();

        // do the length analysis
        let aligned_length = (updated_end - updated_start) as f32;
        let gene_length = *viridiplantae_lengths
            .get(record.target_name().as_str())
            .unwrap() as f32;

        // make the coverage a descriptive tag, with % coverage included
        // if partial
        let coverage = if aligned_length / gene_length < COVERAGE {
            format!("partial:{:.2}", (aligned_length / gene_length) * 100.0)
        } else {
            "full".to_string()
        };

        let attributes = format!(
            "target_name={};ID={}_{};coverage={};e_value={:.2e}",
            genid,
            seqid,
            id,
            coverage,
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
