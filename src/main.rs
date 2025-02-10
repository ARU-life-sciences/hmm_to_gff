// Format the output of hmmscan to a GFF file
// Specifically we are interested in the output from ORFfinder
// fasta -> ORFfinder -> hmmscan -> GFF
// We use the Pfam HMM database to achieve this.
//
// #
// # Program:         hmmscan
// # Version:         3.4 (Aug 2023)
// # Pipeline mode:   SCAN
// # Query file:      ../data/Acaena_ovalifolia_orfs.fasta
// # Target file:     ../pfam_hmm/Pfam-A.hmm
// # Options:         /software/team301/hmmer-3.4/src/hmmscan --tblout ../data/pfam_hmmscans/Acaena_ovalifolia.tbl --cpu 16 ../pfam_hmm/Pfam-A.hmm ../data/Acaena_ovalifolia_orfs.fasta
// # Current dir:     /lustre/scratch123/tol/teams/blaxter/users/mb39/ARU/mito_structural_variation/annotation/src
// # Date:            Fri Feb  7 19 13 12 2025
// # [ok]

use calm_io::*;
use hmm_tblout::Reader;

// example output should be something like:
// u15	ORFfinder	CDS	1350	1589	.	-	.	ID=ORF1766_u15_1350_1589

// program used is ORFfinder
const SOURCE: &str = "ORFfinder";
// feature type is CDS
const FEATURE: &str = "CDS";
// frames are not available
const FRAME: &str = ".";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // we want to read in the HMMER tblout file and convert it to a GFF file
    // this will allow us to use tools like bedtools to compare it to other annotations
    // and to visualize it in a genome browser

    // load in the command line arguments
    let args: Vec<String> = std::env::args().collect();
    // we only expect one arg, the path to the HMM file
    if args.len() != 2 {
        stderrln!("Usage: {} <HMM tblout file>", args[0])?;
        std::process::exit(1);
    }

    let file_path = &args[1];

    // now create a HMMER tblout reader
    let mut reader = Reader::from_path(file_path)?;

    let mut gff_records = Vec::new();

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

        let score = record.e_value_full().unwrap();

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
        // query_name_full=""; accession=""; target_name=""; description="";
        let attr = format!(
            "query_name_full=\"{}\";pfam_accession=\"{}\";target_name=\"{}\";description=\"{}\"",
            record.query_name(),
            record.target_accession(),
            record.target_name(),
            record.description()
        );

        gff_records.push((
            seqname,
            SOURCE,
            FEATURE,
            final_start,
            final_end,
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

    stdoutln!("##gff-version 3")?;

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
        stdoutln!("{}", gff_line)?;
    }

    Ok(())
}
