# -*- coding: utf-8 -*-
"""cpo-pipeline.assembly.parsers.result_parsers

This module provides functions for parsing result files generated by tools
during the QC & Assembly phase of the cpo-pipeline.
"""

import csv
import pandas

def parse_kraken_result(path_to_kraken_result):
    """
    Args:
        path_to_kraken_result (str): Path to the kraken report file.

    Returns:
        dict: Parsed kraken report with species-level results.
        For example:
        { "Escherichia coli": { "fragment_percent": 84.08,
                                "fragment_count_root": 195536,
                                "fragment_count_taxon": 192561,
                                "rank_code": "S",
                                "ncbi_taxon_id": "562",
                                "name": "Escherichia coli",
                              }
          "Another species": { "fragment_percent": 12.1,
                               ...
                             }
        }
        See kraken manual for more detail on report fields:
        http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-report-output-format
    """
    kraken_result = {}

    # Associate a field name with a function for parsing that field
    kraken_report_fields = {
        'fragment_percent': lambda x: float(x.strip()),
        'fragment_count_root': lambda x: int(x),
        'fragment_count_taxon': lambda x: int(x),
        'rank_code': lambda x: x,
        'ncbi_taxon_id': lambda x: x,
        'name': lambda x: x.strip()
    }
    
    with open(path_to_kraken_result) as krakenfile:
        reader = csv.DictReader(krakenfile, delimiter='\t', fieldnames=list(kraken_report_fields))
        kraken_species_record = {}
        for row in reader:
            # Only parse species-level records (not Family, Genus, etc.)
            if row['rank_code'] != 'S':
                continue
            for field_name, parse in kraken_report_fields.items():
                kraken_species_record[field_name] = parse(row[field_name])
            kraken_result[kraken_species_record['name']] = kraken_species_record
    return kraken_result

def parse_fastqc_result(path_to_qc_summary):
    """
    Args:
        path_to_qc_summary (str): Path to the fastqc report summary file.

    Returns:
        dict: Parsed fastqc R1 report.
              All values except are either "PASS", "WARN", or "FAIL".
        For example:
        { "basic_statistics": "PASS",
          "per_base_sequence_quality": "PASS",
          "per_tile_sequence_quality": "PASS",
          "per_sequence_quality_scores": "PASS",
          "per_base_sequence_content": "WARN",
          "per_sequence_gc_content": "FAIL",
          "per_base_n_content": "PASS",
          "sequence_length_distribution": "WARN",
          "sequence_duplication_levels": "PASS",
          "overrepresented_sequences": "WARN",
          "adapter_content": "PASS",
        }
    """
    fastqc_summary = {}
    with open(path_to_qc_summary) as qc_summary:
        reader = csv.reader(qc_summary, delimiter='\t')
        for row in reader:
            field_name = row[1].lower().replace(" ", "_")
            fastqc_summary[field_name] = row[0]
    return fastqc_summary

def parse_mash_result(path_to_mash_screen):
    """
    Args:
        path_to_mash_screen (str): Path to the mash screen report file.

    Returns:
        dict: Parsed mash screen report
        For example:
        { "Citrobacter freundii strain CAV1321": { "identity": 0.996805,
                                                   "shared_hashes": "935/1000",
                                                   "median_multiplicity": 38,
                                                   "p_value": 0.00,
                                                   "query_id": "GCF_001022155.1_ASM102215v1_genomic.fna.gz",
                                                   "query_comment": "[10 seqs] NZ_CP011612.1 Citrobacter freundii strain CAV1321, complete genome [...]"
                                                 },
          "Another species": { "identity": 0.914483,
                               ...
                             }
        }
        See mash docs for more info on mash screen report file:
        https://mash.readthedocs.io/en/latest/tutorials.html#screening-a-read-set-for-containment-of-refseq-genomes
        boolean: true if phiX is present in top hits, false if absent
    
    """
    mash_result = {}

    mash_screen_report_fields = {
        'identity': lambda x: float(x),
        'shared_hashes': lambda x: x,
        'median_multiplicity': lambda x: int(x),
        'p_value': lambda x: float(x),
        'query_id': lambda x: x,
        'query_comment': lambda x: x
    }
    
    # Example mash screen report record (actual report has no header and is tab-delimited):
    # identity    shared-hashes    median-multiplicity    p-value    query-ID                                    query-comment
    # 0.998697    973/1000         71                     0          GCF_000958965.1_matepair4_genomic.fna.gz    [59 seqs] NZ_LAFU01000001.1 Klebsiella pneumoniae strain CDPH5262 contig000001, whole genome shotgun sequence [...]

    with open(path_to_mash_screen) as mashfile:
        reader = csv.DictReader(mashfile, delimiter='\t', fieldnames=list(mash_screen_report_fields))
        mash_record = {}
        for row in reader:
            for field_name, parse in mash_screen_report_fields.items():
                mash_record[field_name] = parse(row[field_name])
            species_name_start = int(mash_record['query_comment'].index(".")) + 3
            species_name_stop = int(mash_record['query_comment'].index(","))
            species = str(mash_record['query_comment'])[species_name_start: species_name_stop]    
            mash_result[species] = mash_record
    return mash_result

def parse_read_stats(path_to_mash_log, path_to_total_bp):
    """
    Args:
        path_to_mash_log (str): Path to the mash log file.
        path_to_total_bp (str): Path to

    Returns:
        dict: Read statistics
        For example:
          { "size": 5185840,
            "depth": 22.28
          }
    """
    total_bp = int([line.rstrip('\n') for line in open(path_to_total_bp)][0])
    mash_log = [line.rstrip('\n') for line in open(path_to_mash_log)]
    read_stats = {}
    for line in mash_log:
        if (line.find("Estimated genome size:") > -1 ):
            size = float(line[line.index(": ") + 2:])
    depth = total_bp / size
    depth = float(format(depth, '.2f'))
    read_stats['size'] = size
    read_stats['depth'] = depth
    return read_stats

def parse_busco_result(path_to_busco_result):
    """
    Args:
        path_to_busco_log (str): Path to the mash log file.
    Returns:
        dict: Parsed busco report
        For example:
        { "complete_single": 778,
          "complete_duplicate": 1,
          "fragmented": 1,
          "missing": 1,
        }
    Raises:
        Exception: When BUSCO output file is empty
    """
    busco_output = [line.rstrip('\n') for line in open(path_to_busco_result)]
    if (len(busco_output) > 0):
        busco_result = {}
        busco_result['complete_single'] = int(busco_output[10].split("\t")[1])
        busco_result['complete_duplicate'] = int(busco_output[11].split("\t")[1])
        busco_result['fragmented'] = int(busco_output[12].split("\t")[1])
        busco_result['missing'] = int(busco_output[13].split("\t")[1])
        busco_result['total'] = int(busco_output[14].split("\t")[1])
    else:
        raise Exception("BUSCO output file has no contents.")
    return busco_result

def parse_quast_result(path_to_quast_result):
    """
    Args:
        path_to_quast_result (str): Path to the QUAST result file.
    Returns:
        dict: Parsed QUAST report
        For example:
        { "contigs_count": 163,
          "largest_contig": 293450,
          "N50": 111197,
          "NG50": 112240,
          "L50": 16,
          "LG50": 15,
          "N75": 57407,
          "NG75": 65034,
          "L75": 32,
          "LG75": 29,
          "total_length": 5432369,
          "reference_length": 5169997,
          "percent_GC": 50.55,
          "reference_percent_GC": 50.65,
          "genome_fraction_percent": 96.296,
          "duplication_ratio": 1.01,
          "data": ["Assembly	BC18-Eco016","# contigs (>= 0 bp)	163", ...]
        }
    """
    quast_output = [line.rstrip('\n') for line in open(path_to_quast_result)]
    quast_result = {}
    quast_result['contigs_count'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("# contigs\t") > -1][-1])].split("\t")[1])
    quast_result['largest_contig'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Largest contig\t") > -1][-1])].split("\t")[1])
    quast_result['N50'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("N50\t") > -1][-1])].split("\t")[1])
    quast_result['NG50'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("NG50\t") > -1][-1])].split("\t")[1])
    quast_result['L50'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("L50\t") > -1][-1])].split("\t")[1])
    quast_result['LG50'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("LG50\t") > -1][-1])].split("\t")[1])
    quast_result['N75'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("N75\t") > -1][-1])].split("\t")[1])
    quast_result['NG75'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("NG75\t") > -1][-1])].split("\t")[1])
    quast_result['L75'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("L75\t") > -1][-1])].split("\t")[1])
    quast_result['LG75'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("LG75\t") > -1][-1])].split("\t")[1])
    quast_result['total_length'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Total length\t") > -1][-1])].split("\t")[1])
    quast_result['reference_length'] = int(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Reference length\t") > -1][-1])].split("\t")[1])
    quast_result['percent_GC'] = float(quast_output[int([i for i,x in enumerate(quast_output) if x.find("GC (%)\t") > -1][-1])].split("\t")[1])
    quast_result['reference_percent_GC'] = float(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Reference GC (%)\t") > -1][-1])].split("\t")[1])
    quast_result['genome_fraction_percent'] = float(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Genome fraction (%)\t") > -1][-1])].split("\t")[1])
    quast_result['duplication_ratio'] = float(quast_output[int([i for i,x in enumerate(quast_output) if x.find("Duplication ratio\t") > -1][-1])].split("\t")[1])
    quast_result['data'] = quast_output

    return quast_result

