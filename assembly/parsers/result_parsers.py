# -*- coding: utf-8 -*-
"""cpo-pipeline.assembly.parsers.result_parsers

This module provides functions for parsing result files generated by tools
during the QC & Assembly phase of the cpo-pipeline.
"""

import csv
import re

def parse_kraken_result(path_to_kraken_result):
    """
    Args:
        path_to_kraken_result (str): Path to the kraken report file.

    Returns:
        list(dict): Parsed kraken report with species-level results.
        For example:
        [
            { "fragment_percent": 84.08,
              "fragment_count_root": 195536,
              "fragment_count_taxon": 192561,
              "rank_code": "S",
              "ncbi_taxon_id": "562",
              "taxon_name": "Escherichia coli",
            },
            { "fragment_percent": 12.1,
              ...
            }
        ]
        See kraken manual for more detail on report fields:
        http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-report-output-format
    """

    # Associate a field name with a function for parsing that field
    kraken_report_fields = {
        'fragment_percent': lambda x: float(x.strip()),
        'fragment_count_root': lambda x: int(x),
        'fragment_count_taxon': lambda x: int(x),
        'rank_code': lambda x: x,
        'ncbi_taxon_id': lambda x: x,
        'taxon_name': lambda x: x.strip()
    }
    
    with open(path_to_kraken_result) as krakenfile:
        reader = csv.DictReader(krakenfile, delimiter='\t', fieldnames=list(kraken_report_fields))
        parsed_kraken_result = []
        kraken_record = {}
        for row in reader:
            for field_name, parse in kraken_report_fields.items():
                kraken_record[field_name] = parse(row[field_name])
            parsed_kraken_result.append(kraken_record.copy())
    return parsed_kraken_result

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
        list(dict): Parsed mash screen report
        For example:
        [
          { "identity": 0.996805,
            "shared_hashes": "935/1000",
            "median_multiplicity": 38,
            "p_value": 0.00,
            "query_id": "GCF_001022155.1_ASM102215v1_genomic.fna.gz",
            "query_comment": "[10 seqs] NZ_CP011612.1 Citrobacter freundii strain CAV1321, complete genome [...]"
          },
          { "identity": 0.914483,
            ...
          }
        ]
        See mash docs for more info on mash screen report file:
        https://mash.readthedocs.io/en/latest/tutorials.html#screening-a-read-set-for-containment-of-refseq-genomes
    """

    mash_screen_report_fields = {
        'identity': lambda x: float(x),
        'shared_hashes': lambda x: x,
        'median_multiplicity': lambda x: int(x),
        'p_value': lambda x: float(x),
        'query_id': lambda x: x,
        'query_comment': lambda x: x
    }
    
    # Example mash screen report record (actual report has no header and is tab-delimited):
    # identity    shared_hashes    median_multiplicity    p_value    query_id                                    query_comment
    # 0.998697    973/1000         71                     0          GCF_000958965.1_matepair4_genomic.fna.gz    [59 seqs] NZ_LAFU01000001.1 Klebsiella pneumoniae strain CDPH5262 contig000001, whole genome shotgun sequence [...]

    parsed_mash_result = []
    with open(path_to_mash_screen) as mashfile:
        reader = csv.DictReader(mashfile, delimiter='\t', fieldnames=list(mash_screen_report_fields))
        mash_record = {}
        for row in reader:
            for field_name, parse in mash_screen_report_fields.items():
                mash_record[field_name] = parse(row[field_name])
            parsed_mash_result.append(mash_record.copy())
    return parsed_mash_result

def parse_read_stats(path_to_mash_log, path_to_total_bp):
    """
    Args:
        path_to_mash_log (str): Path to the mash log file.
        path_to_total_bp (str): Path to

    Returns:
        dict: Read statistics
        For example:
          { "estimated_genome_size": 5185840,
            "estimated_depth_of_coverage": 22.28
          }
    """
    with open(path_to_total_bp, 'r') as totalbp_file:
        total_bp = int(totalbp_file.readline())

    with open(path_to_mash_log, 'r') as mash_log_file:
        for line in mash_log_file:
            if (line.find("Estimated genome size:") > -1 ):
                estimated_genome_size = float(line.split(":")[1].strip())
    
    read_stats = {}
    estimated_depth_of_coverage = total_bp / estimated_genome_size
    read_stats['estimated_genome_size'] = int('%.0f' % estimated_genome_size)
    read_stats['estimated_depth_of_coverage'] = float('%.2f' % estimated_depth_of_coverage)
    return read_stats
	
def parse_total_bp(path_to_total_bp):
    """
    Args:
        path_to_total_bp (str): Path to

    Returns:
        float: total basepairs
        For example:
          5185840
    """
    with open(path_to_total_bp, 'r') as totalbp_file:
        total_bp = float(totalbp_file.readline())

    return total_bp
    
def parse_reference_genome_stats(path_to_reference_genome_stats):
    """
    Args:
        path_to_reference_genome_stats (str): Path to

    Returns:
        float: genome size
        For example:
          5185840
    """
    with open(path_to_reference_genome_stats, 'r') as reference_stats:
        genome_stats = reference_stats.read().splitlines()
    for line in genome_stats:
        if (line.find("all	all	all	all	total-length") > -1): #find the total length stat
            expected_genome_size = float(line.split("\t")[5].strip()) 
    return expected_genome_size

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
    busco_output = []
    with open(path_to_busco_result, 'r') as busco_result:
        for line in busco_result:
            busco_output.append(line)
    
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
        {
            "contigs_count": 72,
            "largest_contig": 692871,
            "N50": 299446,
            "NG50": 299446,
            "L50": 6,
            "LG50": 6,
            "N75": 123167,
            "NG75": 110534,
            "L75": 12,
            "LG75": 14,
            "total_length": 5182695,
            "reference_length": 5489397,
            "percent_GC": 51.75,
            "reference_percent_GC": 51.59,
            "genome_fraction_percent": 91.202,
            "duplication_ratio": 1.002
        }
    """
    quast_output = []
    with open(path_to_quast_result, 'r') as quast_result:
        for line in quast_result:
            quast_output.append(line)

    def parse_quast_report_line(line):
        """
        Takes a line of the quast report and returns the specific data that we're interested in from that line.
        
        Collapse multiple spaces into a tab char ('\t'), then split the line on tabs and take the second item.
        Cast floats to floats and ints to ints.
        '# contigs        751   ' -> '# contigs\t751\t' -> ['# contigs', '751', ''] -> '751' -> 751
        """
        result_data = re.sub(' {2,}', '\t', line).split('\t')[1]
        if re.match('\d+\.\d+', result_data):
            return float(result_data)
        else:
            return int(result_data)

    # Associate a regex that can be used to identify the line of interest in the quast report
    # with a string to use as a key in the output dict.
    quast_report_parsing_regexes = {
        '^# contigs {2,}\d+': 'contigs_count',
        '^Largest contig {1,}\d+': 'largest_contig',
        '^N50 +\d+': 'N50',
        '^NG50 +\d+': 'NG50',
        '^L50 +\d+': 'L50',
        '^LG50 +\d+': 'LG50',
        '^N75 +\d+': 'N75',
        '^NG75 +\d+': 'NG75',
        '^L75 +\d+': 'L75',
        '^LG75 +\d+': 'LG75',
        '^Total length {2,}\d+': 'total_length',
        '^Reference length +\d+': 'reference_length',
        '^GC \(%\) +\d+\.\d+': 'percent_GC',
        '^Reference GC \(%\) +\d+\.\d+': 'reference_percent_GC',
        '^Genome fraction \(%\) +\d+\.\d+': 'genome_fraction_percent',
        '^Duplication ratio +\d+\.\d+': 'duplication_ratio',
    }
    
    quast_result = {}
    for line in quast_output:
        for regex, key in quast_report_parsing_regexes.items():
            if re.match(regex, line):
                quast_result[key] = parse_quast_report_line(line)
    
    return quast_result

