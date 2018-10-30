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

def parse_mash_genome_result(path_to_mash_screen, size, depth):
    """
    Args:
        path_to_mash_screen (str): Path to the mash screen report file.
        size (int):
        depth (int):

    Returns:
        tuple(dict, boolean) where:
        dict: Parsed mash screen report
        For example:
        { "Escherichia coli": { "size": ,
                                "depth": ,
                                "identity": 0.998694,
                                "shared_hashes": "972/1000",
                                "median_multiplicity": 16,
                                "p_value": 0.00,
                                "query_ID": "GCF_001519515.1_ASM151951v1_genomic.fna.gz",
                                "query_comment": "[219 seqs] NZ_LQSN01000001.1 Escherichia coli strain GN02215 GCID_ECOLID_00024_NODE_1.ctg_1, whole genome shotgun sequence [...]",
                                "accession": "[219 seqs] NZ_LQSN01000001.1 Escherichia coli strain GN02215 GCID_ECOLID_00024_NODE_1.ctg_1, whole genome shotgun sequence [...]",
                                "row": "0.998649\t972/1000\t16\t0\tGCF_001519515.1_ASM151951v1_genomic.fna.gz\t[219 seqs] NZ_LQSN01000001.1 Escherichia coli strain GN02215 GCID_ECOLID_00024_NODE_1.ctg_1, whole genome shotgun sequence [...]",
                                "gcf": ["GCF", "001", "519", "515"],
                                "assembly": "GCF_001519515.1_ASM151951v1",
                                "species": "Escherichia coli"
                              }
          "Another species": { "size": ,
                               ...
                             }
        }
        See mash docs for more info on mash screen report file:
        https://mash.readthedocs.io/en/latest/tutorials.html#screening-a-read-set-for-containment-of-refseq-genomes
        boolean: true if phiX is present in top hits, false if absent
    
    """
    mash_hits = {}
    phiX = False
    
    mash_screen_report = pandas.read_csv(path_to_mash_screen, delimiter='\t', header=None)
    # Example mash screen report record (actual report has no header and is tab-delimited):
    # identity    shared-hashes    median-multiplicity    p-value    query-ID                                    query-comment]
    # 0.998697    973/1000         71                     0          GCF_000958965.1_matepair4_genomic.fna.gz    [59 seqs] NZ_LAFU01000001.1 Klebsiella pneumoniae strain CDPH5262 contig000001, whole genome shotgun sequence [...]
    # parse mash result, using winner takes all
    scores = mash_screen_report[1].values
    score_cutoff = int(scores[0][:scores[0].index("/")]) - 300 #find the cut off value to use for filtering the report (300 below max)
    index = 0
    #find hits with score within top 300
    for score in scores:
        parsed_score = int(score[:score.index("/")])
        if parsed_score >= score_cutoff:
            index+=1
        else:
            break

    #parse what the species are.
    for i in range(index):
        mash_record = {}
        mash_record['size'] = size
        mash_record['depth'] = depth
        mash_record['identity'] = float(mash_screen_report.ix[i, 0])
        mash_record['shared_hashes'] = mash_screen_report.ix[i, 1]
        mash_record['median_multiplicity'] = int(mash_screen_report.ix[i, 2])
        mash_record['p_value'] = float(mash_screen_report.ix[i, 3])
        mash_record['query_ID'] = mash_screen_report.ix[i, 4]
        mash_record['query_comment'] = mash_screen_report.ix[i, 5]
        mash_record['accession'] = mash_record['query_comment']
        mash_record['row'] = "\t".join(str(x) for x in mash_screen_report.ix[i].tolist())
        qID = mash_record['query_ID']
        # find gcf accession
        gcf = (qID[:qID.find("_",5)]).replace("_","")
        gcf = [gcf[i:i+3] for i in range(0, len(gcf), 3)]
        mash_record['gcf'] = gcf 
        # find assembly name
        mash_record['assembly'] = qID[:qID.find("_genomic.fna.gz")]

        if (mash_record['query_comment'].find("phiX") > -1): #theres phix in top hits, just ignore
            phiX = True
            mash_record['species'] = "PhiX"
        else: #add the non-phix hits to a list
            species_name_start = int(mash_record['query_comment'].index(".")) + 3
            species_name_stop = int (mash_record['query_comment'].index(","))
            mash_record['species'] = str(mash_record['query_comment'])[species_name_start: species_name_stop]
            mash_hits[mash_record['species']] = mash_record
    return mash_hits, phiX

def parse_mash_plasmid_result(path_to_mash_screen, size, depth):
    """
    Args:
        path_to_mash_screen (str): Path to the mash screen report file.
        size (int):
        depth (int):

    Returns:
        dict: Parsed mash screen report
        For example:
        { "NZ_CP010882.1": { "size": ,
                             "depth": ,
                             "identity": 0.992291,
                             "shared_hashes": "850/1000",
                             "median_multiplicity": 19,
                             "p_value": 0.00,
                             "query_ID": "ref|NZ_CP010882.1|",
                             "query_comment": "Escherichia coli strain MNCRE44 plasmid pMNCRE44_6, complete sequence",
                             "accession": "NZ_CP010882.1",
                             "row": "0.992291\t850/1000\t19\t0\tref|NZ_CP010882.1|\tEscherichia coli strain MNCRE44 plasmid pMNCRE44_6, complete sequence",
                             "species": ""
                           }
          "Another accession": { "size": ,
                                 ...
                               }
        }

        See mash docs for more info on mash screen report file:
        https://mash.readthedocs.io/en/latest/tutorials.html#screening-a-read-set-for-containment-of-refseq-genomes
    """
    mash_screen_report = pandas.read_csv(path_to_mash_screen, delimiter='\t', header=None)

    #parse mash result, using winner takes all
    scores = mash_screen_report[1].values
    score_cutoff = int(scores[0][:scores[0].index("/")]) - 100
    index = 0

    #find hits with score > (max_score - 100)
    for score in scores:
        parsed_score = int(score[:score.index("/")])
        if parsed_score >= score_cutoff:
            index+=1
        else:
            break

    mash_plasmid_hits = {}
    #parse what the species are.
    for i in range(index):
        mash_record = {}
        mash_record['size'] = size
        mash_record['depth'] = depth
        mash_record['identity'] = float(mash_screen_report.ix[i, 0])
        mash_record['shared_hashes'] = mash_screen_report.ix[i, 1]
        mash_record['median_multiplicity'] = int(mash_screen_report.ix[i, 2])
        mash_record['p_value'] = float(mash_screen_report.ix[i, 3])
        mash_record['query_ID'] = mash_screen_report.ix[i, 4] #accession
        mash_record['query_comment'] = mash_screen_report.ix[i, 5]
        mash_record['accession'] = mash_record['query_ID'][4:len(mash_record['query_ID'])-1]
        mash_record['row'] = "\t".join(str(x) for x in mash_screen_report.ix[i].tolist())
        mash_record['species'] = ""
        if (mash_record['identity'] >= 0.97):
            mash_plasmid_hits[mash_record['accession']] = mash_record
    return mash_plasmid_hits

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

