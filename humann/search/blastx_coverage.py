#! /usr/bin/env python

"""
This is a HUMAnN utility function
* Do a first pass on blastx output
* Identify proteins that were well-covered by reads
* Return them as a dict
* When processing blastx output in HUMAnN, consider only these proteins
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""
import math
import sys
import re
import logging
import argparse
from collections import defaultdict

from humann import config
from humann import utilities
from humann import store

# name global logging instance
logger = logging.getLogger(__name__)


def blastx_coverage(blast6out, min_coverage, alignments=None, log_messages=None, apply_filter=None, nucleotide=False,
                    query_coverage_threshold=config.translated_query_coverage_threshold,
                    identity_threshold=config.nucleotide_identity_threshold):
    # create alignments instance if none is passed
    if alignments is None:
        alignments = store.Alignments()

    # store protein lengths
    protein_lengths = {}
    # store unique positions hit in each protein as sets
    protein_hits = defaultdict(set)
    # track proteins with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage = 0
    # parse blast6out file, applying filtering as selected
    for alignment_info in utilities.get_filtered_translated_alignments(blast6out, alignments,
                                                                       apply_filter=apply_filter,
                                                                       log_filter=log_messages,
                                                                       query_coverage_threshold=query_coverage_threshold,
                                                                       identity_threshold=identity_threshold):
        (protein_name, gene_length, queryid, matches, bug, alignment_length,
         subject_start_index, subject_stop_index) = alignment_info

        # divide the gene length by 3 to get protein length from nucleotide length
        if not nucleotide:
            gene_length = gene_length / 3

        # store the protein length
        protein_lengths[protein_name] = gene_length

        # add the range of the alignment to the protein hits
        protein_range = range(subject_start_index, subject_stop_index)
        if protein_range:
            # keep track of unique hit positions in this protein
            protein_hits[protein_name].update(protein_range)
        else:
            no_coverage += 1
    # track proteins without lengths
    no_length = 0
    # compute coverage
    for protein_name, hit_positions in protein_hits.items():
        try:
            # compute coverage, with 50 indicating that 50% of the protein is covered
            coverage = len(hit_positions) / float(protein_lengths[protein_name]) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length += 1

        if coverage >= min_coverage:
            allowed.add(protein_name)

    output_messages = ["Total alignments without coverage information: " + str(no_coverage)]
    output_messages += ["Total proteins in blastx output: " + str(len(protein_lengths))]
    output_messages += ["Total proteins without lengths: " + str(no_length)]
    output_messages += [
        "Proteins with coverage greater than threshold (" + str(min_coverage) + "): " + str(len(allowed))]

    # write out informational messages to log or stdout, depending on input parameters
    if log_messages:
        for message in output_messages:
            logger.info(message)
    else:
        print("\n".join(output_messages))

    return allowed


def blastx_coverage_parallel(alignment_file_tsv, min_coverage, alignments=None,
                             log_messages=None, apply_filter=None, nucleotide=False,
                             query_coverage_threshold=config.translated_query_coverage_threshold,
                             identity_threshold=config.nucleotide_identity_threshold):
    """
    并行计算
    融合原有的blastx_coverage和utilities.get_filtered_translated_alignments
    """
    # create alignments instance if none is passed
    if alignments is None:
        alignments = store.Alignments()

    # store protein lengths
    protein_lengths = {}
    # store unique positions hit in each protein as sets
    protein_hits = defaultdict(set)
    # track proteins with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage = 0
    # parse blast6out file, applying filtering as selected
    # if identity threshold is not set, use the config default
    if identity_threshold is None:
        identity_threshold = config.identity_threshold
    import time
    start = time.time()
    from multiprocessing import Pool
    from itertools import repeat

    log_evalue = False
    large_evalue_count = 0
    small_identity_count = 0
    small_query_coverage_count = 0
    percent_identity_convert_error = 0
    alignment_length_convert_error = 0
    evalue_convert_error = 0
    rapsearch_evalue_convert_error = 0
    with Pool(config.threads) as p, \
            open(alignment_file_tsv, "rt") as file_handle:
        alignments = store.Alignments()
        chunk_size = 1000000 * config.threads
        lines = file_handle.readlines(chunk_size)
        while len(lines) > 0:
            res = p.starmap(line_process, zip(lines,
                                              repeat(alignments),
                                              repeat(identity_threshold),
                                              repeat(query_coverage_threshold),
                                              repeat(apply_filter),
                                              repeat(nucleotide)))

            # no_coverage += sum([r['no_coverage'] for r in filter(None, res)])
            # percent_identity_convert_error += \
            #     sum([r['percent_identity_convert_error'] for r in filter(None, res)])
            # alignment_length_convert_error += \
            #     sum([r['alignment_length_convert_error'] for r in filter(None, res)])
            # large_evalue_count += \
            #     sum([r['large_evalue_count'] for r in filter(None, res)])
            # small_identity_count += \
            #     sum([r['small_identity_count'] for r in filter(None, res)])
            # small_query_coverage_count += \
            #     sum([r['small_query_coverage_count'] for r in filter(None, res)])

            for r in filter(None, res):
                protein_hits[r['protein_name']].update(r['protein_range'])
                protein_lengths[r['protein_name']] = r['gene_length']

            lines = file_handle.readlines(chunk_size)
    # TODO 完善log
    # if log_messages:
    #     logger.debug("Total alignments where percent identity is not a number: " + str(percent_identity_convert_error))
    #     logger.debug("Total alignments where alignment length is not a number: " + str(alignment_length_convert_error))
    #     logger.debug("Total alignments where E-value is not a number: " + str(evalue_convert_error))
    #     if log_evalue:
    #         logger.debug("Total alignments unable to convert rapsearch e-value: " + str(rapsearch_evalue_convert_error))
    #     logger.debug("Total alignments not included based on large e-value: " +
    #                  str(large_evalue_count))
    #     logger.debug("Total alignments not included based on small percent identity: " +
    #                  str(small_identity_count))
    #     logger.debug("Total alignments not included based on small query coverage: " +
    #                  str(small_query_coverage_count))

    # track proteins without lengths
    no_length = 0
    for protein_name, hit_positions in protein_hits.items():
        try:
            # compute coverage, with 50 indicating that 50% of the protein is covered
            coverage = len(hit_positions) / float(protein_lengths[protein_name]) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length += 1

        if coverage >= min_coverage:
            allowed.add(protein_name)

    output_messages = ["Total alignments without coverage information: " + str(no_coverage)]
    output_messages += ["Total proteins in blastx output: " + str(len(protein_lengths))]
    output_messages += ["Total proteins without lengths: " + str(no_length)]
    output_messages += [
        "Proteins with coverage greater than threshold (" + str(min_coverage) + "): " + str(len(allowed))]

    # write out informational messages to log or stdout, depending on input parameters
    if log_messages:
        for message in output_messages:
            logger.info(message)
    else:
        print("\n".join(output_messages))
    return allowed


def line_process(line, alignments, identity_threshold, query_coverage_threshold, apply_filter, nucleotide):
    log_evalue = False
    large_evalue_count = 0
    small_identity_count = 0
    small_query_coverage_count = 0
    percent_identity_convert_error = 0
    alignment_length_convert_error = 0
    evalue_convert_error = 0
    rapsearch_evalue_convert_error = 0

    if line.startswith("#"):
        # Check for the rapsearch2 header to determine if these are log(e-value)
        if re.search(config.blast_delimiter, line):
            data = line.split(config.blast_delimiter)
            if len(data) > config.blast_evalue_index:
                if re.search("log", data[config.blast_evalue_index]):
                    log_evalue = True
    else:
        alignment_info = line.split(config.blast_delimiter)

        # try to obtain the identity value to determine if threshold is met
        identity = alignment_info[config.blast_identity_index]
        try:
            identity = float(identity)
        except ValueError:
            percent_identity_convert_error += 1
            identity = 0.0

        queryid = alignment_info[config.blast_query_index]

        # try converting the alignment length to a number
        alignment_length = alignment_info[config.blast_aligned_length_index]
        try:
            alignment_length = float(alignment_length)
        except ValueError:
            alignment_length_convert_error += 1
            alignment_length = 0.0

        # try converting evalue to float to check if it is a number
        evalue = alignment_info[config.blast_evalue_index]
        try:
            evalue = float(evalue)
        except ValueError:
            evalue_convert_error += 1
            evalue = 1.0

        # try to get the start and end positions for the query
        try:
            query_start_index = int(alignment_info[config.blast_query_start_index])
            query_stop_index = int(alignment_info[config.blast_query_end_index])
        except (ValueError, IndexError):
            query_start_index = 0
            query_stop_index = 0

        # check for query length annotation
        queryid, query_length = utilities.get_length_annotation(queryid)

        # try to get the start and end positions for the subject
        try:
            subject_start_index = int(alignment_info[config.blast_subject_start_index])
            subject_stop_index = int(alignment_info[config.blast_subject_end_index])
        except (ValueError, IndexError):
            subject_start_index = 0
            subject_stop_index = 0

        # convert rapsearch evalue to blastm8 format if logged
        if log_evalue:
            try:
                evalue = math.pow(10.0, evalue)
            except (ValueError, OverflowError):
                rapsearch_evalue_convert_error += 1
                evalue = 1.0

        # compute the number of matches
        matches = identity / 100.0 * alignment_length

        # get the protein alignment information
        protein_name, gene_length, bug = alignments.process_reference_annotation(
            alignment_info[config.blast_reference_index])

        # check if percent identity is less then threshold
        filter = False
        if identity < identity_threshold:
            filter = True
            small_identity_count += 1

        # filter alignments with evalues greater than threshold
        if evalue > config.evalue_threshold:
            filter = True
            large_evalue_count += 1

        # filter alignments that do not meet query coverage threshold
        if utilities.filter_based_on_query_coverage(query_length, query_start_index, query_stop_index,
                                                    query_coverage_threshold):
            filter = True
            small_query_coverage_count += 1
        if apply_filter:
            if not filter:
                alignment_info = (protein_name, gene_length, queryid, matches, bug,
                                  alignment_length, subject_start_index, subject_stop_index)
            # TODO 移除unaligned_reads_store
            # elif unaligned_reads_store:
            #     # remove the read from the unaligned reads store
            #     unaligned_reads_store.remove_id(queryid)
            else:
                return None
        else:
            alignment_info = (protein_name, gene_length, queryid, matches, bug,
                              alignment_length, subject_start_index, subject_stop_index)

        # divide the gene length by 3 to get protein length from nucleotide length
        if not nucleotide:
            gene_length = gene_length / 3
        protein_range = range(subject_start_index, subject_stop_index)
        no_coverage = 0
        if not protein_range:
            protein_range = None
            no_coverage = 1
        return {
            'protein_name': protein_name,
            'gene_length': gene_length,
            'protein_range': protein_range,
            'no_coverage': no_coverage,

            'log_evalue': str(log_evalue),
            'percent_identity_convert_error': percent_identity_convert_error,
            'alignment_length_convert_error': alignment_length_convert_error,
            'evalue_convert_error': evalue_convert_error,
            'rapsearch_evalue_convert_error': rapsearch_evalue_convert_error,
            'large_evalue_count': large_evalue_count,
            'small_identity_count': small_identity_count,
            'small_query_coverage_count': small_query_coverage_count
        }


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="Compute blastx coverage\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the blastx formatted input file\n",
        required=True)
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        help="the subject coverage threshold\n[ DEFAULT : " + str(config.translated_subject_coverage_threshold) + " ]",
        default=config.translated_subject_coverage_threshold)
    parser.add_argument(
        "--print-protein-list",
        action="store_true",
        help="print the list of proteins that meet the coverage threshold")

    return parser.parse_args()


def main():
    # parse the arguments from the user
    args = parse_arguments(sys.argv)

    # run coverage computation
    allowed = blastx_coverage(args.input, args.coverage_threshold)

    if args.print_protein_list:
        print("\n".join(allowed))


if __name__ == "__main__":
    main()
