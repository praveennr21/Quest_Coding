#!/usr/bin/env python
"""
translate.py translates transcript coordinates to genomic coordinates.

Usage: translate.py -t [--tfile] -q [--qfile] -o [--ofolder]
-t [--tfile]: /path/to/transcripts_file
-q [--qfile]: /path/to/queries_file
-o [--ofolder]: /path/to/output_folder
"""

# Built-in/Generic Imports
import sys
import os
import getopt
import re
import logging
from collections import defaultdict

# Third-party Libs
from cigar import Cigar

__author__ = "Praveen Nadukkalam Ravindran"
__version__ = "1.0.0"
__email__ = "praveen.nadukkalamravindran@mcgill.ca"
__status__ = "Development"

logger = logging.getLogger("translate logger")


def configure_logging(log_file_path):
    """
    Configures the translate logger.
    """
    logger.setLevel(logging.DEBUG)
    # Format log lines
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',
                                  datefmt='%m/%d/%Y %I:%M:%S %p')
    # Setup for console logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Setup for file logging
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)


def load_transcripts(transcripts_file_path, transcripts):
    """
    Loads the transcripts from the transcripts_file_path into a dictionary.
    For each transcript load the corresponding chromosome name,
    0-based starting position on the chromosome and
    CIGAR string indicating the mapping.

    :param transcripts_file_path: path/to/transcripts_file
    """

    line_num = 0
    # Read each transcript and load the corresponding chromosome name,
    # 0-based starting position on the chromosome and
    # CIGAR string indicating the mapping.
    with open(transcripts_file_path) as transcripts_file:
        for transcript in transcripts_file:
            try:
                line_num += 1
                words = transcript.split('\t')
                fields = {}
                fields['CHR'] = words[1]
                fields['POS'] = int(words[2])
                fields['CIGAR'] = validate_cigar(words[3].rstrip())
                transcripts[words[0]] = fields
            # Handle exception in case of missing values for transcript.
            except ValueError:
                logger.error('Error reading transcripts file: %s at line %d',
                             transcripts_file_path, line_num)
                sys.exit(2)
    logger.info('Loaded %d transcripts.', len(transcripts))


def validate_cigar(cigar):
    """
    Validates the CIGAR string for the mapping.
    Checks if the CIGAR string complies to standard CIGAR format i.e,
    Operators: D,I,M
    Format: Natural number followed by operator
    Example: 8M7D6M2I2M11D7M

    :param cigar: CIGAR string for the mapping
    """
    line_num = 0
    try:
        line_num += 1
        # Check if the CIGAR string matches the validation regex.
        if not re.match("^([1-9]+[0-9]*[MDI])+$", cigar):
            raise Exception()
        out_cigar = Cigar(cigar)
        out_cigar_list = list(out_cigar.items())
    # Handle exception in case of invalid CIGAR string.
    except Exception:
        logger.error('Error reading CIGAR string : %s at line %s', cigar,
                     line_num)
        sys.exit(1)
    return out_cigar_list


def process_queries(transcripts, queries_file_path, output_file_path):
    """
    Processes the queries from the queries_file_path.
    Generates and writes output to the output_file_path.
    For each query, if found, outputs the corresponding chromosome name and
    the genomic coordinate, otherwise, outputs '-' for the chromosome name
    and the genomic coordinate fields.
    In case of insertion the start position of insertion
    is output as genomic coordinate.

    :param queries_file_path: path/to/queries_file
    :param output_file_path: path/to/output_file
    """

    line_num = 0
    with open(queries_file_path) as queries_file:
        for query in queries_file:
            transPos = -1
            process_queries.outPos = ''
            genomePos = 0
            cigar = []
            op = ''
            tr = ''
            process_queries.chr = ''
            try:
                line_num += 1
                words = query.split('\t')
                tr = words[0]
                queryPos = int(words[1].rstrip())
            # Handle exception in case of missing values in queries file.
            except Exception:
                logger.error('Error reading queries file: %s at line %s',
                             queries_file_path, line_num)
                sys.exit(1)
            # Check if the query transcript present in transcripts.
            if (tr in transcripts):
                cigar = transcripts[tr]['CIGAR']
                genomePos = transcripts[tr]['POS']
                process_queries.chr = transcripts[tr]['CHR']
                for item in cigar:
                    if (item[1] == 'M'):
                        transPos += item[0]
                        genomePos += item[0]
                    elif (item[1] == 'D'):
                        genomePos += item[0]
                    elif (item[1] == 'I'):
                        transPos += item[0]
                    op = item[1]
                    if (transPos >= queryPos):
                        break
                # Check if the query position for the transcript not found.
                # if not found set chr and outPos to '-'.
                if ((transPos - queryPos) < 0):
                    process_queries.chr = '-'
                    process_queries.outPos = '-'
                else:
                    # If query postion is an insertion set outPos to the start
                    # position of the insertion
                    if (op == 'I'):
                        process_queries.outPos = genomePos - 1
                    else:
                        process_queries.outPos = genomePos - (transPos -
                                                              queryPos + 1)
            # If the query transcript not present in transcripts
            # set chr and outPos to '-'.
            else:
                process_queries.chr = '-'
                process_queries.outPos = '-'
            # Write output to the output_file_path
            if (line_num == 1):
                with open(output_file_path, "w") as out:
                    out.write(tr + "\t" + str(queryPos) + "\t" +
                              process_queries.chr + "\t" +
                              str(process_queries.outPos) + "\n")
            else:
                with open(output_file_path, "a") as out:
                    out.write(tr + "\t" + str(queryPos) + "\t" +
                              process_queries.chr + "\t" +
                              str(process_queries.outPos) + "\n")
    logger.info('Successfully processed %d queries.', line_num)


def main(argv):
    """
    Parse input arguments.
    Display usage message and error message in case of missing arguments.
    Call the load_transcripts() and process_queries() functions

    :param argv: Input arguments.
    """
    transcripts_file_path = ''
    queries_file_path = ''
    output_file_path = ''
    log_file_path = ''
    transcripts = defaultdict(dict)

    try:
        opts, args = getopt.getopt(argv, "ht:q:o:",
                                   ["help", "tfile=", "qfile=", "ofolder="])
    except getopt.GetoptError:
        logger.error(__doc__)
        sys.exit(2)
    if not opts:
        logger.error(__doc__)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(__doc__)
            sys.exit()
        elif opt in ("-t", "--tfile"):
            # Exit if the transcripts file doesn't exist.
            if not os.path.isfile(arg):
                logger.error(
                    'Transcripts file path %s does not exist. Exiting...', arg)
                sys.exit(2)
            transcripts_file_path = arg
        elif opt in ("-q", "--qfile"):
            # Exit if the queries file doesn't exist.
            if not os.path.isfile(arg):
                logger.error('Queries file path %s does not exist. Exiting...',
                             arg)
                sys.exit(2)
            queries_file_path = arg
        elif opt in ("-o", "--ofolder"):
            # Exit if the output folder doesn't exist.
            if not os.path.isdir(arg):
                logger.error(
                    'Output folder path %s does not exist. Exiting...', arg)
                sys.exit(2)
            output_file_path = arg + '/output.txt'
            log_file_path = arg + '/translate.log'
    if (transcripts_file_path == ''):
        logger.error('Provide -t or --tfile')
        sys.exit(2)
    if (queries_file_path == ''):
        logger.error('Provide -q or --qfile')
        sys.exit(2)
    if (output_file_path == ''):
        logger.error('Provide -o or --ofolder')
        sys.exit(2)

    configure_logging(log_file_path)
    logger.info('Started processing the queries.')
    load_transcripts(transcripts_file_path, transcripts)
    process_queries(transcripts, queries_file_path, output_file_path)
    logger.info('Finished processing the queries.')


if __name__ == '__main__':
    main(sys.argv[1:])
