#!/usr/bin/env python

"""
test_translate.py performs automated unit testing for translate.py.

"""
# Built-in/Generic Imports
import pytest
import sys
from collections import defaultdict

# Set path
sys.path.append('../translate')

# Import functions from translate.py
from translate.translate import *

__author__ = "Praveen Nadukkalam Ravindran"
__version__ = "1.0.0"
__email__ = "praveen.nadukkalamravindran@mcgill.ca"
__status__ = "Development"


def test_transcripts_file_corrupt():
    transcripts_file_path = '../data/transcripts_bad.txt'
    transcripts = defaultdict(dict)
    with pytest.raises(SystemExit):
        load_transcripts(transcripts_file_path, transcripts)


def test_transcripts_size():
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    assert len(transcripts) == 2


def test_validate_cigar():
    cigar = '8M2I4D'
    out_cigar_list = validate_cigar(cigar)
    assert len(out_cigar_list) == 3


def test_validate_cigar1():
    cigar = '-M'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_validate_cigar2():
    cigar = '0M'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_validate_cigar3():
    cigar = 'M8'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_validate_cigar4():
    cigar = '8MI'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_validate_cigar5():
    cigar = '8M9II'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_validate_cigar6():
    cigar = '8M9S'
    with pytest.raises(SystemExit):
        validate_cigar(cigar)


def test_queries_file_corrupt():
    queries_file_path = '../data/queries_bad.txt'
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    output_file_path = '../data/output/output.txt'
    log_file_path = '../data/output/translate.log'
    configure_logging(log_file_path)
    with pytest.raises(SystemExit):
        process_queries(transcripts, queries_file_path, output_file_path)


def test_queries_valid_query():
    queries_file_path = '../data/query_1.txt'
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    output_file_path = '../data/output/output.txt'
    log_file_path = '../data/output/translate.log'
    configure_logging(log_file_path)
    process_queries(transcripts, queries_file_path, output_file_path)
    assert process_queries.chr == 'CHR1'
    assert process_queries.outPos == 7


def test_queries_insertion_query():
    queries_file_path = '../data/query_2.txt'
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    output_file_path = '../data/output/output.txt'
    log_file_path = '../data/output/translate.log'
    configure_logging(log_file_path)
    process_queries(transcripts, queries_file_path, output_file_path)
    assert process_queries.chr == 'CHR1'
    assert process_queries.outPos == 23


def test_queries_genome_coordinate_not_found():
    queries_file_path = '../data/query_3.txt'
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    output_file_path = '../data/output/output.txt'
    log_file_path = '../data/output/translate.log'
    configure_logging(log_file_path)
    process_queries(transcripts, queries_file_path, output_file_path)
    assert process_queries.chr == '-'
    assert process_queries.outPos == '-'


def test_queries_transcript_not_found():
    queries_file_path = '../data/query_4.txt'
    transcripts_file_path = '../data/transcripts.txt'
    transcripts = defaultdict(dict)
    load_transcripts(transcripts_file_path, transcripts)
    output_file_path = '../data/output/output.txt'
    log_file_path = '../data/output/translate.log'
    configure_logging(log_file_path)
    process_queries(transcripts, queries_file_path, output_file_path)
    assert process_queries.chr == '-'
    assert process_queries.outPos == '-'
