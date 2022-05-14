# translate.py
The python module translate.py translates transcript coordinates to genomic coordinates. The program uses the CIGAR string that indicates the transcript-to-genome mapping to translate the transcript coordinate to a genomic coordinate. 

 
## Usage: 

translate.py -h [--help] -t [--tfile] -q [--qfile] -o [--ofolder] 
-t [--tfile]: /path/to/transcripts_file 
-q [--qfile]: /path/to/queries_file 
-o [--ofolder]: /path/to/output_folder 

## Assumptions: 

The CIGAR string is in standard format, i.e., has only M (Match), I (Insertion) and D (Deletion) operators and does not support extended operators like H, N, P, S, X and =. 

The transcripts file has one mapping entry per transcript. 

## Output: 

For each query (transcript name and query coordinate) the program outputs the corresponding chromosome name and chromosome coordinate (if found). 

If a corresponding chromosome name and chromosome coordinate not found for a query transcript and coordinate then the program outputs ‘-‘ or chromosome name and chromosome coordinate. 

In case of insertion in transcript, for any transcript coordinate that is part of the insertion, the program outputs the chromosome coordinate at the start position of insertion. For example, for the CIGAR mapping string 8M2I with starting position 3, the chromosome coordinates for the transcript coordinates 8 and 9 would be 10. 


## File structure: 

* Translate
   * bin
     * translate
   * translate
     * translate.py
   * data
   * tests
     * test_translate.py
 
## Run the program: 

### Using the executable: 
$ /Translate/bin/translate -t <transcripts_file_path> -q <queries_file_path> -o <output_folder_path> 

### Using the .py file: 

Requires python version >= 3.8 

Requires ‘cigar’ python package [‘pip install cigar’] 

$ python /Translate/translate/translate.py -t <transcripts_file_path> -q <queries_file_path> -o <output_folder_path> 
