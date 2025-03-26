# FILE : main.py
# WRITER : Elal Gilboa , elal.gilboa , 323083188
# EXERCISE : intro2cs final project 2025
# DESCRIPTION: This program execute the biosequence final project
# STUDENTS I DISCUSSED THE EXERCISE WITH: None
# WEB PAGES I USED: https://www.w3schools.com/python/ref_string_split.asp
#                   https://docs.python.org/3/library/argparse.html
#                   https://www.youtube.com/watch?v=6Q56r_fVqgw
#                   https://en.wikipedia.org/wiki/FASTQ_format
#                   https://docs.pytest.org/en/stable/how-to/monkeypatch.html
# NOTES: None

# Bio_sequence_project version 3.0
# changes:
#   1. validate_flag_combinations has been refactored into 4 simpler functions
#   2. ALN file format has been changed - to contain an Alignment object and not reads dict
#   3. __determine_status func in Alignment has been changed - to count the mapped
#      reads statuses in class self.__variables and also to set genome sources
#      for each read using read class variable "sources" which is used in dumpalign
#   4. load_fasta and load_fastq turned into generators functions
# !/usr/bin/env python3
###############################  Imports  #####################################
import argparse
import gzip, pickle
from typing import *
from collections import Counter
import os
import json
from functools import reduce


#############################################################################

class Kmer:
    """This class creates an K-mer object which is represented by a np array
    of strings. each cell in the array contain one of the following letters:
    A,T,C,G which represent the nucleotide bases in the DNA.
    parms:
    data = contain the np array of letters in length K
    length = the K length of this K-MER object
    specific = set to be True/False if the K-mer exist in many genomes on it is
    unique to only onr genome.
    id = a unique reference to identify this K-mer
    locations = a list of locations in which this Kmer appear in a Read
    sources_refs = list of identifiers of the K-mer source/s"""
    __id_counter: int = 0

    def __init__(self, sequence: str, locations: List[int]):
        self.__data: str = sequence
        self.__length: int = len(sequence)
        self.__specific: bool = False  # default data
        self.__id: int = Kmer.__id_counter
        Kmer.__id_counter += 1
        self.__loc_in_read: List[int] = locations
        self.__sources_refs: List[str] = []

    def __eq__(self, other: object) -> bool:
        """This function compare 2 K-mers object by their data values and length
        their id and location may be different"""
        if type(other) != Kmer:
            return False
        elif self.__length != cast(Kmer, other).__length:
            return False
        else:
            for i in range(self.__length):
                if self.__data[i] != cast(Kmer, other).__data[i]:
                    return False
            return True

    def __str__(self):
        """This func return a string to print a Kmer object"""
        return (f"K-mer id: {self.__id},"
                f" data: {self.__data}," +
                f" locations in read: {self.__loc_in_read}," +
                f" matching reference sources: {self.__sources_refs}," +
                f" specific K-mer: {self.__specific}")

    def __repr__(self):
        return (f"K-mer id: {self.__id},"
                f" data: {self.__data}," +
                f" locations in read: {self.__loc_in_read}," +
                f" matching reference sources: {self.__sources_refs}," +
                f" specific K-mer: {self.__specific}")

    def __len__(self) -> int:
        return self.__length

    def get_data(self) -> str:
        return self.__data

    def get_location(self):
        return self.__loc_in_read

    def get_id(self):
        return self.__id

    def get_specific_status(self):
        return self.__specific

    def set_specific(self, status: bool):
        """Func set specific status of Kmer - if it appears in more than one
        it is not specific and if it only appears in one than it's specific"""
        self.__specific = status

    def get_sources(self):
        return self.__sources_refs

    def add_source(self, sources: List[str]):
        """This func adds another sources to the reference sources list"""
        for source in sources:
            self.__sources_refs.append(source)

    def reverse_data(self) -> str:
        """this func return the reverse complement form of its sequence data"""
        new_data: List[str] = []
        for char in self.__data[::-1]:
            if char == "A":
                new_data.append("T")
            elif char == "T":
                new_data.append("A")
            elif char == 'C':
                new_data.append('G')
            elif char == 'G':
                new_data.append('C')
        return ''.join(new_data)


class Read:
    """This class represent a DNA read from an NGS machine, a read contain
    samples of valid and un valid sequences of A,T,C,G nucleotide bases in the
    DNA of an unknown object.
    parms:
    data = is a string of all A,T,C,G data received from the NGS machine
    quality = a string describing the quality of each nucleotide bases in data.
    the quality is represented by a chat in ASCII and the quality score is the
    ASCII VALUE-30 for each character
    status = could be unmapped / ambiguously_mapped / uniquely_mapped
    id = a unique identifier for the data file for each sample
    """

    def __init__(self, header: str, sequence: str, quality: str):
        self.__data: str = sequence
        self.__quality: str = quality
        self.__status: str = ""
        self.__header: str = header
        self.__sources: List[str] = []
        self.__orientation: str = "Forward"
        self.__reversed_data: str = self.__reverse_data(self.__data)

    def __str__(self):
        """This func return a string to print a Read object"""
        if self.__status == "":
            status = "NOT DEFINED"
            return (f"Read id: {self.__header} \n"
                    f"data: {self.__data} \n" +
                    f"quality: {self.__quality} \n" +
                    f"read mapping status: {status}")
        else:
            return (f"Read id: {self.__header} \n"
                    f"data: {self.__data} \n" +
                    f"quality: {self.__quality} \n" +
                    f"read mapping status: {self.__status}")

    def __eq__(self, other: object):
        if type(other) != Read:
            return False
        elif len(self.__data) != len(cast(Read, other).__data):
            return False
        elif self.__data != cast(Read, other).__data:
            return False
        else:
            return True

    def get_data(self) -> str:
        return self.__data

    def get_status(self) -> str:
        return self.__status

    def get_quality(self) -> str:
        return self.__quality

    def add_sources(self, sources: List[str]):
        for source in sources:
            self.__sources.append(source)

    def get_read_sources(self):
        return self.__sources

    def get_header(self):
        return self.__header

    def get_orient(self):
        return self.__orientation

    def set_orient(self, orient: str):
        self.__orientation = orient

    def __reverse_data(self, sequence: str) -> str:
        """this func return the reverse complement form of sequence data"""
        new_data: List[str] = []
        for char in sequence[::-1]:
            if char == "A":
                new_data.append("T")
            elif char == "T":
                new_data.append("A")
            elif char == 'C':
                new_data.append('G')
            elif char == 'G':
                new_data.append('C')
        return ''.join(new_data)

    def set_status(self, status: str) -> bool:
        """This func set Read status to be as received from the mapping process
        and could only be one of 3 valid statuses"""
        VALID_STATUSES = ["Unmapped", "Ambiguous", "Unique"]
        if status not in VALID_STATUSES:
            return False
        else:
            self.__status = status
            return True

    def extract_kmers(self, k):
        """This func return a dictionary of all k-mers in length K from this read
        as Keys and their location/s in data as Value"""
        data = self.__data
        k_mers_dict: Dict[str, List[int]] = {}
        if k > len(data):
            raise ValueError(
                f"K size is smaller than read's length {self.__header}")
        for i in range(len(data) - k + 1):
            kmer: str = data[i:i + k]
            if 'N' in kmer:
                continue
            if kmer not in k_mers_dict:
                k_mers_dict[kmer] = [i]
            else:
                k_mers_dict[kmer].append(i)
        return k_mers_dict


class Reference:
    """This class represent a genome reference from a FASTA file. Each genome
    has a unique header and the content of the reference made of A,T,C,G.
    parms:
    genome = a string representing the reference data
    header = a strig representing the reference identifier"""
    __id_counter: int = 0

    def __init__(self, header: str, sequence: str):
        self.__genome: str = sequence
        self.__header: str = header
        self.__kmers: Dict[str, List[int]] = {}
        self.__id: int = Reference.__id_counter
        Reference.__id_counter += 1

    def __str__(self):
        """This func return a string to print a Reference object"""
        return (f"reference header:{self.__header}, \n"
                f"genome data:{self.__genome}")

    def __len__(self) -> int:
        return len(self.__genome)

    def get_header(self):
        return self.__header

    def get_genome(self):
        return self.__genome

    def get_kmers(self):
        return self.__kmers

    def get_id(self):
        return self.__id

    def extract_kmers(self, k):
        """This func return a dictionary of all k-mers in length K from this
        reference as Keys and their location/s in data as Value"""
        data = self.__genome
        k_mers_dict: Dict[str, List[int]] = {}
        if k > len(data):
            raise ValueError(
                f"K size is smaller than genome's length {self.__header}")
        for i in range(len(data) - k + 1):
            kmer = data[i:i + k]
            if 'N' in kmer:
                continue
            if kmer not in k_mers_dict:
                k_mers_dict[kmer] = [i]
            else:
                k_mers_dict[kmer].append(i)
        self.__kmers = k_mers_dict
        return k_mers_dict


class KmersDB:
    """K-mers DB is built as a dictionary in which contain as keys all Kmers
    string content (from all references together) and as values a dictionary
    which contain a ref_header as key and a list of location in which this
    K-mer appear in it as value.
    parms:
    refs_kmers = the main database, the reason the kmers string are the keys
    it because string are hashable so the search in this major database should
    be more efficient
    """

    def __init__(self):
        self.__refs_kmers: Dict[str, Dict[str, List[int]]] = {}
        #  {k_mers_str_from_all_refs: {ref_head: List[k_mer_loc_in_ref]}
        self.__similarity: Dict[str, Dict[str, Any]] = {}

    def __str__(self):
        """This func return a printable string of the DB"""
        result = ""
        DB = self.__refs_kmers
        for kmer_key, locs_in_refs_dict in DB.items():
            result += f"K-mer: {kmer_key} \n"
            for ref_header, lst_of_locs in locs_in_refs_dict.items():
                result += (
                        f"appear in reference header: {ref_header}" +
                        f" in locations: {lst_of_locs} \n")
        return result

    def inset_ref_data(self, ref_header: str,
                       kmers_dict: Dict[str, List[int]]) -> bool:
        """This func receive a reference_header and all k_mers extracted from
        this reference and insert them to DB"""
        for kmer_data, locs_in_ref in kmers_dict.items():
            if kmer_data not in self.__refs_kmers:
                self.__refs_kmers[kmer_data] = {}
                self.__refs_kmers[kmer_data][ref_header] = locs_in_ref
            else:
                if ref_header in self.__refs_kmers[kmer_data]:
                    raise ValueError(
                        f"Duplicate ref_header '{ref_header}' for K-mer '{kmer_data}'")
                self.__refs_kmers[kmer_data][ref_header] = locs_in_ref
        return True

    def get_ref_db(self) -> Dict[str, Dict[str, List[int]]]:
        """This func return the build DB of reference Kmers"""
        return self.__refs_kmers

    def get_kmer_size(self) -> int:
        """This func return the size of each Kmer key in DB"""
        for key, value in self.__refs_kmers.items():
            return len(key)
        return 0

    def remove_ref_from_inner_dict(self, ref_header: str) -> None:
        """This func receive a ref_header and remove it from all the inner dicts
        in DB. Func return True if action is successful and False otherwise"""
        for kmer in list(self.__refs_kmers.keys()):
            if ref_header in self.__refs_kmers[kmer]:
                del self.__refs_kmers[kmer][ref_header]
                # if no ref_headers remain for this kmer, delete the k-mer key
                # from outer Dict
                if not self.__refs_kmers[kmer]:
                    del self.__refs_kmers[kmer]

    def get_value(self, key: str):
        """This func receive a Key and return the value from DB for this Key.
        If key does not exist raise a KeyError"""
        if key in self.__refs_kmers:
            return self.__refs_kmers[key]
        else:
            raise KeyError(f"key {key} does not exist in Ref Kmers DB")

    def __repr__(self):
        """This func return a printable string of the DB"""
        result = ""
        DB = self.__refs_kmers
        for kmer_key, locs_in_refs_dict in DB.items():
            result += f"K-mer: {kmer_key} \n"
            for ref_header, lst_of_locs in locs_in_refs_dict.items():
                result += (
                        f"appear in reference header: {ref_header}" +
                        f" in locations: {lst_of_locs} \n")
        return result

    def set_similarity(self, similarity) -> None:
        """This func recieve a similarity dictionary and save it in DB. it is
        used when used choose reference command and do not print results"""
        self.__similarity = similarity

    def get_similarity(self):
        """This func return a similarity dictionary - to be printed in dumpref"""
        return self.__similarity


class Alignment:
    """This class represent an alignment file between a DB of references genomes
    and Reads of DNA samples from an unknown organism and applying the Pseudo
    Alignment algorithm.
    parms:
    references = a dictionary of Reference objects {ref_header, Reference Object}
    Kmers_DB = a Kmer_DB object of all kmers extracted from References
    {kmers_data_from_all_refs: {ref_header: List[kmers_locs_in_ref}}
    reads = a dictionary of Read object received in a FASTAQ file
    {read_header: Read object}
    reads_kmers_db = a dictionary of {read_header: Kmers_for_each_read} a list
    of all Kmers extracted from a read.
    total_unmapped = a counter for all Kmers in reads_kmers_db which received a
    status unmapped.
    total_ambiguous = same as unmapped
    total_unique = same as unmapped"""

    def __init__(self):
        self.__references: List[str] = []
        self.__refs_kmers_db: KmersDB = KmersDB()
        self.__reads: Dict[str, Read] = {}
        self.__reads_kmers_db: Dict[str, List[Kmer]] = {}
        self.__total_unmapped: int = 0
        self.__total_ambiguous: int = 0
        self.__total_unique: int = 0
        self.__coverage_sum: Dict[str, Dict[str, Any]] = {}
        self.__quality_stats: List[int] = []

    def set_references(self, references: List[str]) -> None:
        """Func receive a list of genomes' ref headers and set it to as
         self.__references."""
        self.__references = references

    def set_reads(self, reads: Dict[str, Read]) -> None:
        """func receive a dict or reads from FASTAQ file and set it to reads.
        Dict in format {read_header, Read object}"""
        self.__reads = reads

    def set_ref_kmers_db(self, ref_kmers_db: KmersDB):
        self.__refs_kmers_db = ref_kmers_db

    def get_reads(self):
        return self.__reads

    def get_references(self):
        return self.__references

    def get_ref_kmers_DB(self):
        return self.__refs_kmers_db

    def get_reads_kmers_db(self):
        return self.__reads_kmers_db

    def set_quality_stats(self, quality_stats: List[int]):
        self.__quality_stats = quality_stats

    def get_quality_stats(self):
        return self.__quality_stats

    def set_coverage_sum(self, coverage: Dict[str, Dict[str, Any]]):
        self.__coverage_sum = coverage

    def get_coverage_sum(self):
        return self.__coverage_sum

    def get_reads_stats(self) -> Dict[str, int]:
        """This func return the reads statistics data in a dictionary with each
        data and it's value"""
        UNIQUE, AMBIG, UNMAPPED = "unique_mapped_reads", "ambiguous_mapped_reads", "unmapped_reads"
        reads_stats: Dict[str, int] = {UNIQUE: self.__total_unique,
                                       AMBIG: self.__total_ambiguous,
                                       UNMAPPED: self.__total_unmapped}
        return reads_stats

    def __calc_mean_quality(self, sequence: str):
        """This func receive a read and calculates it's mean quality. Func return
        a float which is the mean quality score of this read"""
        sum_ascii_score = lambda total, char: total + (ord(char) - 33)
        result = reduce(sum_ascii_score, sequence, 0) / len(sequence)
        return result

    def build_reads_kmers_db(self, k_size: int, mkq: int) -> int:
        """This func receive a K int and crete a DB of all Kmers in length K
        from all reads in self.__reads. Func create a Dict[str, List[Kmer]] a
        dictionary of kmer objects matching to each read header as key. Func
        filter out of DB low quality reads in case mkq is not None and return
        the number of filtered out Kmers"""
        reads_kmers_db: Dict[str, List[Kmer]] = {}
        filtered_out: int = 0
        for read_header, read in self.__reads.items():
            read_kmers_dict: Dict[str, List[int]] = read.extract_kmers(k_size)
            list_of_kmers: List[Kmer] = []
            quality: str = read.get_quality()
            for kmer_data, locs_in_read in read_kmers_dict.items():
                if mkq is not None:
                    first_appearance_index: int = locs_in_read[0]
                    kmer_quality: str = quality[
                                        first_appearance_index: first_appearance_index + k_size]
                    mean_quality: float = self.__calc_mean_quality(
                        kmer_quality)
                    if mean_quality >= mkq:
                        new_kmer = Kmer(kmer_data, locs_in_read)
                        list_of_kmers.append(new_kmer)
                    else:
                        filtered_out += 1
                else:
                    new_kmer = Kmer(kmer_data, locs_in_read)
                    list_of_kmers.append(new_kmer)
            if read_header not in reads_kmers_db:
                reads_kmers_db[read_header] = list_of_kmers
            else:
                raise KeyError(
                    f"read Header {read_header} already exist in reads_kmers_db")
        self.__reads_kmers_db = reads_kmers_db
        return filtered_out

    def __find_matching_ref(self, sequence: str) -> Dict[str, List[int]]:
        """Func receive a Kmer data sequence and check the ref_kmers_db for
        matches. func return the inner dictionary if there's a match (one or many)
        and return {} in there's No match. The search in ref_kmers_db is O(1)
        because sequence is the key of this DB and it's hashable"""
        try:
            return self.__refs_kmers_db.get_value(sequence)
        except KeyError:
            return {}

    def __count_sources(self, kmers_list: List[Kmer]) -> Counter:
        """This func receive a list of Kmers and return a Counter object which
        contains Kmer sources as Keys and their count from kmer_list"""
        kmers_sources_count: Counter = Counter()
        for kmer in kmers_list:
            for source in kmer.get_sources():
                if source in kmers_sources_count:
                    kmers_sources_count[source] += 1
                else:
                    kmers_sources_count[source] = 1
        return kmers_sources_count

    def __validate_unique_status(self, map_count: int,
                                 combined_counts: Counter,
                                 p: int) -> str:
        """This func is called after mapping process determined a read as
        Unique and validate it using the unspecific Kmers. Func receive the
         assumed source reference header, the number of Kmers pointing to it,
         a combined Counter object of all Kmers results"""
        AMBIGUOUS, UNIQUE = "Ambiguous", "Unique"
        combined_max_source: str = max(combined_counts,
                                       key=lambda k: combined_counts[k])
        combined_max_count: int = combined_counts[combined_max_source]
        if combined_max_count - map_count > p:
            return AMBIGUOUS
        else:
            return UNIQUE

    def __extract_ambiguous_sources(self, combined_count: Counter,
                                    map_count: int, original_source: str) -> \
            List[str]:
        """
        This func receive a Counter of specific and unspecific Kmers for cases
        in which mapping process resulted Unique status and then validated and
        turned into ambiguous. In these cases the read's sources are the original
        mapped genome and all the other genomes received higher count of combined
        Kmers. Func return a list of read's sources"""
        ambiguous_sources = [
            genome for genome, count in combined_count.items()
            if count > map_count]
        if original_source not in ambiguous_sources:
            ambiguous_sources.append(original_source)
        return ambiguous_sources

    def __determine_status(self, specific_kmers: List[Kmer],
                           unspecific_kmers: List[Kmer], m: int,
                           p: int) -> Tuple[str, List[str]]:
        """Func receive a list of Kmers specific and unspecific received from
        a read and determine the read's mapping status. m and p are thresholds
        for minimum difference between highest and second highest mapped genomes.
        Func return the determined status and a list of read's sources for later
         use in dumpalign"""
        UNMAPPED, AMBIGUOUS, UNIQUE = "Unmapped", "Ambiguous", "Unique"
        REF_HEADER_LOC, MATCHES_COUNT = 0, 1
        FIRST, SECOND = 0, 1
        status: str = ""
        combined_count: Counter = Counter()
        ambiguous_sources: List[str] = []
        read_sources: List[str] = []
        if len(specific_kmers) == 0 and len(unspecific_kmers) == 0:
            # none of the K-mers in read match any of the genomes in DB.
            self.__total_unmapped += 1
            return UNMAPPED, read_sources
        specific_source_count: Counter = self.__count_sources(specific_kmers)
        unspecific_source_count: Counter = self.__count_sources(
            unspecific_kmers)
        # {source_ref_header, num_kmers_pointing_ref}
        if len(specific_source_count) == 0:
            #  all specific Kmers are unmapped to any genome - so read is unmapped
            self.__total_unmapped += 1
            return UNMAPPED, read_sources
        elif len(specific_source_count) == 1:
            # all specific Kmers point to the same genome
            unique_genome: str = list(specific_source_count.keys())[0]
            if p >= 1:
                # mapping process determined unique, yet needs to be validated
                map_count1: int = max(specific_source_count.values())
                status = self.__validate_unique_status(map_count1,
                                                       specific_source_count + unspecific_source_count,
                                                       p)
                if status == UNIQUE:
                    read_sources.append(unique_genome)
                    self.__total_unique += 1
                    return UNIQUE, read_sources
                else:
                    # Extract ambiguous sources
                    combined_count = specific_source_count + unspecific_source_count
                    ambiguous_sources = self.__extract_ambiguous_sources(
                        combined_count, map_count1, unique_genome)
                    for source in ambiguous_sources:
                        read_sources.append(source)
                    self.__total_ambiguous += 1
                    return AMBIGUOUS, read_sources
            else:
                read_sources.append(unique_genome)
                self.__total_unique += 1
                return UNIQUE, read_sources
        else:
            # specific genomes point to more than one genome
            top_two: List[
                Tuple[str, int]] = specific_source_count.most_common(2)
            most_common: Tuple[str, int] = top_two[FIRST]
            sec_most_common: Tuple[str, int] = top_two[SECOND]
            if cast(int, most_common[MATCHES_COUNT]) - cast(int,
                                                            sec_most_common[
                                                                MATCHES_COUNT]) >= m:
                if p >= 1:
                    # mapping process determined unique, yet needs to be validated
                    map_count: int = cast(int, most_common[MATCHES_COUNT])
                    max_map_source: str = cast(str,
                                               most_common[REF_HEADER_LOC])
                    status = self.__validate_unique_status(
                        map_count,
                        specific_source_count + unspecific_source_count, p)
                    if status == UNIQUE:
                        read_sources.append(max_map_source)
                        self.__total_unique += 1
                        return UNIQUE, read_sources
                    else:
                        # Extract ambiguous sources
                        combined_count = specific_source_count + unspecific_source_count
                        ambiguous_sources = self.__extract_ambiguous_sources(
                            combined_count, map_count, max_map_source)
                        for source in ambiguous_sources:
                            read_sources.append(source)
                        self.__total_ambiguous += 1
                        return AMBIGUOUS, read_sources
                else:
                    read_source: str = cast(str, most_common[REF_HEADER_LOC])
                    read_sources.append(read_source)
                    self.__total_unique += 1
                    return UNIQUE, read_sources
            else:
                # the difference between most and second common is less than
                # m  which is a threshold for minimum difference, by default it's 1
                ambiguous_sources = list(specific_source_count.keys())
                # in this case all genomes which this read was partially mapped
                # to are countable and added to read's source list
                for source in ambiguous_sources:
                    read_sources.append(source)
                self.__total_ambiguous += 1
                return AMBIGUOUS, read_sources

    def __reverse_lst_of_kmers(self, lst_of_kmers: List[Kmer]) -> List[Kmer]:
        """This func receive a list of Kmer objects and return a reversed list
        of Kmers with reversed data sequence. Func created new reversed Kmer from
        each original Kmer. Notice - the locations in reversed Kmer are not accurate
        they match the original and not the new Kmer."""
        reversed_lst: List[Kmer] = []
        for kmer in lst_of_kmers:
            data: str = kmer.reverse_data()
            locations: List[int] = kmer.get_location()
            rev_kmer = Kmer(data, locations)
            reversed_lst.append(rev_kmer)
        return reversed_lst

    def __sort_specific_kmers(self, lst_of_kmers: List[Kmer], mg: int) -> \
            Tuple[List[Kmer], List[Kmer], int]:
        """This func receive a list of Kmers and sort them into specific and
        unspecific based on their match to a genome using __find_matching_ref().
        func also receive a threshold mg, if it's not none func remove kmers
        which are highly redundant from mapping process. Func return the specific
        and unspecific Kmers list and the num of removed Kmers"""
        SPECIFIC_KMER = True
        UNSPECIFIC_KMER = False
        filtered_hq_kmers: int = 0
        specific_kmers: List[Kmer] = []
        unspecific_kmers: List[Kmer] = []
        sources: List[str] = []
        for kmer in lst_of_kmers:
            matches: Dict[str, List[int]] = self.__find_matching_ref(
                kmer.get_data())
            if len(matches) > 1:  # more than one genome match this k-mer
                sources = list(matches.keys())
                if mg is not None:
                    if len(sources) <= mg:
                        kmer.add_source(sources)
                        kmer.set_specific(UNSPECIFIC_KMER)
                        unspecific_kmers.append(kmer)
                    else:
                        filtered_hq_kmers += 1
                else:
                    kmer.add_source(sources)
                    kmer.set_specific(UNSPECIFIC_KMER)
                    unspecific_kmers.append(kmer)
            elif len(matches) == 1:
                # there is only one genome matching
                sources = list(matches.keys())
                kmer.add_source(sources)
                kmer.set_specific(SPECIFIC_KMER)
                specific_kmers.append(kmer)
        return specific_kmers, unspecific_kmers, filtered_hq_kmers

    def __determine_orientation(self, status_f: str, status_r: str,
                                specific_kmers: List[Kmer],
                                unspecific_kmers: List[Kmer],
                                rev_specific_kmers: List[Kmer],
                                rev_unspecific_kmers: List[Kmer]) -> bool:
        """This func receive a determined orientation for forward and reversed
        read sequence, Func return True if read should be handled reversed and
        False if it should remain Forward"""
        UNIQUE, AMBIG = "Unique", "Ambiguous"
        if status_f == UNIQUE and status_r != UNIQUE:
            # continue with forward
            return False
        elif status_r == UNIQUE and status_f != UNIQUE:
            # continue with reverse
            return True
        elif status_f == UNIQUE and status_r == UNIQUE:
            # chose the one with more unique Kmers, choose forward if tied
            return len(rev_specific_kmers) > len(specific_kmers)
        elif status_f == AMBIG and status_r == AMBIG:
            # chose the orientation with more total matched Kmers
            # sum the unique and ambig kmers for each orientation (if kmer is
            # unmapped it will not be on either lists)
            total_kmers_f = len(specific_kmers) + len(unspecific_kmers)
            total_kmers_r = len(rev_specific_kmers) + len(rev_unspecific_kmers)
            return total_kmers_r > total_kmers_f
        elif status_r == AMBIG:
            # forward is unmapped so reverse is chosen
            return True
        else:
            # choose forward - both orientations unmapped - choose forward
            return False

    def map_reads(self, m: int, p: int, mg: int, rev_comp: bool) -> int:
        """Func execute the pseudo alignment algorithm, and map all reads in
        self.__reads. func set the status of each read to be the result of the
        mapping process. For each read func check for matches in ref_kmers_db,
        and if there is a match updates the Kmer sources and the read status
        parms: p - unique-threshold, m - ambiguous-threshold.
        mg - maximus genomes (for quality ext), rev_comp - flag (for reverse ext)
        Func return the number of highly redundant kmers, filtered out and updates
        the read sources for each read.
        """
        REV: str = "Reverse"
        filtered_hr_kmers: int = 0
        rev_filtered_hr_kmers: int = 0
        sum_filtered_hr_kmers: int = 0
        status: str = ""
        for read_header, lst_of_kmers in self.__reads_kmers_db.items():
            if not rev_comp:
                specific_kmers, unspecific_kmers, filtered_hr_kmers = self.__sort_specific_kmers(
                    lst_of_kmers, mg)
                sum_filtered_hr_kmers += filtered_hr_kmers
                status, read_sources = self.__determine_status(
                    specific_kmers, unspecific_kmers, m, p)
                read = self.__reads[read_header]
                read.set_status(status)
                read.add_sources(read_sources)
            else:
                # sort Kmers in forward orientation:
                specific_kmers, unspecific_kmers, filtered_hr_kmers = self.__sort_specific_kmers(
                    lst_of_kmers, mg)
                # sort Kmers in reverse orientation:
                rev_lst_of_kmers = self.__reverse_lst_of_kmers(lst_of_kmers)
                rev_specific_kmers, rev_unspecific_kmers, rev_filtered_hr_kmers = self.__sort_specific_kmers(
                    rev_lst_of_kmers, mg)
                # determine status for both orientations:
                status_f, read_sources_f = self.__determine_status(
                    specific_kmers, unspecific_kmers, m, p)
                status_r, read_sources_r = self.__determine_status(
                    rev_specific_kmers, rev_unspecific_kmers, m, p)
                if self.__determine_orientation(status_f, status_r,
                                                specific_kmers,
                                                unspecific_kmers,
                                                rev_specific_kmers,
                                                rev_unspecific_kmers):
                    # continue with reverse
                    sum_filtered_hr_kmers += rev_filtered_hr_kmers
                    read = self.__reads[read_header]
                    read.set_status(status_r)
                    read.set_orient(REV)
                    read.add_sources(read_sources_r)
                else:
                    # continue with forward
                    sum_filtered_hr_kmers += filtered_hr_kmers
                    read = self.__reads[read_header]
                    read.set_status(status_f)
                    read.add_sources(read_sources_f)
        return sum_filtered_hr_kmers


###################################################################

def gen_read_from_fastq(filename) -> Generator:
    """This func is a generator of necessary data: header, sequence and quality
    in order to create a read object. Func yield this data from FASTQ file and
    return a Tuple of """
    with open(filename, 'r') as file:
        while True:
            read_header: str = file.readline().strip()
            if not read_header:
                # generator reached EOF
                break
            if not read_header.startswith('@'):
                file.close()
                raise ValueError(
                    f"Received file {filename} format is invalid. This header: {read_header} does not start with '@'")
            sequence: str = file.readline().strip()
            sep: str = file.readline().strip()
            if not sep.startswith('+'):
                file.close()
                raise ValueError(
                    f"Received file {filename} format is invalid. This seperator line: {sep} does not start with '+'")
            quality = file.readline().strip()
            if len(sequence) != len(quality):
                file.close()
                raise ValueError(
                    f"The lengths of sequence and quality don't match in this read: {read_header[1:]}\n"
                    f"Sequence: {sequence}\nQuality: {quality}")
            yield read_header[1:], sequence, quality


def generate_reads_dict(filename: str) -> Dict[str, Read]:
    """Func receive a filename and load the FASTQ file reads Data. Func creates
    Read objects from each read in file and return a dictionary of Read objects
    {read_header, Read_obj}. Func is using a generator for receiving the next:
    reade_header, sequence and quality from FASTQ file"""
    reads: Dict[str, Read] = {}
    read_generator: Generator = gen_read_from_fastq(filename)
    for read_header, sequence, quality in read_generator:
        new_read: Read = Read(read_header, sequence, quality)
        if read_header not in reads:
            reads[read_header] = new_read
        else:
            raise KeyError(
                f"There are two Reads with the same header in FASTQ file {filename}")
    return reads


def generate_ref_dict(filename: str) -> Dict[str, Reference]:
    """Func receive a filename and load the FASTA file of Genomes references.
    Func creates Reference objects from each genome in file and return a dictionary
    of Reference objects {ref_header, Reference_obj}"""
    references: Dict[str, Reference] = {}
    refs_generator: Generator = gen_ref_from_fasta(filename)
    for ref_header, sequence in refs_generator:
        new_ref: Reference = Reference(ref_header, sequence)
        if ref_header not in references:
            references[ref_header] = new_ref
        else:
            raise KeyError(
                f"There are two Genomes with the same header in FASTA file {filename}")
    return references


def gen_ref_from_fasta(filename) -> Generator:
    """func creates a generator for Reference object, for each yield func return
    a ref_header and a sequence. If FASTA format is invalid func raise Exception"""
    with open(filename, 'r') as file:
        ref_header = None
        sequence: List[str] = []
        for line in file:
            # file object is iterable and yield the next line
            if not line:
                continue
            new_line: str = line.strip()
            if new_line.startswith('>'):
                # mark the start of a new ref_header
                if ref_header:
                    if not sequence:
                        raise ValueError(
                            f"Invalid FASTA format: Missing description or identifier after '>'"
                        )
                    yield ref_header, ''.join(sequence)
                ref_header = new_line[1:]
                if not ref_header:
                    raise ValueError(
                        f"Invalid FASTA format: Header line is empty")
                sequence = []
            else:
                # this line is part of the ref_data_sequence
                sequence.append(new_line)
        if ref_header:
            if not sequence:
                raise ValueError(
                    f"Invalid FASTA format: Missing description or identifier after '>'"
                )
            yield ref_header, "".join(sequence)
        else:
            raise ValueError("Invalid FASTA format: Header line is empty")


def dump_alignment_results(aln: Alignment, filename: str) -> None:
    """Func receive an Alignment object and dump it into an ALN file using
    pickle and gzip to compress to bytestream"""
    with gzip.open(filename, "wb") as file:
        pickle.dump(aln, file)


def load_aln_file(aln_file_path: str) -> Alignment:
    """This func receive an aln_file_path and return Alignment object"""
    with gzip.open(aln_file_path, "rb") as file:
        aln: Alignment = pickle.load(file)
    return aln


def dump_ref_kmers_db(ref_kmers_db: KmersDB, genome_bases: Dict[str, int],
                      filename: str) -> None:
    """Func receive a Kmers_db object and serializes this object into bytes.
    Then func save this serialized bytestream as a compressed pickle file
    """
    with gzip.open(filename, "wb") as file:
        pickle.dump(ref_kmers_db, file)
        pickle.dump(genome_bases, file)


def load_kdb_file(filename) -> Tuple[KmersDB, Dict[str, int]]:
    """Func receive a path to a KDB file and return a KmersDB object and genome
    bases dictionary. using unpickling method
    """
    with gzip.open(filename, "rb") as file:
        ref_kmers_db: KmersDB = pickle.load(file)
        genome_bases: Dict[str, int] = pickle.load(file)
    return ref_kmers_db, genome_bases


def dump_aln_file(reads_mapping_results: Dict[str, Read],
                  filename: str) -> None:
    """Func receive a dictionary of reads results from the mapping process and
    save them into an aln file using pickle and gzip modules.
    """
    with gzip.open(filename, "wb") as file:
        pickle.dump(reads_mapping_results, file)


def load_mapping_results(filename: str) -> Dict[str, Read]:
    """Func receive a dictionary of reads results from the mapping process and
    save them into an aln file using pickle and gzip modules.
    """
    with gzip.open(filename, "rb") as file:
        reads_mapping_results: Dict[str, Read] = pickle.load(file)
    return reads_mapping_results


def readargs(args=None) -> argparse.Namespace:
    """This function return a ParserArgs object which contain all the received
    required and optional arguments from command line"""
    parser = argparse.ArgumentParser(prog='Biosequence project', )
    # General arguments
    parser.add_argument('-t', '--task', help="task", required=True)
    parser.add_argument('-g', '--genomefile',
                        help="Genome fasta file (multiple records)", )
    parser.add_argument('-r', '--referencefile',
                        help="kdb file. Can be either input or name for output file", )
    parser.add_argument('-k', '--kmer-size', type=int,
                        help="length of kmers", )
    # Task specific arguments
    # align
    parser.add_argument('-a', '--alignfile',
                        help="aln file. Can be either input or name for output file", )
    parser.add_argument('--reads', help="fastq reads file", )
    parser.add_argument('-m', '--unique-threshold',
                        help="unique k-mer threshold", default=1, type=int)
    parser.add_argument('-p', '--ambiguous-threhold',
                        help="ambiguous k-mer threshold", default=1, type=int)
    # align+rc
    parser.add_argument('--reverse-complement', action='store_true', )
    # align+quality
    parser.add_argument('--min-read-quality', type=int, )
    parser.add_argument('--min-kmer-quality', type=int, )
    parser.add_argument('--max-genomes', type=int, )
    # coverage
    parser.add_argument('--genomes', type=str)
    parser.add_argument('--coverage', action='store_true')
    parser.add_argument('--window-size', type=int, default=100, )
    parser.add_argument('--min-coverage', type=int, default=1, )
    parser.add_argument('--full-coverage', action='store_true', )
    # similarity
    parser.add_argument('--filter-similar', action='store_true', )
    parser.add_argument('--similarity-threshold', type=float, default=0.95, )
    return parser.parse_args(args)


def calc_mean_quality(sequence: str) -> float:
    """This func receive a read and calculates it's mean quality. Func return
    a float which is the mean quality score of this read"""
    sum_ascii_score = lambda total, char: total + (ord(char) - 33)
    result = reduce(sum_ascii_score, sequence, 0) / len(sequence)
    return result


def filter_lq_reads(reads: Dict[str, Read], mrq: int) -> Tuple[
    Dict[str, Read], int]:
    """This func receive a Reads dictionary and minimum_read_quality threshold
    and return filtered Reads according to threshold"""
    filtered_reads: Dict[str, Read] = {}
    filtered_out: int = 0
    for read_header, read in reads.items():
        mean_quality: float = calc_mean_quality(read.get_quality())
        if mean_quality >= mrq:
            filtered_reads[read_header] = read
        else:
            filtered_out += 1
    return filtered_reads, filtered_out


def build_ref_kmers_db(refs_dict: Dict[str, Reference], k: int) -> Tuple[
    KmersDB, Dict[str, int]]:
    """
    Func Iterates over references and use extract_kmers(k) to create a dict
    with all Kmers from each Reference and saves it to Reference Kmers DB.
    If Reference was incorrectly saved to DB an Exception will be raised
    """
    ref_kmers_DB: KmersDB = KmersDB()
    genome_bases: Dict[str, int] = {}
    for ref_header, ref in refs_dict.items():
        ref_kmers_dict: Dict[str, List[int]] = ref.extract_kmers(k)
        genome_bases[ref_header] = len(ref)
        result = ref_kmers_DB.inset_ref_data(ref_header, ref_kmers_dict)
        if not result:
            raise ValueError(
                f"{ref_header} data was not saved correctly to Reference Kmers DB")
    return ref_kmers_DB, genome_bases


def ascending_sort_refs(refs_dict: Dict[str, Reference],
                        refs_kmers_db: KmersDB) -> List[Reference]:
    """This func receive a dictionary of Reference objects and return a list
    of sorted References according to the following rules:
    1. num of specific Kmers - only count once each Kmer
    2. total number of Kmers
    3. genome length
    4. insertion to db order (id)
    """
    sorted_refs: List[Reference] = sorted(refs_dict.values(), key=lambda ref: (
        count_unique_kmers(ref, refs_kmers_db), len(ref.get_kmers()),
        len(ref.get_genome()), ref.get_id()))
    return sorted_refs


def calc_sim_score(ref_kmers: Set[str], compare_ref_kmers: Set[str]) -> float:
    """This func receive 2 sets of Kmers strings, intersect them and calc the
    similarity score of these sets. If both sets are equal so the score will be
    1. in all other cases score will be less than 1. If they are strangers 0"""
    min_size: int = min(len(ref_kmers), len(compare_ref_kmers))
    return len(ref_kmers & compare_ref_kmers) / min_size


def count_unique_kmers(ref: Reference, ref_kmers_db: KmersDB) -> int:
    """this func count how many unique kmers are in this specific reference
    and uses the ref_kmers_db to count. Func return an int which is this number.
    This function count Kmer as unique multiple times if it appears multiple times
    in a genome"""
    kmers: List[str] = list(ref.get_kmers().keys())
    unique_counter: int = 0
    for key in kmers:
        appearances: Dict[str, List[int]] = ref_kmers_db.get_value(key)
        if len(appearances) == 1:
            # this Kmer is unique to this genome
            appr_in_gen: int = len(list(appearances.values())[0])
            # if kmer is unique and appear multiple times in genome it should be
            # counted few times
            unique_counter += appr_in_gen
    return unique_counter


def find_original_stats(refs_dict: Dict[str, Reference],
                        ref_kmers_db: KmersDB) -> Dict[str, Dict[str, int]]:
    """This function receive a reference dictionary and return a Dict[str, Dict[str, int]]
    mapping each ref_header to it's total Kmer number and it unique Kmer's number"""
    UNIQUE, TOTAL = "unique_kmers", "total_kmers"
    original_stats: Dict[str, Dict[str, int]] = {
        ref.get_header(): {UNIQUE: 0, TOTAL: 0}
        for ref in refs_dict.values()}
    for ref_header, ref in refs_dict.items():
        ref_kmers: Dict[str, List[int]] = ref.get_kmers()
        for kmer, appearances in ref_kmers.items():
            original_stats[ref_header][TOTAL] += len(appearances)
            if len(ref_kmers_db.get_value(kmer)) == 1:
                original_stats[ref_header][UNIQUE] += len(appearances)
    return original_stats


def filter_ref_kmers_db(refs_dict: Dict[str, Reference],
                        refs_kmers_db: KmersDB, genome_bases: Dict[str, int],
                        sim_thresh: float) -> Dict[str, Dict[str, Any]]:
    """This func receive a dictionary of Reference objects and filter it according
    to similarity threshold. Func sort the references in ascending order and
    then filter the highly similar genomes. Func return change the original DB
    and return the Similarity dict with all data"""
    KEPT, UNIQUE, TOTAL = "kept", "unique_kmers", "total_kmers"
    LENGTH, SIMILAR, SCORE = "genome_length", "similar_to", "similarity_score"
    similarity: Dict[str, Dict[str, Any]] = {}
    original_stats: Dict[str, Dict[str, int]] = find_original_stats(refs_dict,
                                                                    refs_kmers_db)
    sorted_refs: List[Reference] = ascending_sort_refs(refs_dict,
                                                       refs_kmers_db)
    unique_kmers_count: Dict[str, int] = {
        ref.get_header(): count_unique_kmers(ref, refs_kmers_db) for ref in
        sorted_refs}
    for index, ref in enumerate(sorted_refs):
        ref_kmers: Set[str] = set(ref.get_kmers().keys())
        filtered_flag: bool = False
        for compare_ref in sorted_refs[index + 1:]:
            compare_ref_kmers: Set[str] = set(compare_ref.get_kmers().keys())
            sim_score: float = calc_sim_score(ref_kmers, compare_ref_kmers)
            if sim_score > sim_thresh:
                similarity[ref.get_header()] = {
                    KEPT: 'no',
                    UNIQUE: original_stats[ref.get_header()][UNIQUE],
                    TOTAL: original_stats[ref.get_header()][TOTAL],
                    LENGTH: len(ref.get_genome()),
                    SIMILAR: compare_ref.get_header(),
                    SCORE: sim_score}
                refs_kmers_db.remove_ref_from_inner_dict(ref.get_header())
                # if ref_header is deleted then it needs to be deleted from
                # genome_bases as well. Change conducted in-place
                del genome_bases[ref.get_header()]
                filtered_flag = True
                break
        if not filtered_flag:
            similarity[ref.get_header()] = {
                KEPT: 'yes',
                UNIQUE: original_stats[ref.get_header()][UNIQUE],
                TOTAL: original_stats[ref.get_header()][TOTAL],
                LENGTH: len(ref.get_genome()),
                SIMILAR: "NA",
                SCORE: "NA"}
    return similarity


def execute_reference(genome_file_path: str, k_size: int,
                      output_file: str, filter_sim: bool,
                      sim_thresh: float) -> None:
    """This function receive necessary arguments to run reference command:
     genome ref FASTA file path, K-mer size, output_file path to save REF_KMER_DB
    in KDB file format. Func return None and create the KDB file in current folder"""
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_file_path)
    tup = build_ref_kmers_db(refs_dict, k_size)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    if filter_sim:
        similarity = filter_ref_kmers_db(refs_dict, refs_kmers_db,
                                         genomes_bases, sim_thresh)
        similarity_summary: Dict[str, Dict[str, Any]] = {
            "Similarity": similarity}
        refs_kmers_db.set_similarity(similarity_summary)
    dump_ref_kmers_db(refs_kmers_db, genomes_bases, output_file)


def combine_gen_sum_data(specific_summary: Dict[str, Tuple[int, int]],
                         genome_bases: Dict[str, int]) -> Dict[
    str, Dict[str, int]]:
    """This func receive 2 dictionaries: one mapping ref header to the genome's
    length and the other mapping ref headers to the number of specific and unspecific
    Kmers they have. Func return a combined dictionary for each ref header mapped
    to all three values"""
    TOTAL_BASES: str = "total_bases"
    SPECIFIC_KMERS: str = "unique_kmers"
    UNSPECIFIC_KMERS: str = "multi_mapping_kmers"
    if len(specific_summary) != len(genome_bases):
        raise KeyError(
            "genome bases and specific summery dicts have different lengths and can't be combined")
    else:
        genome_summary: Dict[str, Dict[str, int]] = {}
        for key in specific_summary.keys():
            if key not in genome_bases:
                raise KeyError(
                    f"ref_header {key} does not exist in genome_bases dictionary")
            else:
                specific_kmers, unspecific_kmers = specific_summary[key]
                base_length: int = genome_bases[key]
                if key in genome_summary:
                    raise KeyError(
                        f"key: {key} appear in both genome_bases and and specific_summary twice")
                else:
                    genome_summary[key] = {TOTAL_BASES: base_length,
                                           SPECIFIC_KMERS: specific_kmers,
                                           UNSPECIFIC_KMERS: unspecific_kmers}
    return genome_summary


def generate_gen_sum(ref_kmers_db: KmersDB,
                     genome_bases: Dict[str, int]) -> Dict[
    str, Dict[str, int]]:
    """ This func receive KmersDB object, and a dictionary of genome lengths dictionary
    mapping ref_headers to their lengths.
    Func generates a ref_summery dictionary mapping each ref header to another
    dictionary with 3 keys: "total_bases","unique_kmers","multi_mapping_kmers"
    a kmer is unique if it appears only in one genome in DB.
    """
    specific_summary: Dict[str, Tuple[int, int]] = {}
    for kmer_data, ref_locs_dict in ref_kmers_db.get_ref_db().items():
        is_unique_kmer: bool = len(ref_locs_dict) == 1
        # if Kmer is specific in a genome it has only one source meaning
        # only one key in the inner dictionary - meaning it appears only in one ref
        for ref_header, appearances in ref_locs_dict.items():
            if ref_header not in specific_summary:
                specific_summary[ref_header] = (0, 0)
            specific_count, unspecific_count = specific_summary[ref_header]
            if is_unique_kmer:
                specific_summary[ref_header] = (
                    specific_count + len(appearances), unspecific_count)
            else:
                specific_summary[ref_header] = (
                    specific_count, unspecific_count + len(appearances))
    return combine_gen_sum_data(specific_summary, genome_bases)


def execute_dumpref_and_reference(genome_file_path: str, k_size: int,
                                  filter_sim: bool, sim_thresh: float) -> None:
    """This func receive a genome_file path, a Kmer_size and creates a reference
    Kmers DB and print the data of the DB like execute_dumpref()"""
    # building a Kmer_ref_DB:
    similarity: Dict[str, Dict[str, Any]] = {}
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_file_path)
    tup = build_ref_kmers_db(refs_dict, k_size)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    if filter_sim:
        similarity = filter_ref_kmers_db(refs_dict, refs_kmers_db,
                                         genomes_bases, sim_thresh)
    genome_summary: Dict[str, Dict[str, int]] = generate_gen_sum(refs_kmers_db,
                                                                 genomes_bases)
    ref_DB_dump = {"Kmers": refs_kmers_db.get_ref_db(),
                   "Summary": genome_summary}
    json_object = json.dumps(ref_DB_dump)
    print(json_object)
    if filter_sim:
        similarity_summary = {"Similarity": similarity}
        sim_json = json.dumps(similarity_summary)
        print(sim_json)


def execute_dumpref(ref_kmers_db_path: str) -> None:
    """This func receive a path to a KDB file with the reference Kmers database
    and print to stdout(screen) a summary of the Data in json format:
    parms printed:
    'Kmers' - for each Kmer a Dict of the references in which it appears and the
    location in each reference.
    'genome-summary'-
    for each genome_reference in refs_kmers_DB:
        1. total_bases: Total length of this genome's sequence
        2. unique_kmers: Number of k-mers that appear ONLY in this genome
        3. multi_mapping_kmers: Number of k-mers that appear in this genome AND
        at least one other genome
    """
    tup: Tuple[KmersDB, Dict[str, int]] = load_kdb_file(ref_kmers_db_path)
    ref_kmers_db: KmersDB = tup[0]
    genome_bases: Dict[str, int] = tup[1]
    genome_summary: Dict[str, Dict[str, int]] = generate_gen_sum(
        ref_kmers_db, genome_bases)
    ref_DB_dump = {"Kmers": ref_kmers_db.get_ref_db(),
                   "Summary": genome_summary}
    json_object: str = json.dumps(ref_DB_dump)
    print(json_object)
    similarity_summary = ref_kmers_db.get_similarity()
    if similarity_summary != {}:
        sim_json = json.dumps(similarity_summary)
        print(sim_json)


def execute_align_and_reference(genome_file_path: str, k_size: int,
                                aln_output_path: str, reads_path: str,
                                unique_thresh: int, ambig_thresh: int,
                                mrq: int, mkq: int, mg: int, rev_comp: bool,
                                coverage: bool, full_cov: bool,
                                min_cov: int,
                                gen_to_cov: Union[List[str], None]) -> None:
    """This func execute pseudo alignment on reads mapping it to ref_db, the
    output of the mapping process will be saved into an ALN file. parameters
    m and p are thresholds for minimum unique K-mers and minimum ambiguous K-mers
    for the mapping process logic. By default, m=p=1"""
    # building a Kmer_ref_DB:
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_file_path)
    tup = build_ref_kmers_db(refs_dict, k_size)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    filtered_q_reads: int = 0
    if mrq is not None:
        reads, filtered_q_reads = filter_lq_reads(reads, mrq)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(refs_dict.keys()))
    aln.set_ref_kmers_db(refs_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k_size, mkq)
    filtered_hr_kmers: int = aln.map_reads(unique_thresh, ambig_thresh, mg,
                                           rev_comp)
    aln.set_quality_stats(
        [mrq, mkq, mg, filtered_q_reads, filtered_q_kmers, filtered_hr_kmers])
    cov_dict: Dict[str, Tuple[Counter, Counter]] = {}
    if coverage:
        genomes_to_cover = list(
            genomes_bases.keys()) if not gen_to_cov or full_cov else gen_to_cov
        cov_dict = calc_coverage(aln, genomes_to_cover, genomes_bases)
        cov_stats_sum: Dict[str, Dict[str, Any]] = generate_cov_stats_sum(
            cov_dict, min_cov)
        if full_cov:
            det_pos_cov: Dict[str, Dict[str, Any]] = generate_det_pos_cov(
                cov_dict)
            dump_coverage = {"Coverage": cov_stats_sum, "Details": det_pos_cov}
        else:
            dump_coverage = {"Coverage": cov_stats_sum}
        aln.set_coverage_sum(dump_coverage)
    dump_alignment_results(aln, aln_output_path)


def execute_align(ref_db_path: str, aln_output_path: str, reads_path: str,
                  unique_thresh: int, ambig_thresh: int, mrq: int, mkq: int,
                  mg: int, rev_comp: bool, coverage: bool, full_cov: bool,
                  min_cov: int, gen_to_cov: Union[List[str], None]) -> None:
    """This func execute pseudo alignment on reads mapping it to ref_db, the
     output of the mapping process will be saved into an ALN file. parameters
     m and p are thresholds for minimum unique K-mers and minimum ambiguous K-mers
     for the mapping process logic. By default, m=p=1"""
    tup: Tuple[KmersDB, Dict[str, int]] = load_kdb_file(ref_db_path)
    ref_kmers_db: KmersDB = tup[0]
    genome_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    filtered_q_reads: int = 0
    if mrq is not None:
        reads, filtered_q_reads = filter_lq_reads(reads, mrq)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genome_bases.keys()))
    aln.set_ref_kmers_db(ref_kmers_db)
    k_size: int = find_kmer_size(ref_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k_size, mkq)
    filtered_hr_kmers: int = aln.map_reads(unique_thresh, ambig_thresh, mg,
                                           rev_comp)
    aln.set_quality_stats(
        [mrq, mkq, mg, filtered_q_reads, filtered_q_kmers, filtered_hr_kmers])
    cov_dict: Dict[str, Tuple[Counter, Counter]] = {}
    if coverage:
        genomes_to_cover = list(
            genome_bases.keys()) if not gen_to_cov or full_cov else gen_to_cov
        cov_dict = calc_coverage(aln, genomes_to_cover, genome_bases)
        cov_stats_sum: Dict[str, Dict[str, Any]] = generate_cov_stats_sum(
            cov_dict, min_cov)
        if full_cov:
            det_pos_cov: Dict[str, Dict[str, Any]] = generate_det_pos_cov(
                cov_dict)
            dump_coverage = {"Coverage": cov_stats_sum, "Details": det_pos_cov}
        else:
            dump_coverage = {"Coverage": cov_stats_sum}
        aln.set_coverage_sum(dump_coverage)
    dump_alignment_results(aln, aln_output_path)


def find_kmer_size(ref_kmers_db: KmersDB) -> int:
    """This func receive a Kmers_DB dictionary and find the size of Kmers in this
    which each read contain"""
    return ref_kmers_db.get_kmer_size()


def extract_gen_map_sum_rev(reads: Dict[str, Read], refs: List[str]) -> Dict[
    str, Dict[str, int]]:
    """This func receive a dictionary of mapped reads and for each read extract its
    sources and add them to a count for each genome in references. Func creates
    for each genome a 'forward' and 'reverse' key with all the unique and
    ambiguous reads mapped to it"""
    gen_map_sum: Dict[str, Dict[str, int]] = {}
    FORWARD, REV = "_F", "_R"
    REV_ORIENT = "Reverse"
    UNIQUE, AMBIG = "unique_reads", "ambiguous_reads"
    for ref_header in refs:
        gen_map_sum[ref_header + FORWARD] = {UNIQUE: 0, AMBIG: 0}
        gen_map_sum[ref_header + REV] = {UNIQUE: 0, AMBIG: 0}
    for read_header, read in reads.items():
        sources: List[str] = read.get_read_sources()
        if len(sources) == 1:
            source: str = sources[0]
            if read.get_orient() == REV_ORIENT:
                gen_map_sum[source + REV][UNIQUE] += 1
            else:
                gen_map_sum[source + FORWARD][UNIQUE] += 1
        else:
            if read.get_orient() == REV_ORIENT:
                for source in sources:
                    gen_map_sum[source + REV][AMBIG] += 1
            else:
                for source in sources:
                    gen_map_sum[source + FORWARD][AMBIG] += 1
    return gen_map_sum


def extract_gen_map_sum(reads: Dict[str, Read], refs: List[str]) -> Dict[
    str, Dict[str, int]]:
    """This func receive a dictionary of mapped reads and for each read extract its
    sources and add them to a count for each genome in references"""
    gen_map_sum: Dict[str, Dict[str, int]] = {}
    UNIQUE, AMBIG = "unique_reads", "ambiguous_reads"
    for ref_header in refs:
        gen_map_sum[ref_header] = {UNIQUE: 0, AMBIG: 0}
    for read_header, read in reads.items():
        sources: List[str] = read.get_read_sources()
        if len(sources) == 1:
            source: str = sources[0]
            if source in gen_map_sum:
                gen_map_sum[source][UNIQUE] += 1
            else:
                raise KeyError(
                    f"reads and references in Alignment object doesnt "
                    f"don't match. {source} doesnt exist in references")
        else:
            for source in sources:
                gen_map_sum[source][AMBIG] += 1
    return gen_map_sum


def execute_dumpalign(aln_results_path) -> None:
    """This func receive an alignment result ALN file output path, and load it,
    and print it. Func print the following parameters:
     {"Statistics":reads-stats, "Summary": genome-mapping-summary}"""
    MRQ, MKQ, MG, FILTERED_READS, FILTERED_Q_KMER, FILTERED_HR_KMERS = 0, 1, 2, 3, 4, 5
    aln: Alignment = load_aln_file(aln_results_path)
    reads_stats: Dict[str, int] = aln.get_reads_stats()
    gen_map_sum: Dict[str, Dict[str, int]] = extract_gen_map_sum(
        aln.get_reads(), aln.get_references())
    quality_stats = aln.get_quality_stats()
    if not quality_stats:
        dump_align: Dict[str, Dict[str, Any]] = {"Statistics": reads_stats,
                                                 "Summary": gen_map_sum}
    else:
        mrq, mkq, mg = quality_stats[MRQ], quality_stats[MKQ], quality_stats[
            MG]
        filtered_reads, filtered_q_kmers, filtered_hr_kmers = quality_stats[
            FILTERED_READS], quality_stats[FILTERED_Q_KMER], quality_stats[
            FILTERED_HR_KMERS]
        dump_align = generate_dump_align(aln, reads_stats, mrq, mkq, mg,
                                         filtered_reads,
                                         filtered_q_kmers, filtered_hr_kmers,
                                         False)
    json_object: str = json.dumps(dump_align)
    print(json_object)
    coverage_sum = aln.get_coverage_sum()
    if coverage_sum != {}:
        json_cov: str = json.dumps(coverage_sum)
        print(json_cov)


def add_filtered_stats(reads_stats: Dict[str, int], mrq, mkq, mg,
                       filtered_q_reads, filtered_q_kmers,
                       filtered_hr_kmers) -> None:
    """This func receive filtered objects data from this mapping process and
    add it to the reads_stats dictionary. Func change the original dictionary"""
    F_Q_READS, F_Q_KMERS, F_HR_KMERS = "filtered_quality_reads", "filtered_quality_kmers", "filtered_hr_kmers"
    if mrq is not None:
        reads_stats[F_Q_READS] = filtered_q_reads
    if mkq is not None:
        reads_stats[F_Q_KMERS] = filtered_q_kmers
    if mg is not None:
        reads_stats[F_HR_KMERS] = filtered_hr_kmers


def calc_coverage(aln: Alignment, gen_to_cover: List[str],
                  genome_bases: Dict[str, int]) -> Dict[
    str, Tuple[Counter, Counter]]:
    """This func receive an alignment object, genome to calc coverage for and
     genome bases lengths. Func initialize for each genome requested 2 counters
     for unique and ambiguous reads mapped to it and count how many Kmers cover
     each index of the original reference sequence data, and return a Dictionary of
    {ref_header, (unique_cov_counter, ambig_cov_counter)}.
    a kmer Uniquely covers a position if it belongs to a unique read and ambig
    cover a position if it belongs to an ambiguous read """
    COUNT_UNIQUE, COUNT_AMBIG = 0, 1
    cov_dict: Dict[str, Tuple[Counter, Counter]] = {}
    reads_kmers_db: Dict[str, List[Kmer]] = aln.get_reads_kmers_db()
    reads: Dict[str, Read] = aln.get_reads()
    ref_kmers_db: KmersDB = aln.get_ref_kmers_DB()
    for ref_header in gen_to_cover:
        # initialize for each genome requested to_cover 2 Counters:
        # both in the length of the genome and contain 0 for each position
        gen_length: int = genome_bases[ref_header]
        unique_cov: Counter = Counter({i: 0 for i in range(gen_length)})
        ambig_cov: Counter = Counter({i: 0 for i in range(gen_length)})
        cov_dict[ref_header] = (unique_cov, ambig_cov)

    for read_header, read in reads.items():
        if any(gen in read.get_read_sources() for gen in gen_to_cover):
            # checking if one of the read's sources is requested, if not - pass
            read_unique_counter: Counter = Counter()
            read_ambig_counter: Counter = Counter()
            for kmer in reads_kmers_db[read_header]:
                try:
                    kmer_locs_in_gens: Dict[
                        str, List[int]] = ref_kmers_db.get_value(
                        kmer.get_data())
                except KeyError:
                    continue
                for source_gen, locs_in_source in kmer_locs_in_gens.items():
                    for loc in locs_in_source:
                        # kmer may appear few times in the same genome
                        if source_gen in gen_to_cover:
                            for i in range(loc, loc + len(kmer.get_data())):
                                # add 1 to all positions: loc is the starting index
                                if read.get_status() == "Unique":
                                    if read_unique_counter[i] == 0:
                                        read_unique_counter[i] = 1
                                else:
                                    if read_ambig_counter[i] == 0:
                                        read_ambig_counter[i] = 1
            for source_gen in gen_to_cover:
                # adding reads' coverage to genome total coverage
                cov_dict[source_gen][COUNT_UNIQUE].update(read_unique_counter)
                cov_dict[source_gen][COUNT_AMBIG].update(read_ambig_counter)
    return cov_dict


def count_covered(unique_cov: Counter, ambig_cov: Counter, min_cov: int) -> \
        Tuple[int, int]:
    """This func receive two counters per genome and return the number of uniquely
    and ambiguously covered bases out of genome's length"""
    unique_cov_bases: int = 0
    ambig_cov_bases: int = 0
    for uq, am in zip(unique_cov.values(), ambig_cov.values()):
        if uq >= min_cov:
            unique_cov_bases += 1
        if am >= min_cov:
            ambig_cov_bases += 1
    return unique_cov_bases, ambig_cov_bases


def generate_cov_stats_sum(cov_dict: Dict[str, Tuple[Counter, Counter]],
                           min_cov: int):
    """This func receive a coverage dict for each requested genome and generates
    a coverage statistics summary containing: uniquely covered bases, ambiguously
     covered bases, mean_cov_unique and mean_cov_ambig"""
    COUNT_UNIQUE, COUNT_AMBIG = 0, 1
    UNIQUE, AMBIG = "covered_bases_unique", "covered_bases_ambiguous"
    MEAN_UNIQUE, MEAN_AMBIG = "mean_coverage_unique", "mean_coverage_ambiguous"
    cov_stat_sum: Dict[str, Dict[str, Any]] = {}
    for ref_header, tup_counters in cov_dict.items():
        uniquely_cov: int = sum(tup_counters[COUNT_UNIQUE].values())
        ambiguously_cov: int = sum(tup_counters[COUNT_AMBIG].values())
        gen_length: int = len(tup_counters[COUNT_AMBIG].keys())
        unique_cov_bases, ambig_cov_bases = count_covered(
            tup_counters[COUNT_UNIQUE], tup_counters[COUNT_AMBIG], min_cov)
        cov_stat_sum[ref_header] = {UNIQUE: unique_cov_bases,
                                    AMBIG: ambig_cov_bases,
                                    MEAN_UNIQUE: round(
                                        uniquely_cov / gen_length, ndigits=1),
                                    MEAN_AMBIG: round(
                                        ambiguously_cov / gen_length,
                                        ndigits=1)}
    return cov_stat_sum


def generate_det_pos_cov(cov_dict: Dict[str, Tuple[Counter, Counter]]) -> Dict[
    str, Dict[str, List[int]]]:
    """this func receive a coverage dict and return a dictionary of all genomes
    requested to cover and the cover: unique and ambiguous per each position in
    the genome"""
    COUNT_UNIQUE, COUNT_AMBIG = 0, 1
    INDEX, COUNT = 0, 1
    UNIQUE, AMBIG = "unique_cov", "ambiguous_cov"
    det_pos_cov: Dict[str, Dict[str, List[int]]] = {}
    for ref_header, count_tup in cov_dict.items():
        uq_counter = count_tup[COUNT_UNIQUE]
        am_counter = count_tup[COUNT_AMBIG]
        lst_unique_tups = sorted(uq_counter.items(), key=lambda x: x[INDEX])
        lst_ambig_tups = sorted(am_counter.items(), key=lambda x: x[INDEX])
        det_pos_cov[ref_header] = {
            UNIQUE: [tup[COUNT] for tup in lst_unique_tups],
            AMBIG: [tup[COUNT] for tup in lst_ambig_tups]}
    return det_pos_cov


def generate_dump_align(aln: Alignment, reads_stats: Dict[str, int], mrq: int,
                        mkq: int, mg: int, filtered_q_reads: int,
                        filtered_q_kmers: int,
                        filtered_hr_kmers: int, rev_comp: bool) -> Dict[
    str, Dict[str, Any]]:
    """This func generates a dumpalign summary considering reverse and quality
    extensions. Func return the finalized dictionary to be converted to json"""
    gen_map_sum: Dict[str, Dict[str, int]] = {}
    quality_stats: List[int] = aln.get_quality_stats()
    if any(x is not None for x in [mrq, mkq, mg]):
        add_filtered_stats(reads_stats, mrq, mkq, mg, filtered_q_reads,
                           filtered_q_kmers, filtered_hr_kmers)
    elif quality_stats:
        # quality_stats contain data from a previous align task
        MRQ, MKQ, MG, FILTERED_READS, FILTERED_Q_KMER, FILTERED_HR_KMERS = 0, 1, 2, 3, 4, 5
        add_filtered_stats(reads_stats, quality_stats[MRQ], quality_stats[MKQ],
                           quality_stats[MG], quality_stats[FILTERED_READS],
                           quality_stats[FILTERED_Q_KMER],
                           quality_stats[FILTERED_HR_KMERS])
    if not rev_comp:
        gen_map_sum = extract_gen_map_sum(aln.get_reads(),
                                          aln.get_references())
    else:
        gen_map_sum = extract_gen_map_sum_rev(aln.get_reads(),
                                              aln.get_references())
    dump_align: Dict[str, Dict[str, Any]] = {"Statistics": reads_stats,
                                             "Summary": gen_map_sum}
    return dump_align


def execute_align_and_dump(kdb_file_path: str, reads_path: str,
                           unique_thresh: int, ambig_thresh: int, mrq: int,
                           mkq: int, mg: int, rev_comp: bool, coverage: bool,
                           full_cov: bool, min_cov: int,
                           gen_to_cov: Union[List[str], None]) -> None:
    """This func execute align mapping process and then dump the results in json
    format like execute_dumpalign()"""
    filtered_q_reads: int = 0
    tup: Tuple[KmersDB, Dict[str, int]] = load_kdb_file(kdb_file_path)
    ref_kmers_db: KmersDB = tup[0]
    genome_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    if mrq is not None:
        reads, filtered_q_reads = filter_lq_reads(reads, mrq)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genome_bases.keys()))
    aln.set_ref_kmers_db(ref_kmers_db)
    k_size: int = find_kmer_size(ref_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k_size, mkq)
    filtered_hr_kmers: int = aln.map_reads(unique_thresh, ambig_thresh, mg,
                                           rev_comp)
    reads_stats: Dict[str, int] = aln.get_reads_stats()
    dump_align = generate_dump_align(aln, reads_stats, mrq, mkq, mg,
                                     filtered_q_reads, filtered_q_kmers,
                                     filtered_hr_kmers, rev_comp)
    json_object: str = json.dumps(dump_align)
    print(json_object)
    cov_dict: Dict[str, Tuple[Counter, Counter]] = {}
    coverage_sum = aln.get_coverage_sum()
    if coverage:
        genomes_to_cover = list(
            genome_bases.keys()) if not gen_to_cov or full_cov else gen_to_cov
        cov_dict = calc_coverage(aln, genomes_to_cover, genome_bases)
        cov_stats_sum: Dict[str, Dict[str, Any]] = generate_cov_stats_sum(
            cov_dict, min_cov)
        if full_cov:
            det_pos_cov: Dict[str, Dict[str, Any]] = generate_det_pos_cov(
                cov_dict)
            dump_coverage = {"Coverage": cov_stats_sum,
                             "Details": det_pos_cov}
        else:
            dump_coverage = {"Coverage": cov_stats_sum}
        json_cov = json.dumps(dump_coverage)
        print(json_cov)
    elif coverage_sum != {}:
        json_cov = json.dumps(coverage_sum)
        print(json_cov)


def execute_build_ref_and_align_and_dump(genome_file_path: str,
                                         k_size: int,
                                         reads_input_path: str,
                                         unique_thresh: int,
                                         ambig_thresh: int,
                                         mrq: int, mkq: int, mg: int,
                                         rev_comp: bool,
                                         coverage: bool, full_cov: bool,
                                         min_cov: int,
                                         gen_to_cov: Union[List[str], None]
                                         ) -> None:
    """This func receive a genome reference FASTA file, k_size and create a
    reference DB object. The load the reads from FASTQ file and create Alignment
    object and map the reads. Then, func print the mapping results in Json format"""
    filtered_q_reads: int = 0
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_file_path)
    tup = build_ref_kmers_db(refs_dict, k_size)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_input_path)
    if mrq is not None:
        reads, filtered_q_reads = filter_lq_reads(reads, mrq)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genomes_bases.keys()))
    aln.set_ref_kmers_db(refs_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k_size, mkq)
    filtered_hr_kmers: int = aln.map_reads(unique_thresh, ambig_thresh, mg,
                                           rev_comp)
    reads_stats: Dict[str, int] = aln.get_reads_stats()
    dump_align = generate_dump_align(aln, reads_stats, mrq, mkq, mg,
                                     filtered_q_reads, filtered_q_kmers,
                                     filtered_hr_kmers, rev_comp)
    json_object: str = json.dumps(dump_align)
    print(json_object)
    cov_dict: Dict[str, Tuple[Counter, Counter]] = {}
    coverage_sum = aln.get_coverage_sum()
    if coverage:
        genomes_to_cover = list(
            genomes_bases.keys()) if not gen_to_cov or full_cov else gen_to_cov
        cov_dict = calc_coverage(aln, genomes_to_cover, genomes_bases)
        cov_stats_sum: Dict[str, Dict[str, Any]] = generate_cov_stats_sum(
            cov_dict, min_cov)
        if full_cov:
            det_pos_cov: Dict[str, Dict[str, Any]] = generate_det_pos_cov(
                cov_dict)
            dump_coverage = {"Coverage": cov_stats_sum,
                             "Details": det_pos_cov}
        else:
            dump_coverage = {"Coverage": cov_stats_sum}
        json_cov = json.dumps(dump_coverage)
        print(json_cov)
    elif coverage_sum != {}:
        json_cov = json.dumps(coverage_sum)
        print(json_cov)


def check_valid_paths(file_path: str, valid_endings: List[str]) -> bool:
    """this func receive a file path and a list of valid endings for this file,
    if file doesn't have valid ending, or it doesn't exist, func return False and
    True otherwise"""
    full_path = os.path.abspath(file_path)
    if not os.path.isfile(full_path):
        return False
    for ending in valid_endings:
        if file_path.lower().endswith(ending.lower()):
            return True
    return False


def check_output_file(file_path: str,
                      required_extension: List[str]) -> bool:
    """This"""
    for req_ending in required_extension:
        if not file_path.lower().endswith(req_ending.lower()):
            raise ValueError(
                f"The file must end with '{required_extension}'")
    directory = os.path.dirname(file_path) or "."
    if not os.access(directory, os.W_OK):
        raise ValueError(f"The directory '{directory}' is not writable")
    return True


def check_valid_flag_combination(args) -> bool:
    """This func receive a parser with agrs from user and check if the received
     args are a valid combination. If not func raises matching Exception"""
    VALID_TASKS = ["reference", "dumpref", "align", "dumpalign"]
    if args.task not in VALID_TASKS:
        raise ValueError(f"Invalid command: {args.task}. Try again!")
    if args.task == "reference":
        return validate_reference(args)
    elif args.task == "dumpref":
        return validate_dumpref(args)
    elif args.task == "align":
        return validate_align(args)
    elif args.task == "dumpalign":
        return validate_dumpalign(args)
    return False


def validate_reference(args) -> bool:
    """This func validate arguments combination for the 'reference' task"""
    KDB_FILE = [".kdb"]
    GENOME_FILE = [".fa", ".fa.gz"]

    if not (20 <= args.kmer_size <= 31):
        raise ValueError(
            f"Kmer size should be 20-31, not {args.kmer_size}.")
    if not check_valid_paths(args.genomefile, GENOME_FILE):
        raise ValueError(f"Invalid genome file path: {args.genomefile}")
    if not check_output_file(args.referencefile, KDB_FILE):
        raise ValueError(
            f"Invalid reference file output path: {args.referencefile}")
    return True


def validate_dumpref(args) -> bool:
    """This func validate arguments combination for the 'dumpref' task."""
    KDB_FILE = [".kdb"]
    GENOME_FILE = [".fa", ".fa.gz"]
    if args.referencefile is not None:
        if args.genomefile is not None or args.kmer_size is not None:
            raise ValueError(
                "Flags (-g, -k) and (-r) can't be used together for 'dumpref'.")
        if not check_valid_paths(args.referencefile, KDB_FILE):
            raise ValueError(
                f"Invalid reference file path: {args.referencefile}")
    elif args.genomefile is not None and args.kmer_size is not None:
        if not (20 <= args.kmer_size <= 31):
            raise ValueError(
                f"Kmer size should be 20-31, not {args.kmer_size}.")
        if not check_valid_paths(args.genomefile, GENOME_FILE):
            raise ValueError(
                f"Invalid genome file path: {args.genomefile}")
    else:
        raise ValueError(
            "Either (-r) or both (-g, -k) must be provided for 'dumpref'.")
    return True


def validate_align(args) -> bool:
    """This func validate arguments combination for the 'align' task."""
    KDB_FILE = [".kdb"]
    GENOME_FILE = [".fa", ".fa.gz"]
    ALN_FILE = [".aln"]
    READS_FILE = [".fq", ".fq.gz"]
    if args.referencefile is not None:
        if args.genomefile or args.kmer_size:
            raise ValueError(
                "Flags (-g, -k) and (-r) can't be used together for 'align'.")
        if not check_valid_paths(args.referencefile, KDB_FILE):
            raise ValueError(
                f"Invalid reference file path: {args.referencefile}")
    elif args.genomefile is not None and args.kmer_size is not None:
        if not (20 <= args.kmer_size <= 31):
            raise ValueError(
                f"Kmer size should be 20-31, not {args.kmer_size}.")
        if not check_valid_paths(args.genomefile, GENOME_FILE):
            raise ValueError(
                f"Invalid genome file path: {args.genomefile}")
    if not check_output_file(args.alignfile, ALN_FILE):
        raise ValueError(
            f"Invalid alignment output file path: {args.alignfile}")
    if args.reads is None:
        raise ValueError(f"missing argument --reads")
    elif not check_valid_paths(args.reads, READS_FILE):
        raise ValueError(f"Invalid reads file path: {args.reads}")
    return True


def validate_dumpalign(args) -> bool:
    """This func validate arguments combination for the 'dumpalign' task."""
    KDB_FILE = [".kdb"]
    GENOME_FILE = [".fa", ".fa.gz"]
    READS_FILE = [".fq", ".fq.gz"]
    ALN_FILE = [".aln"]
    if args.alignfile is not None:
        # Align file cannot be combined with other flags
        if any(x is not None for x in
               [args.referencefile, args.genomefile, args.kmer_size,
                args.reads]):
            raise ValueError(
                "Task 'dumpalign' with flag (-a) should not receive any other parameter.")
        elif not check_valid_paths(args.alignfile, ALN_FILE):
            raise ValueError(f'ALN file path {args.alignfile} is invalid')
    elif args.referencefile:
        if args.genomefile is not None or args.kmer_size is not None:
            raise ValueError(
                "Flags (-g, -k) and (-r) can't be used together for 'dumpalign'.")
        if not check_valid_paths(args.referencefile, KDB_FILE):
            raise ValueError(
                f"Invalid reference file path: {args.referencefile}")
        if args.reads is None:
            raise ValueError(f"missing argument --reads")
        elif not check_valid_paths(args.reads, READS_FILE):
            raise ValueError(f"Invalid reads file path: {args.reads}")
    elif args.genomefile is not None and args.kmer_size is not None:
        if not (20 <= args.kmer_size <= 31):
            raise ValueError(
                f"Kmer size should be 20-31, not {args.kmer_size}.")
        if not check_valid_paths(args.genomefile, GENOME_FILE):
            raise ValueError(
                f"Invalid genome file path: {args.genomefile}")
        if args.reads is None:
            raise ValueError(f"missing argument --reads")
        elif not check_valid_paths(args.reads, READS_FILE):
            raise ValueError(f"Invalid reads file path: {args.reads}")
    return True


def extract_gen_to_cover(coverage: bool, full_cov: bool, genomes: str) -> \
        Union[List[str], None]:
    """This func extract genome to cover for the coverage extension"""
    if full_cov and coverage:
        return []
    if genomes is None and coverage:
        return []
    elif genomes is not None and coverage:
        gen_to_cov: List[str] = genomes.split(',')
        return gen_to_cov
    else:
        return None


def execute_command(args: argparse.Namespace) -> None:
    """This func execute command received from the command line"""
    try:
        check_valid_flag_combination(args)
        command: str = args.task
        genome_file_path: str = args.genomefile
        k_size: int = args.kmer_size
        kdb_output_path, kdb_input_path = args.referencefile, args.referencefile
        aln_output_path, aln_input_path = args.alignfile, args.alignfile
        reads_input_path: str = args.reads
        unique_thresh: int = args.unique_threshold
        ambig_thresh: int = args.ambiguous_threhold
        mrq: int = args.min_read_quality
        mkq: int = args.min_kmer_quality
        mg: int = args.max_genomes
        rev_comp: bool = args.reverse_complement
        filter_sim: bool = args.filter_similar
        sim_thresh: float = args.similarity_threshold
        full_cov: bool = args.full_coverage
        min_cov: int = args.min_coverage
        coverage: bool = args.coverage
        genomes_str: str = args.genomes
        gen_to_cov = extract_gen_to_cover(coverage, full_cov, genomes_str)
        if command == "reference":
            execute_reference(genome_file_path, k_size, kdb_output_path,
                              filter_sim, sim_thresh)
        elif command == "dumpref":
            if kdb_input_path is not None:
                execute_dumpref(kdb_input_path)
            else:
                execute_dumpref_and_reference(genome_file_path, k_size,
                                              filter_sim, sim_thresh)
        elif command == "align":
            if kdb_input_path is not None:
                execute_align(kdb_input_path, aln_output_path,
                              reads_input_path,
                              unique_thresh, ambig_thresh, mrq, mkq, mg,
                              rev_comp, coverage, full_cov,
                              min_cov, gen_to_cov)
            else:
                execute_align_and_reference(genome_file_path, k_size,
                                            aln_output_path,
                                            reads_input_path,
                                            unique_thresh, ambig_thresh,
                                            mrq,
                                            mkq, mg, rev_comp, coverage,
                                            full_cov,
                                            min_cov, gen_to_cov)
        elif command == "dumpalign":
            if kdb_input_path is None and genome_file_path is None and k_size is None:
                execute_dumpalign(aln_input_path)
            elif kdb_input_path is not None:
                execute_align_and_dump(kdb_input_path, reads_input_path,
                                       unique_thresh, ambig_thresh, mrq,
                                       mkq,
                                       mg, rev_comp, coverage, full_cov,
                                       min_cov, gen_to_cov)
            else:
                execute_build_ref_and_align_and_dump(genome_file_path,
                                                     k_size,
                                                     reads_input_path,
                                                     unique_thresh,
                                                     ambig_thresh, mrq,
                                                     mkq,
                                                     mg, rev_comp,
                                                     coverage,
                                                     full_cov, min_cov,
                                                     gen_to_cov)
    except ValueError as ve:
        print(ve)
        print(
            "Try again running the program with valid arguments. Use --help for help")


def main():
    args = readargs()
    try:
        execute_command(args)
    except Exception as ve:
        print(ve)


if __name__ == "__main__":
    main()
