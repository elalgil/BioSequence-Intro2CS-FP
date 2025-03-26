# FILE : final_project_test.py
# WRITER : Elal Gilboa , elal.gilboa , 323083188
# EXERCISE : intro2cs final project 2025
# DESCRIPTION: This program is a test file using pytest, monkeypatch and capsys,
# to test the main.py program in the final project.
# STUDENTS I DISCUSSED THE EXERCISE WITH: None
# WEB PAGES I USED: https://www.geeksforgeeks.org/monkey-patching-in-python-dynamic-behavior/
# NOTES: None
import os

###########################################################################
from main import *
import tempfile
import sys
import random
import pytest


############################################################################

def generate_ref_quality():
    """This func generate a FASTA file for test_quality() function"""
    ref_path = "test_quality.fa"
    ref_content = """>MGYG000004169_21|test1:taxid|51124
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
>MGYG000004169_22|test2:taxid|51315
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
>MGYG000004129_23|test3:taxid|51316
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGAAAAAA
ACTGGACACGGCCCA
>MGYG000004179_33|test4:taxid|51316
GACTGATATTGTACATGCTATAGTTCTACTGGCTGGATTTACTATTGCCGTTCCATTTGCC
ATTCATAACGCTGG
>MGYG000004119_36|test5:taxid|51316
GAGACGGAGAGGGTTATCAGGTCTTGATGGAACTACTGAGTACTGATCCGGATACCATTGG
AAACTTTATGTCCG"""
    with open(ref_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    return ref_path


def generate_reads_quality():
    """This func generate a FASTQ file for test_quality() function"""
    reads_path = "test_quality.fq"
    reads_content = """@MGYG000002096_86|kraken:taxid|3024_read_0
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCAGCTCCGGCTGACAGC
+
DC@@B>>A>BBA??B=@CB@==?A><@??=B>><==A??>=@><;;>?>=@;?>?99>=>>:;>:8>8;87;97:
@MGYG000004832_79|kraken:taxid|6050_read_1
GCCCATTGTCATGAGCATTATGGCGGCGGCCTATCCAGGTTTTGAACTGGACACGTGCCAGCTCCGGCTGACAGT
+
B@AC??BC>>?CA@B@?A@?==BB@<@A?=@>=A?<AAA?@A;@;;><?=<=;=;9;>:<?<98;<==987887<
@MGYG000002096_86|kraken:taxid|3024_read_2
ACCCGCATCCGGAATGCGGACGATGTGGAGAAGTTACTGGATCTGCCGGTTTTGGGGCAGATTCCGCGGGAAGAT
+
B@?CCD>C>BD>ACCB?A>@B>?A?@=BABB@@A>==;=>;<@==>==><<:??;9?;?=:9=<9<:<>;<;<87
@MGYG000000975_36|kraken:taxid|1603_read_3
CAGCACTATGACCATAAACGCTTCCATAAACAATTAATGCACCTTCAACTTCAGGTGTATATGAACTCCATTTAT
+
DDE@DB@@AA>BA?=?=B=>CA?CAB=??B<B@?><;;@;@@;;??@<;><?>?=>::<>?8<=;<=;9879:8;
@MGYG000001129_3|kraken:taxid|1782_read_4
CGCTGGAAGTGTCCATGACGGTACCAAGTGGCAGACGCGGTAGCTTGGCAAGTAATTTTTCTCGAATGGTTTTCA
+
@?BBC@A@>D@@BA@@=CBA@?CB?>>?=?<BA?A=@=>==><=;=<@;>?;:>>>9>?>9>=;8<;89<9::8:
@MGYG000004787_66|kraken:taxid|6005_read_5
GGGTCGCCCTCCACCGCCGCCAGCAAGGGCAGCACGTGCTCGTATTTGGTGCTCTTCATGCCGGAGGGCAGCAGC
+
@DDEEDDD@@A>@B@>>AB?C>=@=?@<A<?<<?;A@?<A?A?=@<?=:;<;=;9>>?9<=>8>;<:8>>=9;9<
@MGYG000004787_66|kraken:taxid|6005_read_6
ACGAAAATTCGTCACGCCGCTGGTGTTCCGGCTAAAATCAGCCCCGCTCACGCATCTTGGCGAAGCAGTCGTTCA
+
BD?CCBCD>>>ADB=CC@B>@>BAB>A><>=A@?<>==A<??>?<=@;;@=>9==>><9<98<:<>;88>=9::<
@MGYG000000125_5|kraken:taxid|343_read_7
AGACTCAGATTATTTCTGTACAGTTTGGACTACATCTGCGGCATCTTTTATATTATCATAGGCAGGAGCCAACAA
+
?DABEA@BADB?A>@==>CBCBA?@<=@B=>?>?A=?<;A<=@;><<;=:<?<=?9;?>:9:9:=<8=><;9:98
@MGYG000004169_19|kraken:taxid|5353_read_8
TCAGTACCTTGGCAGAAACTTTGACGAAGGTATCCGGAACAATCATGAAGTCTGACGGAGCAGCGGAGTTAGGCA
+
BBEBDC>D?>B?@D=B@=?A@CAA@@@AA>A=><?A<><;=?@;>@;@;<:9<<==;:;<9>98;9<;;8=;<77
@MGYG000002096_86|kraken:taxid|3024_read_9
GCGGATGGGGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGG
+
?EC@@?B?CB@D@CCBBCB?=?CAB?<>B?A>A@;@>A<??=?;<@;@<@;>?9:;<:<=98=9=;8;8>;9;=;
@MGYG000004169_19|kraken:taxid|5353_read_10
ATGAATTTTTGAGAAAAGCAGGTGAACAGGTTGACTGAAAATTGTAAGAAATATGCAGTAGAGCTGCATGAACTC
+
A@?CE>A?@?B>A>BA?>=>B>?C=>B@B=>??>><?<=<;@:<;?>:=:<?=>9:9:>;<>;=9>8<:=7;=:;
@MGYG000002096_86|kraken:taxid|3024_read_11
AGAAAATGGGCAACAATATAGGAAAGGATGAACACTGTTTGAAGAGAGGAGCACGCGGATGGGGAAAAATATGGA
+
?@E@@BBB>ABCB?CB????=>>?BBA?=@A?<?<A@@@;>A??;:;@?=@?9<999<9==>>9;==>8:888<;
@MGYG000003260_60|kraken:taxid|4334_read_12
GGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGGCGACGAGG
+
DBA?CBCABCA@A?@>AB==AB>?AB>@A=<<B=A?;=;;=?;=;;>>;@<>:><>;<=?><<<9;>9>;:;=77
@MGYG000002096_86|kraken:taxid|3024_read_13
GAGGGACTGGAAATTATCATTGCCCATCCCGAGCGATATCAATATGTGCCCTGGGGATCGCTCTGGACTGGCAGG
+
ABC@D@ABBCDD>B=>@C=@@>CCABB>A@@?@@=<@A=;;>@;;@=??@;?9=>??=;;:9:8:<<=;8=8<9:
@MGYG000001009_274|kraken:taxid|1639_read_14
CCTGTGGCCGGGCTTGGCGTACGAGCGCCTGGACGGTCTGCGCGAGTGGTTCTTTGGCGACTTTGAGGCCGAGCG
+
DC?DCABCDD@DA@B@B??A?C@CAAAA@AB<A@==A?;@A><@;;>?<;::;=?;=><9?=>=>9;9<9==:78
@MGYG000004832_79|kraken:taxid|6050_read_15
TGGCCCTGATACAGGCCGGTCAAAAACTCCGCCTCGGTGAGGTTGGGCACGATCACGTCGGCCACGCCGATGAGA
+
CDCC@@AD>AC>C?=>BB@?CB@=B=AA@B@BA<?<@;?<A<:==<>;:??:<9=>9=9;:>>9;98>8=<;<;9
@MGYG000000125_5|kraken:taxid|343_read_16
TCCAAAAGTGCATCCATCAACCCATAACTACCCTCCCGACAAGATGTCACAAATAGCCTTTAAATGTAGCGTTTA
+
EA?@C??DD@>CDB>@=C@>B@@A==A?A<B<@;=><=;=;A?=>>@<?@<>=<=:::;9<9>>==<>:8<8=;:
@MGYG000002096_86|kraken:taxid|3024_read_17
TGGCCGGATTTTTGGGGCTTGTCATTAGCCGCATCCTTGGCCTGGGGCTGTCCTGGGCGGTGCTTTGCATTGCCT
+
DA?C@@AC?CD?BDACCAC@B@@>B<AA<A@AB@=?A@;?>;>;??<><>;:=>??<?:999<=<8>8:9<=<;;
@MGYG000004832_79|kraken:taxid|6050_read_18
GATCTGGGCCAGACGGCTGTCGTCCCGCCCGGCGGACAGGGTGTCCACCAGCTGGGCGAAGGGCCCCTGGGGGTC
+
ECB??AD@?CDBBBAB=A?BBC??=B<AA>>@=><A<<=>A=<@==;<;:@<>?;>:==??989:=:=8:8;9=:
@MGYG000000975_36|kraken:taxid|1603_read_19
CTAATCATAACACTATCTTCGATTTCACAAATAAATCCATAATGTGAACCTAAATCAATAACTTCAACAACATTT
+
@@BCC?D>>B>@BCC?=@A>CA?>@<<=<@AB?A@A=;>;<A??:?>@;=<>>>?;;<::;9:89>8:<<;;:;8
@MGYG000002096_86|kraken:taxid|3024_read_20
GCCCATTGTTAAGAGCATTATGGCGGCGGCCTATCCCGGTTTTGTACTGGACACGGCCCATCTCCGGCTGACTGC
+
EE?@?D>@A@@A?>CBB@C?@=@C?@==<A=>A<?=@A;A??@<;@?:?=>:99?<;:=:=;=8>:;9;;;8;<:
@MGYG000004169_19|kraken:taxid|5353_read_21
CGAGAAGGAAGACGGTACAGTGGTCAACATGGGGGATCAGGTGCTGGATAATCTTCATGAGAATACCAGTCTGGG
+
AD?@D@BADC@>BC=@@BB@A?@?==A@=<>?A>A=@?<<><;=:;??>=:?<?;:=?=9<::9:<>8=978=9;
@MGYG000004832_79|kraken:taxid|6050_read_22
ATGGTCAGCACATAGACCCGCCCCTCCTCCCGCTGGATGAGGGAGGTACCCAGGTCCACCAGGTTCTGCGTGTTG
+
DAD@CB?BC>BD?DC?>A=?CA=>?B<A=<?BA<<;=>>?;>>:=@;@:>=?>:>?>=9=>:>==:=9;;<87;<
@MGYG000000125_5|kraken:taxid|343_read_23
GTTTTTTGAGAAGAGATTCTCCCATTGCCGAAGATCCCATTCCTGCATCGCAGGCAAACACTACTAGATTAATAT
+
EE?ACC?CC@DBA>>=A==?CAAA>??>B=<=AA<;=@;>=?:=:@===:?>;<?;=:9>?=>:<=9:<;<9<:9
@MGYG000004787_66|kraken:taxid|6005_read_24
CAATGCCGAGGATTTCGCCGCCGTACGCCGAAAACCCGACCTTGTCAAGCAGCGTGACGCCGTCGTCGCTGATGC
+
?@CED>BBBD@?CCBAA@ABBCB?><BA@<>?=@;<<>;<A=@;<;:@;>=9;?=>::::9=>9<;9><999879
@MGYG000000975_36|kraken:taxid|1603_read_25
CTAGGAATATTAAATAATTTTAATAATTCTTCATCATAGTCAAGTGTGTTGATATTAAAAAGCATTGTTCTAGAA
+
E?B?A>AC@C@BBCCA>=?A=>>@AB<<>@<A><>A><=?AA:;=:@@>@=<<>:9=:?;>9:><:<><=7:8:=
@MGYG000000125_5|kraken:taxid|343_read_26
GACTGATATTGTACATGCTATAGTTCTACTGGCTGGATTTACTATTGCCGTTCCATTTGCCATTCATAACGCTGG
+
@DA@@@D@>>BDDBAC=CBCB>CA>==@?BB<>A<@@<;@A><=><;<<>;9?<=?=;=<==8>98>>:99=;7;
@MGYG000004169_19|kraken:taxid|5353_read_27
GAGACGGAGAGGGTTATCAGGTCTTGATGGAACTACTGAGTACTGATCCGGATACCATTGGAAACTTTATGTCCG
+
ADACBCACDDDCAB@CA>A??=?@??@?@BA?A@A<@@;A>A=>>:;:@:=9;><=><>99;<>>==99>;<7;=
@MGYG000002720_43|kraken:taxid|3759_read_28
AGTGCGGAGATCAGCAGGGAACCCACAAGGTTCATCGCCAGCACGATGATGACCGCAACGACGACGGCGATCAGC
+
E@DBABCB>B>C?@?>?AA?AC?C=>@@=?@=<?=<<<A>@;;>;<?;>>::;9<>?9???;<9;<8:=87:8;:
@MGYG000002720_43|kraken:taxid|3759_read_29
TTTCCGCGATCCTGTGCTGTGCGCCCTCGATGGTCAGCACGCAGGGAAGTCCCAGCTCATCCACCTTCTTAGCAA
+
B?D@?CDCABD@D?@C?=CBBA>C=@>BA@B???=??A@A?><>@?:?;><?>9:>=>;?;;9;=:<8:9==<=<
"""
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    return reads_path


def test_quality_extension():
    """This func test quality extension on map_reads() function to see if low
    quality reads got filtered out. to do so func generates FASTA/Q files to
    create an alignment object and run map_reads()"""
    genome_path = generate_ref_quality()
    reads_path = generate_reads_quality()
    k = 31
    filtered_q_reads: int = 0
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_path)
    tup = build_ref_kmers_db(refs_dict, k)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    mrq = 29
    if mrq is not None:
        reads, filtered_q_reads = filter_lq_reads(reads, mrq)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genomes_bases.keys()))
    aln.set_ref_kmers_db(refs_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k, None)
    filtered_hr_kmers: int = aln.map_reads(1, 1, 2, False)
    assert filtered_q_reads == 10
    os.remove(genome_path)
    os.remove(reads_path)


def test_quality_highly_redundant():
    """This func test quality extension on map_reads() function to see if highly
     redundant Kmers got filtered out. to do so func generates FASTA/Q files to
    create an alignment object and run map_reads()"""
    genome_path = generate_ref_quality()
    reads_path = generate_reads_quality()
    k = 31
    filtered_q_reads: int = 0
    refs_dict: Dict[str, Reference] = generate_ref_dict(genome_path)
    tup = build_ref_kmers_db(refs_dict, k)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genomes_bases.keys()))
    aln.set_ref_kmers_db(refs_kmers_db)
    filtered_q_kmers: int = aln.build_reads_kmers_db(k, None)
    filtered_hr_kmers: int = aln.map_reads(1, 1, 2, False)
    assert filtered_hr_kmers == 45
    os.remove(genome_path)
    os.remove(reads_path)


def test_kmer_class():
    """This function test Kmer Class and it's attributes to see logic of the class
    remain and comparison between Kmer objects works correctly """
    kmer1 = Kmer("ATTC", [17])
    kmer2 = Kmer("GGGC", [12])
    kmer3 = Kmer("TTTC", [114])
    kmer4 = Kmer("AAA", [1])
    kmer5 = Kmer("ATTC", [17])
    assert kmer1.get_data() == "ATTC"
    assert len(kmer2) == 4
    assert kmer1 != kmer2
    assert kmer3.get_location() == [114]
    assert kmer1 != {1, 2, 3}
    assert kmer1 != kmer4
    assert kmer1 == kmer5

    def print_kmer(kmer: Kmer):
        return (f"K-mer id: {kmer.get_id()},"
                f" data: {kmer.get_data()}," +
                f" locations in read: {kmer.get_location()}," +
                f" matching reference sources: {kmer.get_sources()}," +
                f" specific K-mer: {kmer.get_specific_status()}")

    assert print_kmer(kmer1) == str(kmer1)
    assert str(kmer1) == kmer1.__repr__()


def test_read_class():
    """This function test Read class and it's attributes, func it testing
    comparison between Read object, extracting Kmers and printing Reads"""
    read1 = Read("read1", "AAAC", "!!!!")
    read2 = Read("read2", "CCCC", "????")
    read3 = Read("read3", "CCCC", "????")
    read4 = Read("read3", "CCCCT", "???!?")
    assert read2 == read3
    dict_kmers = {'AAA': [0], 'AAC': [1]}
    assert dict_kmers == read1.extract_kmers(3)
    assert read1 != {1, 2, 3}
    assert read1 != read4
    read5 = Read("read5", "AAANAAA", "IIIIIII")
    read_kmers = read5.extract_kmers(3)
    assert read_kmers == {"AAA": [0, 4]}
    assert read1.set_status("Unmapped") == True
    read1 = Read("MGYG000004169_19|kraken:taxid|5353_read_10",
                 "ATGAATTTTTGAGAAAAGCAGGTGAACAGGTTGACTGAAAATTGTAAGAAATATGCAGTAGAGCTGCATGAACTC",
                 "A@?CE>A?@?B>A>BA?>=>B>?C=>B@B=>??>><?<=<;@:<;?>:=:<?=>9:9:>;<>;=9>8<:=7;=:;")

    read2 = Read("MGYG000000125_5|kraken:taxid|343_read_7",
                 "AGACTCAGATTATTTCTGTACAGTTTGGACTACATCTGCGGCATCTTTTATATTATCATAGGCAGGAGCCAACAA",
                 "?DABEA@BADB?A>@==>CBCBA?@<=@B=>?>?A=?<;A<=@;><<;=:<?<=?9;?>:9:9:=<8=><;9:98")

    assert read1 != read2
    assert read1.get_data() == "ATGAATTTTTGAGAAAAGCAGGTGAACAGGTTGACTGAAAATTGTAAGAAATATGCAGTAGAGCTGCATGAACTC"
    assert read2.get_quality() == "?DABEA@BADB?A>@==>CBCBA?@<=@B=>?>?A=?<;A<=@;><<;=:<?<=?9;?>:9:9:=<8=><;9:98"
    dict_of_kmers = {'ATGAATTTTT': [0], 'TGAATTTTTG': [1], 'GAATTTTTGA': [2],
                     'AATTTTTGAG': [3], 'ATTTTTGAGA': [4], 'TTTTTGAGAA': [5],
                     'TTTTGAGAAA': [6], 'TTTGAGAAAA': [7], 'TTGAGAAAAG': [8],
                     'TGAGAAAAGC': [9], 'GAGAAAAGCA': [10], 'AGAAAAGCAG': [11],
                     'GAAAAGCAGG': [12], 'AAAAGCAGGT': [13],
                     'AAAGCAGGTG': [14], 'AAGCAGGTGA': [15],
                     'AGCAGGTGAA': [16], 'GCAGGTGAAC': [17],
                     'CAGGTGAACA': [18], 'AGGTGAACAG': [19],
                     'GGTGAACAGG': [20], 'GTGAACAGGT': [21],
                     'TGAACAGGTT': [22], 'GAACAGGTTG': [23],
                     'AACAGGTTGA': [24], 'ACAGGTTGAC': [25],
                     'CAGGTTGACT': [26], 'AGGTTGACTG': [27],
                     'GGTTGACTGA': [28], 'GTTGACTGAA': [29],
                     'TTGACTGAAA': [30], 'TGACTGAAAA': [31],
                     'GACTGAAAAT': [32], 'ACTGAAAATT': [33],
                     'CTGAAAATTG': [34], 'TGAAAATTGT': [35],
                     'GAAAATTGTA': [36], 'AAAATTGTAA': [37],
                     'AAATTGTAAG': [38], 'AATTGTAAGA': [39],
                     'ATTGTAAGAA': [40], 'TTGTAAGAAA': [41],
                     'TGTAAGAAAT': [42], 'GTAAGAAATA': [43],
                     'TAAGAAATAT': [44], 'AAGAAATATG': [45],
                     'AGAAATATGC': [46], 'GAAATATGCA': [47],
                     'AAATATGCAG': [48], 'AATATGCAGT': [49],
                     'ATATGCAGTA': [50], 'TATGCAGTAG': [51],
                     'ATGCAGTAGA': [52], 'TGCAGTAGAG': [53],
                     'GCAGTAGAGC': [54], 'CAGTAGAGCT': [55],
                     'AGTAGAGCTG': [56], 'GTAGAGCTGC': [57],
                     'TAGAGCTGCA': [58], 'AGAGCTGCAT': [59],
                     'GAGCTGCATG': [60], 'AGCTGCATGA': [61],
                     'GCTGCATGAA': [62], 'CTGCATGAAC': [63],
                     'TGCATGAACT': [64], 'GCATGAACTC': [65]}
    assert read1.extract_kmers(10) == dict_of_kmers


def test_kerms_db():
    """This func test KmerDB class and it's attributes to make sure inserting
    a new Reference to DB works correctly"""
    db = KmersDB()
    text1 = "TTTCATACGCACTTGCGGCATTCGAAGCGATGAATGTCGGAGCAAGCAGCAACGCAGCAA\
ACTTTACACGTAAAAAGCTTTGCGTTTTCAAATGCATAACACTCTTCTTTCTCATAAAAT\
TTGATGTCTTTCAACAGTGTGCAAGCCTGATCCAGGTCGGCATTGATTTTGTTATATACC\
TCAGCCACTGAGCTGCGGGCCAGTTCTTCGGCCACGGCTTTTTCACGGTAAGGAATACCA\
TCCTGGCTGTTGTCTGTTCCAGCTTTATACGGTTTGGCATAGAGCTGAACCAGATTGAAG"

    text2 = "TTTCATACGCACTTGCGGCATTCGAAGCGATGAATGTCGGAGCAAGCAGCAACGCAGCAA\
ATGTTTGCAAAACCGCACATCAGCTGGTCGCCTTCCATCAAGGCATGGCCCTGTGCTTCT\
GCTAAATCAATAGCAGTCTTAGCATATTGAGCTGCATTCGGATAATCCTGCATGGTCAGT\
TCCTGGCTGTTGTCTGTTCCAGCTTTATACGGTTTGGCATAGAGCTGAACCAGATTGAAG"

    ref1 = Reference("MGYG000002614_64|kraken:taxid|1201", text1)
    ref2 = Reference("MGYG000332214_62|kraken:taxid|3633", text2)
    db.inset_ref_data(ref1.get_header(), ref1.extract_kmers(30))
    db.inset_ref_data(ref2.get_header(), ref2.extract_kmers(30))
    inner_dict = db.get_value("TTTCATACGCACTTGCGGCATTCGAAGCGA")
    assert len(inner_dict) == 2, ("inserting references to DB should have "
                                  "matched the same string to 2 refs")


def test_map_reads1():
    """This func test function map_reads in Alignment Class to make sure it
    returns the expected result. In this example there are 2 Genomes and the read
    is uniquely mapped to genome snail. Func checks read got mapped Uniquely and
    that it's source is set to snail only"""
    align = Alignment()
    k = 4
    read = Read("@MGYG000004787_66|kraken:taxid|6005_read_103", "ATGGCTAT",
                "AT!!ABCD")
    align.set_reads({"@MGYG000004787_66|kraken:taxid|6005_read_103": read})
    align.build_reads_kmers_db(k, None)
    ref1 = Reference(">MGYG000002614_62|snail:taxid|3633", "ATGGCTA")
    ref2 = Reference(">MGYG000002614_62|bird:taxid|2211", "CTAT")
    align.set_references([">MGYG000002614_62|snail:taxid|3633",
                          ">MGYG000002614_62|bird:taxid|2211"])
    ref_dict = {ref1.get_header(): ref1, ref2.get_header(): ref2}
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    align.set_ref_kmers_db(ref_kmers_db)
    align.map_reads(1, 1, None, False)
    read = list(align.get_reads().values())[0]
    assert read.get_status() == "Unique"
    assert read.get_read_sources() == ['>MGYG000002614_62|snail:taxid|3633']


def test_map_reads2():
    """This func test function map_reads in Alignment Class to make sure it
    returns the expected result. In this example there are 4 Genomes and the read
    is ambiguously mapped to genome snail, dog, bird. Func checks read got
    mapped Ambiguous and that it's source is set to dog, snail, bird"""
    align = Alignment()
    k = 4
    read = Read("@MGYG087_66|kraken:taxid|6005_read_103", "ATGCCTAT",
                "AT!!ABCD")
    align.set_reads({"@MGYG000004787_66|kraken:taxid|6005_read_103": read})
    align.build_reads_kmers_db(k, None)
    ref1 = Reference(">MGYG014_62|snail:taxid|3633", "ATGCC")
    ref2 = Reference(">MGYG061_62|bird:taxid|2211", "GCCTA")
    ref3 = Reference(">MGYG061_41|dog:taxid|2213", "TGCCTAT")
    ref4 = Reference(">MGYG061_11|cat:taxid|2214", "ATGC")
    align.set_references([">MGYG014_62|snail:taxid|3633",
                          ">MGYG061_62|bird:taxid|2211",
                          ">MGYG061_41|dog:taxid|2213",
                          ">MGYG061_11|cat:taxid|2214"])
    ref_dict = {ref1.get_header(): ref1, ref2.get_header(): ref2,
                ref3.get_header(): ref3, ref4.get_header(): ref4}
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    align.set_ref_kmers_db(ref_kmers_db)
    align.map_reads(1, 1, None, False)
    read = list(align.get_reads().values())[0]
    assert read.get_status() == "Ambiguous"
    assert read.get_read_sources() == ['>MGYG061_41|dog:taxid|2213',
                                       '>MGYG014_62|snail:taxid|3633',
                                       '>MGYG061_62|bird:taxid|2211']


def test_map_reads3():
    """This func test function map_reads in Alignment Class to make sure it
    returns the expected result. In this example there are 5 Genomes and the read
    is ambiguously mapped to genome snail, dog, bird. In this case read originally
    mapped Unique but in the validation process resulted ambiguous.
    Func checks read got mapped Ambiguous and that it's source is set to dog, snail, bird"""
    align = Alignment()
    k = 4
    read = Read("@MGYG087_66|kraken:taxid|6005_read_103", "ATGCCGGGGCTAA",
                "AT!!ABCDI?III")
    align.set_reads({"@MGYG000004787_66|kraken:taxid|6005_read_103": read})
    align.build_reads_kmers_db(k, None)
    ref1 = Reference(">MGYG014_62|snail:taxid|3633", "ATGCCGGGG")
    ref2 = Reference(">MGYG061_62|bird:taxid|2211", "GCCGGGGCTAA")
    ref3 = Reference(">MGYG061_41|dog:taxid|2213", "CCGG")
    ref4 = Reference(">MGYG061_11|cat:taxid|2214", "GGGCTAA")
    ref5 = Reference(">MGYG061_34|owl:taxid|4455", "GCTAA")
    align.set_references([">MGYG014_62|snail:taxid|3633",
                          ">MGYG061_62|bird:taxid|2211",
                          ">MGYG061_41|dog:taxid|2213",
                          ">MGYG061_11|cat:taxid|2214",
                          ">MGYG061_34|owl:taxid|4455"])
    ref_dict = {ref1.get_header(): ref1, ref2.get_header(): ref2,
                ref3.get_header(): ref3, ref4.get_header(): ref4,
                ref5.get_header(): ref5}
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    align.set_ref_kmers_db(ref_kmers_db)
    align.map_reads(1, 1, None, False)
    read = list(align.get_reads().values())[0]
    assert read.get_status() == "Ambiguous"
    assert read.get_read_sources() == ['>MGYG014_62|snail:taxid|3633',
                                       '>MGYG061_62|bird:taxid|2211',
                                       '>MGYG061_11|cat:taxid|2214']


def test_load_dump_kdb():
    """Func test the loading of a KDB file and converting it to a KmerDB object.
    the loading results 2 pickled items from KDB file: KmerDB and genome bases.
    Func checks they are the same as original"""
    align = Alignment()
    k = 4
    read = Read("@MGYG087_66|kraken:taxid|6005_read_103", "ATGCCGGGGCTAA",
                "AT!!ABCDI?III")
    align.set_reads({"@MGYG000004787_66|kraken:taxid|6005_read_103": read})
    align.build_reads_kmers_db(k, None)
    ref1 = Reference(">MGYG014_62|snail:taxid|3633", "ATGCCGGGG")
    ref2 = Reference(">MGYG061_62|bird:taxid|2211", "GCCGGGGCTAA")
    ref3 = Reference(">MGYG061_41|dog:taxid|2213", "CCGG")
    ref4 = Reference(">MGYG061_11|cat:taxid|2214", "GGGCTAA")
    ref5 = Reference(">MGYG061_34|owl:taxid|4455", "GCTAA")
    align.set_references([">MGYG014_62|snail:taxid|3633",
                          ">MGYG061_62|bird:taxid|2211",
                          ">MGYG061_41|dog:taxid|2213",
                          ">MGYG061_11|cat:taxid|2214",
                          ">MGYG061_34|owl:taxid|4455"])
    ref_dict = {ref1.get_header(): ref1, ref2.get_header(): ref2,
                ref3.get_header(): ref3, ref4.get_header(): ref4,
                ref5.get_header(): ref5}
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    align.set_ref_kmers_db(ref_kmers_db)
    original_kmer_db: KmersDB = ref_kmers_db
    genome_bases = {">MGYG014_62|snail:taxid|3633": 9,
                    ">MGYG061_62|bird:taxid|2211": 11,
                    ">MGYG061_41|dog:taxid|2213": 4,
                    ">MGYG061_11|cat:taxid|2214": 7,
                    ">MGYG061_34|owl:taxid|4455": 5
                    }
    dump_ref_kmers_db(original_kmer_db, genome_bases, "reference_file.kdb")
    tup = load_kdb_file("reference_file.kdb")
    ref_kmers_db: KmersDB = tup[0]
    genome_bases1 = tup[1]
    assert original_kmer_db.__repr__() == ref_kmers_db.__repr__()
    assert genome_bases1 == genome_bases


def generate_random_fasta(filename: str, num_genomes: int):
    """This is a helper function to generate a randon FASTA file"""

    def gen_random_seq(length: int) -> str:
        return ''.join(random.choices("ACGT", k=length))

    def gen_random_header(index: int) -> str:
        return f">ELAL{random.randint(100000, 999999)}_{index}|test{index}:gen-id|{random.randint(1000, 99999)}"

    def generate_fasta(num_genomes: int, sequence_length: int) -> str:
        fasta_content = ""
        for i in range(1, num_genomes + 1):
            header = gen_random_header(i)
            sequence = gen_random_seq(sequence_length)

            formatted_sequence = '\n'.join(
                sequence[j:j + 60] for j in range(0, len(sequence), 60))

            fasta_content += f"{header}\n{formatted_sequence}\n"
        return fasta_content

    def save_fasta_file(filename: str, num_genomes: int,
                        sequence_length: int):
        fasta_content = generate_fasta(num_genomes, sequence_length)
        with open(filename, "w") as file:
            file.write(fasta_content)

    save_fasta_file(filename, num_genomes=num_genomes, sequence_length=150)


def generate_random_fastq(filename: str, num_reads: int):
    """This is a helper function to generate a randon FASTQ file"""

    def gen_random_seq(length: int) -> str:
        return ''.join(random.choices("ACGT", k=length))

    def gen_random_quality(length: int) -> str:
        return ''.join(
            random.choices("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI",
                           k=length))

    def gen_random_fastq_header(index: int) -> str:
        return f"@ELAL{random.randint(100000, 999999)}_{index}|starfish:taxid|{random.randint(1000, 99999)}_read_{index}"

    def generate_fastq(num_reads: int, read_length: int) -> str:
        fastq_content = ""
        for i in range(1, num_reads + 1):
            header = gen_random_fastq_header(i)
            sequence = gen_random_seq(read_length)
            quality = gen_random_quality(read_length)

            fastq_content += f"{header}\n{sequence}\n+\n{quality}\n"
        return fastq_content

    def save_fastq_file(filename: str, num_reads: int, read_length: int):
        fastq_content = generate_fastq(num_reads, read_length)
        with open(filename, "w") as file:
            file.write(fastq_content)

    save_fastq_file(filename, num_reads=num_reads, read_length=70)


def test_load_fasta_and_fastq():
    """This func test the loading of a FASTA/Q file"""
    ref_path = "random_ref.fa"
    reads_path = "random_reads.fq"
    generate_random_fasta(ref_path, 7)
    generate_random_fastq(reads_path, 25)
    refs_dict = generate_ref_dict(ref_path)
    assert len(refs_dict) == 7
    reads_dict = generate_reads_dict(reads_path)
    assert len(reads_dict) == 25
    os.remove(ref_path)
    os.remove(reads_path)


def test_execute_dumpref(monkeypatch, capsys):
    """This func tests execute_dumpref(), using monkeypatch and capsys to mock
    arguments passed from command line and to capture the printed output to screen.
    func generates a KDB file and send it to dumpref"""
    align = Alignment()
    k = 4
    ref_file_path = "reference12_file.kdb"
    read = Read("@MGYG087_66|kraken:taxid|6005_read_103", "ATGCCGGGGCTAA",
                "AT!!ABCDI?III")
    align.set_reads({"@MGYG000004787_66|kraken:taxid|6005_read_103": read})
    align.build_reads_kmers_db(k, None)
    ref1 = Reference(">MGYG014_62|snail:taxid|3633", "ATGCCGGGG")
    ref2 = Reference(">MGYG061_62|bird:taxid|2211", "GCCGGGGCTAA")
    ref3 = Reference(">MGYG061_41|dog:taxid|2213", "CCGG")
    ref4 = Reference(">MGYG061_11|cat:taxid|2214", "GGGCTAA")
    ref5 = Reference(">MGYG061_34|owl:taxid|4455", "GCTAA")
    align.set_references([">MGYG014_62|snail:taxid|3633",
                          ">MGYG061_62|bird:taxid|2211",
                          ">MGYG061_41|dog:taxid|2213",
                          ">MGYG061_11|cat:taxid|2214",
                          ">MGYG061_34|owl:taxid|4455"])
    genome_bases = {">MGYG014_62|snail:taxid|3633": 9,
                    ">MGYG061_62|bird:taxid|2211": 11,
                    ">MGYG061_41|dog:taxid|2213": 4,
                    ">MGYG061_11|cat:taxid|2214": 7,
                    ">MGYG061_34|owl:taxid|4455": 5
                    }
    ref_dict = {ref1.get_header(): ref1, ref2.get_header(): ref2,
                ref3.get_header(): ref3, ref4.get_header(): ref4,
                ref5.get_header(): ref5}
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    align.set_ref_kmers_db(ref_kmers_db)
    original_kmer_db: KmersDB = ref_kmers_db
    dump_ref_kmers_db(original_kmer_db, genome_bases, ref_file_path)
    mock_argv = ["main.py", "-t", "dumpref", '-r', ref_file_path]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    args = readargs()
    execute_command(args)
    genome_summary: Dict[str, Dict[str, int]] = generate_gen_sum(
        ref_kmers_db, genome_bases)
    ref_DB_dump = {"Kmers": ref_kmers_db.get_ref_db(),
                   "Summary": genome_summary}
    json_output = json.loads(capsys.readouterr().out)
    expected_json = json.dumps(ref_DB_dump)
    assert json_output == json.loads(expected_json)
    os.remove(ref_file_path)


def test_dumpref_with_similarity(monkeypatch, capsys):
    """This function creates a KDB file using execute_reference() with similarity
    flag, and then execute_dumpref() - the similarity output should be printed
    from the previous command"""
    ref_content = """>MGYG000004169_19|kraken:taxid|5353
    TGGTGGGAGACCGTGCAGGTTCAAGTCCTGTTAACCGCAGTAAATTGTGAAAAAGTGGCG
    GTATCGTGAATGATACCGTCACTTTTTTCGTATACCTTTTTGTCCGGTTCTTGGTTAAGT
    TGGATTGGCAGAGGGCTTGTATATAGAAGCTTTGCAACGATACGTTGCGCAGATTCCTAC
    ATGGGCTGGTAGATGTGAAATAGTTTAGGCGTTATGATGATTCAAAATTAATTGCAGGAG
    GAGTGGTCGAATGGGCAAAGGGTTAATTCGGTTCGTAATCATTTTTAGGAGTATGGGATG
    CTTCTGAAGTGCTAAGGAGTCCTTCTTTAAATGTTTTTACCCAGCAGTAAATCGTGCCAT
    TCAGGACATTCTCCTTTCTGCTCCGATGGCTTGACTGGTATCAAGGATGCAATCTCTACA
    GCATTTCCAAAAACGGAGCAACAGCGTTGTATTGTACATATGGTGAGAAACACGCTTAAA
    TACGTTGCAAACAAGGACATGAAAGCATTTGCAAAAGATTTAAAGACAATCTATACCGCT
    GCAGATGAAGAAGCCGCCAGAAAGCAGTTGGAGTCTGTAACAGAAAAATGGTCTGCCCAG
    TATCCAAGTGCGATGAATCGCTGAATTCCTGCTATCGTAGATTGAACAAGCAACGCAGTG
    TATTCCAAAGCTCTCAGGCGCTTATGAAAGCCCTATATCTGGGAACCTTCGAGATTGCAA
    >ELAL0000000_131|filter_me_please:taxid|3322
    TGGTGGGAGACCGTGCAGGTTCAAGTCCTGTTAACCGCAGTAAATTGTGAAAAAGTGGCG
    GTATCGTGAATGATACCGTCACTTTTTTCGTATACCTTTTTGTCCGGTTCTTGGTTAAGT
    TGGATTGGCAGAGGGCTTGTATATAGAAGCTTTGCAACGATACGTTGCGCAGATTCCTAC
    ATGGGCTGGTAGATGTGAAATAGTTTAGGCGTTATGATGATTCAAAATTAATTGCAGGAG
    GAGTGGTCGAATGGGCAAAGGGTTAATTCGGTTCGTAATCATTTTTAGGAGTATGGGATG
    CTTCTGAAGTGCTAAGGAGTCCTTCTTTAAATGTTTTTACCCAGCAGTAAATCGTGCCAT
    TCAGGACATTCTCCTTTCTGCTCCGATGGCTTGACTGGTATCAAGGATGCAATCTCTACA
    TACGTTGCAAACAAGGACATGAAAGCATTTGCAAAAGATTTAAAGACAATCTATACCGCT
    GCAGATGAAGAAGCCGCCAGAAAGCAGTTGGAGTCTGTAACAGAAAAATGGTCTGCCCAG
    TATCCAAGTGCGATGAATCGCTGAATTCCTGCTATCGTAGATTGAACAAGCAACGCAGTG
    TATTCCAAAGCTCTCAGGCGCTTATGAAAGCCCTATATCTGGGAACCTTCGAGATTGCAA
    >MGYG000004832_79|kraken:taxid|6050
    CCGCCCGGGTGAGCCGGTCCCGGGGGTCCAGGGCCGTCTGGGCCACCAGGGCCTCGTAGG
    CCCCGCAGATGAGCCCTAAGTCGGTGAGCTTGTCCCCGCTGGGGCCCTCGGCCTCCTCCC
    CCGCCCGGATGAGCTGCTCCGGGCTTACCCGGCAGCTTTTCAGCTCGTCCACCGTGGCCA
    GCAGGGATTGGAGGAAGGCCGGACGCTTGGAGGGGCGGCCATAGACCTTCAGCTGCTGGG
    CGGTCTCCTGCACCGCCCGATACATCAGCAGCAGCCGTCCGCCCCCGTCCAGCTCCTCCT
    GGCCCAGTCCGCCCGCGGCCTGGAACACCCGGTTGGCCAGGCGGCTGAAGGACAGCACCT
    CTGCCCGCAGGGAGATCCCTGGGCCGCCAGCCTTACATAAGGCCCGCTCGGCCTCGTGGG
    ACTGCTGCTCCGGGGTCATCAGCACCTGGGGCCGCTCCCCGGCGGCCTGGCACAGGCGGT
    TTAGCACGGCGGTGGTCTTGCCGCTCCCGGCCCGTCCCATGAGAATGCGCAGCATCCGGC
    AGCCCTGCCTCCACGTCTCCCCCGGCGTTGTTCAGCAGGATGAGCAGGCCGTGGATGTCC
    GGGTCGCCCTCCACCGCCGCCAGCAAGGGCAGCACGTGCTCGTATTTGGTGCTCTTCATG
    CCGGAGGGCAGCAGCTGGTGCCCCTCGATCTGCCCCACGATGGTCAGCACATAGACCCGC
    CCCTCCTCCCGCTGGATGAGGGAGGTACCCAGGTCCACCAGGTTCTGCGTGTTGGCCGGG
    GCCTCCGGCTGTGTGTGCTCCTGTTCCATGTCTGCCATTTCTTGCGCCTCCTAAACCGGA
    ATTTCTGTTAGGGTGCGCGGAATAGGGCCGGTTTACTCACAAAAATGGGCGCGCGGCGTG
    ACGAAAATTCGTCACGCCGCTGGTGTTCCGGCTAAAATCAGCCCCGCTCACGCATCTTGG
    CGAAGCAGTCGCTGCAATAGACAGGACGGTCAGACTTGGGCTCGAAGGGCACCTTGGCCT
    CGCCGCCGCAGGCAGCGCAGACAGCGGTGAAGTACTCACGGGGGCCACGAGCCGCGTTCT
    TGCGGGCATCACGGCAGGCCTTGCAGCGCTGGGGCTCGTTCTGGAAGCCGCGCTCGGC
    >MGYG000002096_86|kraken:taxid|3024
    GGGGCATAGACATAGCCATACTGCGCTTTTCCGCTCAAAGAATATCCTTTGATAAAACGC
    CCTTCCCAGCGGCCGTCTTTCCGTTTGCGAATATTTTCTCCTCTTCTTGGCATGCGAACT
    CCTCCTTACAAGTATATGTTGTGTTTTGACCCCTGCATTTGACCACGAATAGACCCCCGC
    ATGCAAACAGTCCCGTATATTGCCGATATCATCGGGAGCAGAATTGACAAATATAGCGGG
    CAAGGCTTATAATTAGCATAGCTTTTGCATAGGGCAACGCATAAATAAAAGAGGTTTCAG
    CAGGCAAATATCCTACCTATAGCAATATACTTTCGGAAAGAAAATTGAGCGGCTGGGTAA
    TTTTGTTCTGTGGCCAAAAAAGGGAGAGGATAAAATCAGCGGTCGCTGCAAAAAGAAAAA
    TGGCGCAGTGGTATTATACCATGCGGAAACAAACGGATTTTCAGGGCTCATAGTTTGCAG
    AAAATGGGCAACAATATAGGAAAGGATGAACACTGTTTGAAGAGAGGAGCACGCGGATGG
    GGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCA
    TCGGCGGCGACGCTTCATATTCTGGGAACGGCCCCGGGCAGCAATGATCTTCAGGAAAAT
    TGTGAAATTTTATTGATTTCAGATAACAAAACCGGCGCGGTTTTGACGGAAAAGAACGCG
    GGCGAAAAGATGACTCTCATCGGCGGCACGGTCCGCTGGATGGTGGCTCTTACCGCCGTC
    GATCACCTGGAATTATCCGAGGAGGTGACTCTGGAGGCGGCAGATCTCGAGCCCTTTCAG
    GGAAACCGCAAAATCGGCCTTCGGGCGGGACAGACTCTGACGGTCAGGGACCTCCTGGCG
    GCCATGATGGTGGACTGCGCCCAGGATGCCGCCGTGGCGCTGGCAAACAAGGCGGCGCAG
    AAGGCGGGAGCGGATGATTTTGTCGCTCTGATGAATGAGAAGGCAGCCGAACTGGGAATG
    AAGGATACAACATTTAAAAATGCCACCGGCGCGCTGGACAGCGGGCAGGTGACAACGGCG
    AAATTGTTCCCTCCTATCAATGGACGGGGGAGAGCGATACCTCCCAGAAACAAGAGGCAT
    AAATGCGGCAGGCTGGACCAATTGCCGGCTGCAAATGTGGCTAAAATGAAGAAAGCAGCA
    ATATTTATGTGTTCTCATACAGAAAAAGAGAACCTGAATCAACGGCATTTGGAAGCAATA
    TGCAATGGAGTAAAGCAGTGTGAATCTTATGTCCGCCTGTGTAGCCGCAGCCTTAAAAAG
    GCCAAGTGGATTTAACCGCGGGAAGAACTGCACAAGCAGTGAGCATAGTTGGGGCCGGAT
    TATGTTGGGCTATGCTGGGCGGCTAATGAACGGCGGGGGAAAATGAAAAATTTTCTATAT
    AGCAGCATCTATTCGATGTAAAGAGGCGGTTTCTCAATTTCCCAATTAGGTATGTAATGA
    TCTGTGGGCCAATGGCATGAAATAGTTATAAAGCGAAAATAAACTTGAGGGAAAGCAATT
    CTTATGCCTCATCAGGGAAATAATTGATCTGGCGCAGTGGGGGAAAATAGCTGCAGCCAA
    AGCTGGATAAATTATATAAATAAATGTACTCCGGGATGTTGGGATA
    >MGYG000002720_43|kraken:taxid|3759
    GTTCTCCGGCCGAGCTGCTGAAAAAGAACGGCCTGTTCCGCCATATGGCGCAGCTTCAGA
    GCGAGAGCATGGAGTGGACGGCGTAGAAAGTGCTGAGATAAACAGGATTCCAAGTGCAAT
    ACTGGAAAACGAAGTGTCGAAAACAAGGAGGCGTTTATGTAAATAAATGGTCTGTCAGGG
    ATGTGATCACCACGGTATTGCTCTCCGCAGTACTTATCGTGATCCAGCTCGTTGTCAATA
    TGGTCTGCATGGCCAACGACTTTGTCAGCATGGTGCTGTCGGTGGGGATCACCATGTTCC
    TCTGCGCGCCGGTCTATATGCTCATGGTCAGCCGGATCGGCAAGCGGTTTGTGACGCTGA
    TCTATATGACGCTGCTTGGCGTGATTTTCCTGCTGATGGGAAACTGGTTTTTGCTTCCCT
    ACTTTATCGTTACAGGCGTCATCTGCGAGGCGATTCTCTGGAAGGAAGGCTCCTGCCAAA
    AGCCGAAACGGCTGACCGCCGCCTGGACGGTGGCGAGCCTGCTGTACAACGGGGTCAATC
    TGCTCCCCATCTGGTTCTTCTGGGACACCTACTACGATTTTGCACTGGCAAGCGGCATGG
    AGCAGAGCTATATCGATTCCTATGTGCGTTACTACACCTCTCCCGGCTGGCTGGCCTTTA
    TTCTGCTGTTTACGACGCTGATGGGCTTTTTAGGCTGCATGGTGGGCAGTCGGCTGATCC
    GCAGGCATTTTAAGAAGGCCGGCGTTCTATGAGGGCGGACGCATCCTTTGCGGTGCCGGT
    CAAGCTGTGGGCGCTGCTGTGCGTCTTTGCCGGAGTAACCATCGGCGGGAATGTGCTGCT
    CACCTGCATCCTGACCGGCGGGGCGCTTCTTTATCTCGTCCTGCAGCGGAACTTCCGCCT
    TGCCGCGTCCTATGGCTGCTTTTATCTGCTGCTGGCGCTGCTGCTGTACGGCATCCGCTT
    CCACGGGCTGCACATGCCGGTCTTTTCGGAATTTTATGTTCTGATGTTCTGGAATCTGTC
    ACCGATCTTTCTGGTGTCATGGGATCTGATCACCACGCCGCCGGGGATGCTCTCCGCGTT
    TTTATCCCGCCTCCGTATGCCCACCCCGTTCATTCTCGGGCTTCTTGTGGTCTTTCGCTT
    CTTCCCCACCATGCGGGCGGAGCTAAAGGGTGTAGGCAGGTCCATGAAAAACCGCGGGCT
    GACCGCGGCGGGGCAGCTGCTCACGCATCCCGTGCAGAGCATGGAGTATGTGCTGGTTCC
    GTTTCTCCTTCGGGTGTTGCAGCTTGCCGATCAGCTTTCCGTTTCCGCTGTTGCAAGAGG
    CGCGGAGCGTCCGGGCGTGCGCGGCAGCTATTATGAGAAAAGGGGCGGGGCGCGGGATCA
    CATCGCGGCAGCCGCATGTGCGCTGGTCACGGCTTCCTATTTGGTTTTGGAAAGGAGCAT
    GGCATGATCGAGCTTACACGCGCTTCTTTTCAATACGAAAACTCAGACAGAGGCGTGCGG
    GACATTTCCCTTTCGGTAAAAAGCGGCGAATGCGTCGTGCTGACGGGGCTCTCCGGCTGC
    GGCAAGACCACCGTGACAAGGCTTGTCAACGGGCTTGCGCCATTCTACTATCCCGGAGTT
    TTCAGCGGGAGCGTGAGGATAGACGGCAAGGACATATCAAGGCTTTCCACATGGGAGATC
    GGGCGTTTGGTGGGAAGCGTTTTTCAGGACCCCAAAAGCCAGTTCTTTTCCTCCGAGCTT
    CATAACACGATGCCGCCCCGATGTCAATATGCAAATTCGTCGCAAATTAGCCGCCAACGC
    GGATAGTGCCCGCTCCGCTGAAATGGGCTGCGGCAGAGCTGAATGATCTGCGCTCCGTTT
    TCGGCGGCAAAAAGGCAGCGGCAAATGCCAGGTGAAAAAGCGGTAGAAAAACCTGAAAAG
    ATATGCTTCACACTTATCAAACCCCAAAAGCTCTTTCAGAGCAAAAGGATTTCATTACGA
    ACAGCAGCTGACCCTCCGTATGTACTGAAGCATACGGAGGGTCAATCGGCTTGTGTGGCC
    TGCTCAGTCGAGCAGGGGCTTTCTTTTCTGATGCGTCAACATATCCGTGCCCGGCCTATC
    CCTCAAACAAATTCAGGAGGTATCCGTAGCCCTCGGCCTCCATTTTCGCGTAAGGCACGA
    AGCGGAGGGCGGCGCTGTTGATGCAGTAGCGCACACCGTTGGGCGACTCCGGGTCGCCGG
    TGAACACATGACCCAAATGGGAGTCACCGGCGCGGCTGCGCACCTCTGTGCGGCGCATGC
    CGTGGCTCAAGTCCTCCCGTTCCACCACAGCAGGACTCTCGATGGGCTTTGTGAAGGCAG
    GCCAGCCGCAGCCGCTTTCAAACTTATCCGTGGAGGAGAACAGCGGCTCGCCGGTGACGA
    TGTCCACATAGATGCCCTTTTCAAACTTGTCCCAGAGCTCACCGGTAAAGGCGCGCTCCG
    TGCCGCTCTCCTGCGTGATGCGGTACTGCTCCGCCGTCAGCGTGTCCCGGATGGTCTCCG
    CCGCGGGCTTCCGGTAGTCGCCCGGATCGATGCGAAGCCTGGAGAACAGCTCCATCTCCG
    CCCGGGGGATGTGGCAGTAGCCGTTGGGGTTCTTTTCTAGGTAGTTTTGGTGATATTCCT
    CGGCGGGATAGTAATTCTTCAGCGGACCGATCTCCACAAAGAATTTCTCGCTGCGGCCCC
    GCTCAATGCCTGCGATGCGCTCCACCGTTTCCCTGACCTTGTCGTTGGTGTAGTACACGC
    CGGTCTGATACTGGCTGCCCCGGTCGTTGCCCTGCCGGTTCTGCACCGTGGGGTCGATCA
    CATAGAAGTAAGCCAGCAGCAGCGCGTCGAGACTCACCTGCTCGGGGTCGTATTCCACCC
    GCACGGTCTCCCGGAAGCCGGTGTCGCCCTTGCAGACGGTCTGGTAGGTGGCGTCTGCCT
    CGCAGGTGCCGTTGGCATAGCCGCTTTGGGCGTCGATGACGCCGGGGATGGATTGCATGA
    GCTGCTCCATGCCCCAGAAGCAGCCGCCTGCCAGATAGATGACATTTTCCGTGGTGTCCA
    TATATATACCCTCCTCGCCGGACATATCCGAGGACTGTCCGGTTTTCTCCGCTTCCTGCG
    TGTCCTTTTCCCGCGGTACACGCTGCCGTTCTGCGGCGCAGCCGCTGAGCAGCAGTACCG
    AACTTCGGCAGGCTGCTGTCATCCGTCTTGCTGCCCGTGCTGCCGCTATCGGTCTTATTG
    CCTGTGTTTCCGCTGTCGGTTTTATTGTCCGTACTGCCGCTGTCGGTTTTGTTGTCCTTA
    TTGCCGCCGTCGGGCATCGTGTCCATATCCTTGTCCGTCGGCATATTCATATCGCTGTCC
    AGAGATTTTTGCATGATCTCGGGATACTTTTCCTCCAGTGCGGCCAGCTTGTTTCCGATG
    TCCACGCCGTGGATATATGCCTTTACCAGCTCCCGCGCCGCCTGCCGCATCCGGCCCGCA
    GATGAACGGCAGCCAGAGCAACCACATGGCGATGGGAGGGCTTTCACGATGATCATTCAA
    GCAGAACTGAAGTGTAAGCAGACCGGGTGCGAGGCAGGCCCCTGTGCCGTGGATAAGGTC
    ATCGAGCTGCCAAGCCCGCGGTTTCAGCAGTTCAGCCGCGCACTGCTGGCTGATTACGAT
    TTCATCGCAGAGAATAAAAACGCCATCCGGCACGATGAGGATGCCAGGCACTGTCTGCTC
    ATCCTCGACGCGGATGGAACGGACGGTTTCCTTGTTGACCCACAGGGGCACAACTACGCC
    CGGTACAGTGCCTTTGTCCCCAACGCCCGCAGTTTGCTGACGCCGGATATGGGGATCGAC
    CGCAGCTATCTTTCGCCGGCAGAGCCTT
    >MGYG000000975_36|kraken:taxid|1603
    TGACGTTGTATTCTTTGTTTGTACATCTATTTAGTTTTCAAAGATCTCTTGCACTAATTG
    ATTCTATCAATCTTAACGAAGAATTAATCCTTCGAACAAGTGCTTAATAATTATATCTAA
    GAGATATAATCATGTCAACAAAGAAGTCTCAAATCTTGCGACAACTTCAAAATGATACAC
    CTGTTCGATATTGTAAATCAAGCACAAACGAAATAGTATCTCTATCCTTATTATATATAT
    AAGAGTTTTTTATTCCTGAAAAGTCCTCCAAAAAAGACAAATATTCATTTTCAATTTCTC
    ATACATAATCCATGTGTCAACCGTTCCAAACATAAGTTCGCCTTTTTCAGCTCTTTCTTG
    TCCATCTTCAATATTATCTAAAATAAATCTTATCTTTGAAGCACTGAAATAAGGGTTAAT
    TATAAGTCCAGGCAAATTGGAGCAGATTGCCTGGATTGCCAAACAATCGCATTATAAACA
    GGCATACCAGTTTCTTTTTCCCAAATAATTGTTGTTTCTCTTTGGTTTGTGATACCAACA
    GAATCGATGTCATCATATGTTAAATTAGTTTTTAATAAAATATCATTGATAACATAAATT
    ACACTTAAATATAAATCCAAAGCGTTTTGTTCAACCCAACCACTTCGTGGAAAAAGACAT
    TCTACTTCTTGTTGAGCTTTGAAGAAAACTTTTCCTTCATGATCAAAAAGCATCGCACGA
    GTTGATGTCGTGCCTTGGTCAATTGTTAAGATGTACTTAATAGTTAATTAAAACTTTAAA
    TCTGAATTTAAAAACGAATAATAAGCTAAATTCCCACCTTTAGAATGTCCTAAAACTATA
    CTCTACGAGGATCTTTGAGCTCTTGCTTAACATATTCATGTCCTCTTTTACATGTGTAAC
    CAGTAATTTCATCATCTTTAACATGAATATGGCAGCCTTTAGGGCAAACAATACAAATAA
    AATCTTTTTCCATATTAAAGTTCCTCAATTTCTAAACTTAAAGAATTCAAAGATTTATAA
    TCAACATCTTTTACTGAAATCTTTTGCATTTCTGATGGAATTAAATCTACCTTCATAATT
    TGCATCACTTAAATCATCAATATTATTAAGCAAAACAAAAATGTCTTCTCTATCAATAGG
    AGTAATAAATTCTTTTGCTAATTGTTCTTCAACATCATCAGCTTGATGTTCAATTGCATG
    AACACTATTTTTTAATTCTAATAATTTAGAATTATCAAAATTCAAAAGTCCTTCGCTTAA
    TACGTCAATTTCGTTTTTAGAATATTCAATCAAATCTTTAAATGATTGATAGAAATAATT
    CTTTTCTTTTTTCTTAAAAAACATATAAAACTCCTTACCTTGCTATAGCAATAAATATAA
    AAAATGGCCTTGTCCTATTTCTATTCATTTTTCTACAAATTAATTGAATAGTTTTAGTAA
    CTAGCCATCCTACAAAAAATCCAAATACACAACTAAGAATTAAACCATAAATTGTAAATG
    CCCAGCTTGAATCTGGACCAAATTTTATAGTTGCATCTAATCCCATTGTCATTGCAGCAA
    GTCCTGCTCCTGTAATTCCAGCAATTAATGCATGTGATTCACTACTTGGAATACCAAACC
    ACCAAGCAAGTACTCCCCAAACGATAACACCAACCATGCCTGCACAAAGTGCACTTAAAG
    CAAGATTGTAATTAGTACCAAAATCAACCATGTTAGAAATAGTTTCAGCAGTTTGCTTAG
    ATATAAAAAACATAGTTAAAACACCTAAGAAATTGAATATTGCCGCCATTATAATAGCTG
    TTTTTGGCTTGATACATCTTGTTGAAATGCAAGTTGCAATTGCGTTAGGTGCATCTGTCC
    AACCGTTAATAAAAATTACACCAATCACTAGGATAGTTGCAACGCTAAGCATTGGATAGG
    AACCAAATAACGCAAATATTTGTTCTAAGCTCATACCTAAAAATTATAAATCAGTGCTTT
    TTAGAATGTAAAGTTTTTCCTTTATATTTGATAGCGTTTCATTATTGTAGGTAACAAAAC
    AAAAGAGAACTCGTATAATTAAGC
    >MGYG000001129_3|kraken:taxid|1782
    CTCTTTAAACGATCGGCAATGGAAATGCGGCGATCAGCTGCCATGAAATAGCCGGCATGG
    GTTTGTAACATTTTGGAAAGGCTTGCCGTTACAGCACCGCCGAAAATAGTAATCAATTCA
    CAGCCAAGTGCAAGCCATGCCGTTTTGCTGCCTTGTTTTCCATCGGTCAGTGCCAGTACA
    ACAAAATAGATTGCTGCCACTTGGAACATTTGCAGGATTGCATTTAAAAATCCAACTGCC
    ACAGACTGATTGATGTTTTTTGTTTCTTTTCCGGCAAACTGCCAGATTTTACGGAATGTT
    TCAATCATTTTGTGTCACCGTCCTTTGCTCCAATATGCGCCTGCCACATTTCATGGTAAA
    GTGCACAATTTTTCAATAGCTCTTCATGCGTTCCCTGCGCTGCAATATGACCATCCTGTA
    GGACAACAATGTTGTTGGCATCCGTGATAGTGGAAAGGCGGTGCGCAATCACAATGACGG
    CGGATGTGCCGCCGTCGGAATATTCCGGCATCTGAGCAGCCTGTGCAGCACGGTCCATTT
    GGTGAAAAGTTTCTTGAATCTTGAGAAAGCGTTCGATCATTTCACTGAGCACAGGCTCAA
    CGATTGCGGAAAAAAGCCCTTCCTTATCTCCAAAGCGGACATAAATGGAGTTGGTACTGG
    TATTGGCAGCTGTTGCAATGGTGCGCAAAGAGGCATCTACATAACCTTTTTCTAAAAATT
    CCTGTTCGGCAGCTGTAAGAATGCGTTCTGCAACACCTTCAATTTGTTTTGCCATTTGTT
    TCGCCTCCTTTCATTCATAACAGCGTTTCATAACACTGTTTTTATTATATACCAGCAAAT
    TTAAAAAGTCAATACCTCAGGATGGATTTCGTTGGAAATCAGTTGAAAATCCGTATATGC
    TATGATATAATGAGAAAAA"""

    ref_path = "test_extsim_dumpref.fa"
    kdb_file_path = "test_extsim_dumpref.kdb"
    with open(ref_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    mock_argv = ["main.py", "-t", "reference", '-g', ref_path, '-k', "31",
                 '-r', kdb_file_path, '--filter-similar']
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    mock_argv = ["main.py", "-t", "dumpref", '-r', kdb_file_path]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)
    if len(outputs) != 2:
        pytest.fail(
            "Expected two separate JSON outputs, but got something else.")

    similarity = {'Similarity': {
        'ELAL0000000_131|filter_me_please:taxid|3322': {'kept': 'no',
                                                        'unique_kmers': 29,
                                                        'total_kmers': 630,
                                                        'genome_length': 660,
                                                        'similar_to': 'MGYG000004169_19|kraken:taxid|5353',
                                                        'similarity_score': 0.953968253968254},
        'MGYG000004169_19|kraken:taxid|5353': {'kept': 'yes',
                                               'unique_kmers': 89,
                                               'total_kmers': 690,
                                               'genome_length': 720,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000001129_3|kraken:taxid|1782': {'kept': 'yes',
                                              'unique_kmers': 889,
                                              'total_kmers': 889,
                                              'genome_length': 919,
                                              'similar_to': 'NA',
                                              'similarity_score': 'NA'},
        'MGYG000004832_79|kraken:taxid|6050': {'kept': 'yes',
                                               'unique_kmers': 1108,
                                               'total_kmers': 1108,
                                               'genome_length': 1138,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000002096_86|kraken:taxid|3024': {'kept': 'yes',
                                               'unique_kmers': 1636,
                                               'total_kmers': 1636,
                                               'genome_length': 1666,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000000975_36|kraken:taxid|1603': {'kept': 'yes',
                                               'unique_kmers': 2034,
                                               'total_kmers': 2034,
                                               'genome_length': 2064,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000002720_43|kraken:taxid|3759': {'kept': 'yes',
                                               'unique_kmers': 3838,
                                               'total_kmers': 3838,
                                               'genome_length': 3868,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'}}}

    summary_output = json.loads(outputs[0])
    similarity_output = json.loads(outputs[1])
    assert similarity_output == similarity
    os.remove(ref_path)
    os.remove(kdb_file_path)


def generate_random_kdb(genomes_path: str, kdb_filename: str, num_genomes: int,
                        k: int):
    """This is a helper function to generate a random KDB file using generate
    FASTA() function"""
    generate_random_fasta(genomes_path, num_genomes)
    ref_dict = generate_ref_dict(genomes_path)
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    dump_ref_kmers_db(ref_kmers_db, genome_bases, kdb_filename)


def test_dump_and_align(monkeypatch, capsys):
    """This func tests execute_dump_and_align(), using monkeypatch and capsys to mock
    arguments passed from command line and to capture the printed output to screen.
    func generates FASTA/Q files, cretes a KDB file and send them to dumpalign
    and then check the expected result was printed"""
    ref_path = "random_ref.fa"
    reads_path = "random_reads.fq"
    kdb_filename = "ref_file123.kdb"
    generate_random_kdb(ref_path, kdb_filename, 100, 31)
    generate_random_fastq(reads_path, 50)
    mock_argv = ["main.py", "-t", "dumpalign", '-r', kdb_filename, '--reads',
                 reads_path]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    tup: Tuple[KmersDB, Dict[str, int]] = load_kdb_file(kdb_filename)
    ref_kmers_db: KmersDB = tup[0]
    genome_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genome_bases.keys()))
    aln.set_ref_kmers_db(ref_kmers_db)
    k_size: int = 31
    aln.build_reads_kmers_db(k_size, None)
    aln.map_reads(1, 1, None, False)
    reads_stats: Dict[str, int] = aln.get_reads_stats()
    gen_map_sum = extract_gen_map_sum(aln.get_reads(), aln.get_references())
    dump_align: Dict[str, Dict[str, Any]] = {"Statistics": reads_stats,
                                             "Summary": gen_map_sum}
    json_output = json.loads(capsys.readouterr().out)
    expected_json = json.dumps(dump_align)
    assert json_output == json.loads(expected_json)
    os.remove(ref_path)
    os.remove(reads_path)
    os.remove(kdb_filename)


def test_execute_dumpalign(monkeypatch, capsys):
    """This func test execute_dumpalign() to check if it prints additional data
    from a previous align mapping process. func test quality and coverage"""
    ref_path = "test_dumpalign.fa"
    reads_path = "test_dumpalign.fq"
    aln_path = "test_align.aln"
    generate_random_fasta(ref_path, 50)
    generate_random_fastq(reads_path, 100)
    mock_argv = ["main.py", "-t", "align", '-g', ref_path, '-k', "31", '-a',
                 aln_path, '--reads', reads_path, '--min-read-quality', "21",
                 '--min-kmer-quality', "33"]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    mock_argv = ["main.py", "-t", "dumpalign", '-a', aln_path]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    dumpalign_output: str = capsys.readouterr().out.strip()
    mock_argv = ["main.py", "-t", "dumpalign", '-g', ref_path, '-k', "31",
                 '--reads', reads_path, '--min-read-quality', "21",
                 '--min-kmer-quality', "33"]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    dump_and_align_output: str = capsys.readouterr().out.strip()
    assert dumpalign_output == dump_and_align_output

    mock_argv = ["main.py", "-t", "align", '-g', ref_path, '-k', "31", '-a',
                 aln_path, '--reads', reads_path, '--coverage',
                 '--full-coverage']
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    mock_argv = ["main.py", "-t", "dumpalign", '-a', aln_path]
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    dumpalign_output_cov: str = capsys.readouterr().out.strip()
    mock_argv = ["main.py", "-t", "dumpalign", '-g', ref_path, '-k', "31",
                 '--reads', reads_path, '--coverage',
                 '--full-coverage']
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    dump_and_align_output_cov: str = capsys.readouterr().out.strip()
    assert dumpalign_output_cov == dump_and_align_output_cov
    os.remove(reads_path)
    os.remove(ref_path)
    os.remove(aln_path)


def test_build_ref_align_and_dump_reverse(monkeypatch, capsys):
    """This func tests execute_dump_and_align(), using monkeypatch and capsys to mock
    arguments passed from command line and to capture the printed output to screen.
    func generates FASTA/Q files and send them to dumpalign and then check the
    expected result was printed.
    Func is testing function with --reverse-compliment flag"""
    ref_path = "random_ref.fa"
    reads_path = "random_reads.fq"
    k_size = 31
    unique_thresh: int = 1
    ambig_thresh: int = 1
    mrq: float = None
    mkq: int = None
    mg: int = None
    generate_random_fasta(ref_path, 100)
    generate_random_fastq(reads_path, 200)
    mock_argv = ["main.py", "-t", "dumpalign", '-g', ref_path, '-k', "31",
                 '--reads', reads_path, '--reverse-complement']
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    refs_dict: Dict[str, Reference] = generate_ref_dict(ref_path)
    tup = build_ref_kmers_db(refs_dict, k_size)
    refs_kmers_db: KmersDB = tup[0]
    genomes_bases: Dict[str, int] = tup[1]
    reads: Dict[str, Read] = generate_reads_dict(reads_path)
    aln: Alignment = Alignment()
    aln.set_reads(reads)
    aln.set_references(list(genomes_bases.keys()))
    aln.set_ref_kmers_db(refs_kmers_db)
    k_size: int = 31
    aln.build_reads_kmers_db(k_size, None)
    aln.map_reads(1, 1, None, True)
    reads_stats: Dict[str, int] = aln.get_reads_stats()
    gen_map_sum = extract_gen_map_sum_rev(aln.get_reads(),
                                          aln.get_references())
    dump_align: Dict[str, Dict[str, Any]] = {"Statistics": reads_stats,
                                             "Summary": gen_map_sum}
    json_output = json.loads(capsys.readouterr().out)
    expected_json = json.dumps(dump_align)
    assert json_output == json.loads(expected_json)
    os.remove(ref_path)
    os.remove(reads_path)


def test_reverse_data():
    """Func test reverse_data function of Kmer class"""
    kmer = Kmer("ATCGCCAT", [1])
    assert kmer.reverse_data() == "ATGGCGAT"


def test_calc_mean_quality():
    """Func test clac_mean_quality() to make sure calculation is done properly"""
    quality = "ABCD22@@!!PBVC"
    score = calc_mean_quality(quality)
    expected_score = 397 / 14
    assert score == expected_score


def test_print_kmers_db():
    """Func is testing the printing of a KmerDB object"""

    def print_db(db: KmersDB):
        """This func return a printable string of the DB"""
        result = ""
        DB = db.get_ref_db()
        for kmer_key, locs_in_refs_dict in DB.items():
            result += f"K-mer: {kmer_key} \n"
            for ref_header, lst_of_locs in locs_in_refs_dict.items():
                result += (
                        f"appear in reference header: {ref_header}" +
                        f" in locations: {lst_of_locs} \n")
        return result

    db = KmersDB()
    assert str(db) == print_db(db)


def test_determine_status():
    """This func is testing determine status function in class Alignment to make
    sure the pseudo alignment algorithm is done properly, the func is mocking
    a result from the function __validate_mapping_results() from class Alignment
     to create a situation in which a read's results were Unique and changed
     to ambiguous and all sources of the read are extracted. Func also checks
     class alignment updates the counters it has for each mapping result"""
    k_size = 3

    # Mock reference and read content
    genome_content = """>MGYG000004832_79|kraken:taxid|6050
AAACTGTTGCC
>MGYG000004832_78|kraken:taxid|1111
ACCCTGAAAGG
    """
    reads_content = """@MYREADS0000|read_1
    AAA
    +
    !!!
    @MYREADS0000|read_2
    CCC
    +
    ???
    @MYREADS0000|read_3
    AACCT
    +
    ?!?!!"""
    genome_path = "test_genome.fa"
    reads_path = "test_reads.fq"
    with open(genome_path, "w") as file:
        file.write(genome_content)
    with open(reads_path, "w") as file:
        file.write(reads_content)
    refs_dict = generate_ref_dict(genome_path)
    ref_kmers_db, genome_bases = build_ref_kmers_db(refs_dict, k_size)
    reads_dict = generate_reads_dict(reads_path)
    aln = Alignment()
    aln.set_reads(reads_dict)
    aln.set_references(list(genome_bases.keys()))
    aln.set_ref_kmers_db(ref_kmers_db)
    aln.build_reads_kmers_db(k_size, None)
    read = Read("MYREADS0000|read_3", "AACCT", "?!?!!")
    all_kmers = aln.get_reads_kmers_db()["MYREADS0000|read_3"]
    tup = aln._Alignment__sort_specific_kmers(all_kmers, None)

    def mock_validate_unique_status(map_count, combined_counts, p):
        return "Ambiguous"

    aln._Alignment__validate_unique_status = mock_validate_unique_status

    result_status, read_sources = aln._Alignment__determine_status(tup[0],
                                                                   tup[1],
                                                                   m=1, p=1)
    assert result_status == "Ambiguous"
    assert read_sources != [], "The read should have ambiguous sources added."
    assert aln._Alignment__total_ambiguous > 0
    os.remove(genome_path)
    os.remove(reads_path)


def test_add_filtered_stats():
    """This func test the add_filtered_stats() function used for quality extension
    to create a summary printed in dumpalign"""
    F_Q_READS, F_Q_KMERS, F_HR_KMERS = "filtered_quality_reads", "filtered_quality_kmers", "filtered_hr_kmers"
    reads_stats = {F_Q_READS: 0, F_Q_KMERS: 0, F_HR_KMERS: 0}
    add_filtered_stats(reads_stats, 1, 1, 1, 1, 2, 3)
    assert reads_stats == {F_Q_READS: 1, F_Q_KMERS: 2, F_HR_KMERS: 3}


def test_execute_align():
    """This func test the entire execute align function without creating an Alignment
    object in the test function. Arguments are supplied to execute align and
    FASTA/Q files are created. Func make sure the mapping results are accurate"""
    # Define thresholds for the test
    unique_thresh = 1
    ambig_thresh = 1
    mrq = None
    mkq = None
    mg = None
    rev_comp = False

    ref_content = """>MGYG000004169_21|test1:taxid|51124
    GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCAG
    CTCCGGCTGACAGCAAAACCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
    >MGYG000004169_22|test2:taxid|51315
    GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
    GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
    >MGYG000004129_23|test3:taxid|51316
    GCCCATTGTCAAGAGCATTACGGCGTGTCCATGACGGTACCAAGTGGCAGACGCGGTAGCTTGGCAAGT
    GCATTACGGCGGCGGCC
    GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGAAAAAA
    ACTGGACACGGCCCACGTGCTCGTATTTGGTGCTCTTCATGCCGGAG
    >MGYG000004179_33|test4:taxid|51316
    GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
    GCTCCGGCTGACAGC
    >MGYG000004119_36|test5:taxid|51316
    GAGACGGAGAGGGTTATCAGGTCTTGATGGAACTACTGAGTACTGATCCGGATACCATTGG
    AAACTTTATGTCCGCGTGCTCGTATTTGGTGCTCTTCATGCCGGAGGGACGGCGCTCATTG
    CGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGGAGGGACTGGAAATTATCATTGCCCATC
    CCGAGCGATATCAATATGTGCCCTGGGGATCGCTCTGGACTGGCAGG
    """

    reads_content = """@MGYG000002096_86|kraken:taxid|3024_read_0
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCAGCTCCGGCTGACAGC
+
DC@@B>>A>BBA??B=@CB@==?A><@??=B>><==A??>=@><;;>?>=@;?>?99>=>>:;>:8>8;87;97:
@MGYG000004832_79|kraken:taxid|6050_read_1
GCCCATTGTCATGAGCATTATGGCGGCGGCCTATCCAGGTTTTGAACTGGACACGTGCCAGCTCCGGCTGACAGT
+
B@AC??BC>>?CA@B@?A@?==BB@<@A?=@>=A?<AAA?@A;@;;><?=<=;=;9;>:<?<98;<==987887<
@MGYG000002096_86|kraken:taxid|3024_read_2
ACCCGCATCCGGAATGCGGACGATGTGGAGAAGTTACTGGATCTGCCGGTTTTGGGGCAGATTCCGCGGGAAGAT
+
B@?CCD>C>BD>ACCB?A>@B>?A?@=BABB@@A>==;=>;<@==>==><<:??;9?;?=:9=<9<:<>;<;<87
@MGYG000000975_36|kraken:taxid|1603_read_3
CAGCACTATGACCATAAACGCTTCCATAAACAATTAATGCACCTTCAACTTCAGGTGTATATGAACTCCATTTAT
+
DDE@DB@@AA>BA?=?=B=>CA?CAB=??B<B@?><;;@;@@;;??@<;><?>?=>::<>?8<=;<=;9879:8;
@MGYG000001129_3|kraken:taxid|1782_read_4
CGCTGGAAGTGTCCATGACGGTACCAAGTGGCAGACGCGGTAGCTTGGCAAGTAATTTTTCTCGAATGGTTTTCA
+
@?BBC@A@>D@@BA@@=CBA@?CB?>>?=?<BA?A=@=>==><=;=<@;>?;:>>>9>?>9>=;8<;89<9::8:
@MGYG000004787_66|kraken:taxid|6005_read_5
GGGTCGCCCTCCACCGCCGCCAGCAAGGGCAGCACGTGCTCGTATTTGGTGCTCTTCATGCCGGAGGGCAGCAGC
+
@DDEEDDD@@A>@B@>>AB?C>=@=?@<A<?<<?;A@?<A?A?=@<?=:;<;=;9>>?9<=>8>;<:8>>=9;9<
@MGYG000004787_66|kraken:taxid|6005_read_6
ACGAAAATTCGTCACGCCGCTGGTGTTCCGGCTAAAATCAGCCCCGCTCACGCATCTTGGCGAAGCAGTCGTTCA
+
BD?CCBCD>>>ADB=CC@B>@>BAB>A><>=A@?<>==A<??>?<=@;;@=>9==>><9<98<:<>;88>=9::<
@MGYG000000125_5|kraken:taxid|343_read_7
AGACTCAGATTATTTCTGTACAGTTTGGACTACATCTGCGGCATCTTTTATATTATCATAGGCAGGAGCCAACAA
+
?DABEA@BADB?A>@==>CBCBA?@<=@B=>?>?A=?<;A<=@;><<;=:<?<=?9;?>:9:9:=<8=><;9:98
@MGYG000004169_19|kraken:taxid|5353_read_8
TCAGTACCTTGGCAGAAACTTTGACGAAGGTATCCGGAACAATCATGAAGTCTGACGGAGCAGCGGAGTTAGGCA
+
BBEBDC>D?>B?@D=B@=?A@CAA@@@AA>A=><?A<><;=?@;>@;@;<:9<<==;:;<9>98;9<;;8=;<77
@MGYG000002096_86|kraken:taxid|3024_read_9
GCGGATGGGGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGG
+
?EC@@?B?CB@D@CCBBCB?=?CAB?<>B?A>A@;@>A<??=?;<@;@<@;>?9:;<:<=98=9=;8;8>;9;=;
@MGYG000004169_19|kraken:taxid|5353_read_10
ATGAATTTTTGAGAAAAGCAGGTGAACAGGTTGACTGAAAATTGTAAGAAATATGCAGTAGAGCTGCATGAACTC
+
A@?CE>A?@?B>A>BA?>=>B>?C=>B@B=>??>><?<=<;@:<;?>:=:<?=>9:9:>;<>;=9>8<:=7;=:;
@MGYG000002096_86|kraken:taxid|3024_read_11
AGAAAATGGGCAACAATATAGGAAAGGATGAACACTGTTTGAAGAGAGGAGCACGCGGATGGGGAAAAATATGGA
+
?@E@@BBB>ABCB?CB????=>>?BBA?=@A?<?<A@@@;>A??;:;@?=@?9<999<9==>>9;==>8:888<;
@MGYG000003260_60|kraken:taxid|4334_read_12
GGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGGCGACGAGG
+
DBA?CBCABCA@A?@>AB==AB>?AB>@A=<<B=A?;=;;=?;=;;>>;@<>:><>;<=?><<<9;>9>;:;=77
@MGYG000002096_86|kraken:taxid|3024_read_13
GAGGGACTGGAAATTATCATTGCCCATCCCGAGCGATATCAATATGTGCCCTGGGGATCGCTCTGGACTGGCAGG
+
ABC@D@ABBCDD>B=>@C=@@>CCABB>A@@?@@=<@A=;;>@;;@=??@;?9=>??=;;:9:8:<<=;8=8<9:
"""

    ref_file_path = "ref_test.fa"
    reads_file_path = "reads_test.fq"
    reference_path = "kmers_db.kdb"
    aln_output_path = "alignment.aln"
    k = 31
    with open(ref_file_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")

    with open(reads_file_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    ref_dict = generate_ref_dict(ref_file_path)
    ref_kmers_db, genome_bases = build_ref_kmers_db(ref_dict, k)
    with gzip.open(reference_path, "wb") as file:
        pickle.dump(ref_kmers_db, file)
        pickle.dump(genome_bases, file)
    execute_align(reference_path, aln_output_path, reads_file_path,
                  unique_thresh, ambig_thresh, mrq, mkq, mg, rev_comp, False,
                  False, None, None)
    with gzip.open(aln_output_path, "rb") as file:
        aln: Alignment = pickle.load(file)

    results = []
    for read in aln.get_reads().values():
        results.append(cast(Read, read).get_status())
    assert results == ['Unmapped', 'Unmapped',
                       'Unmapped', 'Unmapped', 'Unique', 'Ambiguous',
                       'Unmapped', 'Unmapped', 'Unmapped', 'Unique',
                       'Unmapped', 'Unmapped', 'Unique', 'Unique']
    os.remove(ref_file_path)
    os.remove(reads_file_path)
    os.remove(reference_path)


def test_read_print():
    """This func test the print of a Read object"""
    read = Read("MY_READ123", "AAAAAAAAAAAA", "AAAAAAAAAAAA")

    def create_read_print(read: Read):
        """This func return a string to print a Read object"""
        if read.get_status() == "":
            status = "NOT DEFINED"
            return (f"Read id: {read.get_header()} \n"
                    f"data: {read.get_data()} \n" +
                    f"quality: {read.get_quality()} \n" +
                    f"read mapping status: {status}")
        else:
            return (f"Read id: {read.get_header()} \n"
                    f"data: {read.get_data()} \n" +
                    f"quality: {read.get_quality()} \n" +
                    f"read mapping status: {read.get_status()}")

    assert str(read) == create_read_print(read)
    read.set_status("Unmapped")
    assert str(read) == create_read_print(read)


def test_ascending_sort_refs():
    """This func test the function ascending_sort_refs() used for the similarity
    extension, this test is making sure the genomes provided are sorted as expected
    and return the list of sorted genomes.
    In this example snail has the largest num of unique kmers so its last, after his
    bird because of smaller unique Kmers,cat is the longest sequence from the
    remains and dog and owl have the same length yet dog was inserted to DB first
    so it's the first and owl second"""
    ref1 = Reference("MGYG014_62|snail:taxid|3633", "ATGCCGGGG")
    ref2 = Reference("MGYG061_62|bird:taxid|2211", "GCCGGGGCTAA")
    ref3 = Reference("MGYG061_41|dog:taxid|2213", "CCGGG")
    ref4 = Reference("MGYG061_11|cat:taxid|2214", "GGGCTAA")
    ref5 = Reference("MGYG061_34|owl:taxid|4455", "GCTAA")

    refs_dict: Dict[str, Reference] = {ref1.get_header(): ref1,
                                       ref2.get_header(): ref2,
                                       ref3.get_header(): ref3,
                                       ref4.get_header(): ref4,
                                       ref5.get_header(): ref5
                                       }
    ref_kmers_db, genome_bases = build_ref_kmers_db(refs_dict, 4)
    sorted_refs = ascending_sort_refs(refs_dict, ref_kmers_db)
    assert sorted_refs == [ref3, ref5, ref4, ref2, ref1]


def test_similarity_extension(monkeypatch, capsys):
    """This function creates FASTA file with two identical genomes and supposed
    to filter out 'ELAL0000000_131|filter_me_please:taxid|3322' genome from
    Reference Kmers DB, func check the expected output is printed to screen"""
    ref_content = """>MGYG000004169_19|kraken:taxid|5353
TGGTGGGAGACCGTGCAGGTTCAAGTCCTGTTAACCGCAGTAAATTGTGAAAAAGTGGCG
GTATCGTGAATGATACCGTCACTTTTTTCGTATACCTTTTTGTCCGGTTCTTGGTTAAGT
TGGATTGGCAGAGGGCTTGTATATAGAAGCTTTGCAACGATACGTTGCGCAGATTCCTAC
ATGGGCTGGTAGATGTGAAATAGTTTAGGCGTTATGATGATTCAAAATTAATTGCAGGAG
GAGTGGTCGAATGGGCAAAGGGTTAATTCGGTTCGTAATCATTTTTAGGAGTATGGGATG
CTTCTGAAGTGCTAAGGAGTCCTTCTTTAAATGTTTTTACCCAGCAGTAAATCGTGCCAT
TCAGGACATTCTCCTTTCTGCTCCGATGGCTTGACTGGTATCAAGGATGCAATCTCTACA
GCATTTCCAAAAACGGAGCAACAGCGTTGTATTGTACATATGGTGAGAAACACGCTTAAA
TACGTTGCAAACAAGGACATGAAAGCATTTGCAAAAGATTTAAAGACAATCTATACCGCT
GCAGATGAAGAAGCCGCCAGAAAGCAGTTGGAGTCTGTAACAGAAAAATGGTCTGCCCAG
TATCCAAGTGCGATGAATCGCTGAATTCCTGCTATCGTAGATTGAACAAGCAACGCAGTG
TATTCCAAAGCTCTCAGGCGCTTATGAAAGCCCTATATCTGGGAACCTTCGAGATTGCAA
>ELAL0000000_131|filter_me_please:taxid|3322
TGGTGGGAGACCGTGCAGGTTCAAGTCCTGTTAACCGCAGTAAATTGTGAAAAAGTGGCG
GTATCGTGAATGATACCGTCACTTTTTTCGTATACCTTTTTGTCCGGTTCTTGGTTAAGT
TGGATTGGCAGAGGGCTTGTATATAGAAGCTTTGCAACGATACGTTGCGCAGATTCCTAC
ATGGGCTGGTAGATGTGAAATAGTTTAGGCGTTATGATGATTCAAAATTAATTGCAGGAG
GAGTGGTCGAATGGGCAAAGGGTTAATTCGGTTCGTAATCATTTTTAGGAGTATGGGATG
CTTCTGAAGTGCTAAGGAGTCCTTCTTTAAATGTTTTTACCCAGCAGTAAATCGTGCCAT
TCAGGACATTCTCCTTTCTGCTCCGATGGCTTGACTGGTATCAAGGATGCAATCTCTACA
TACGTTGCAAACAAGGACATGAAAGCATTTGCAAAAGATTTAAAGACAATCTATACCGCT
GCAGATGAAGAAGCCGCCAGAAAGCAGTTGGAGTCTGTAACAGAAAAATGGTCTGCCCAG
TATCCAAGTGCGATGAATCGCTGAATTCCTGCTATCGTAGATTGAACAAGCAACGCAGTG
TATTCCAAAGCTCTCAGGCGCTTATGAAAGCCCTATATCTGGGAACCTTCGAGATTGCAA
>MGYG000004832_79|kraken:taxid|6050
CCGCCCGGGTGAGCCGGTCCCGGGGGTCCAGGGCCGTCTGGGCCACCAGGGCCTCGTAGG
CCCCGCAGATGAGCCCTAAGTCGGTGAGCTTGTCCCCGCTGGGGCCCTCGGCCTCCTCCC
CCGCCCGGATGAGCTGCTCCGGGCTTACCCGGCAGCTTTTCAGCTCGTCCACCGTGGCCA
GCAGGGATTGGAGGAAGGCCGGACGCTTGGAGGGGCGGCCATAGACCTTCAGCTGCTGGG
CGGTCTCCTGCACCGCCCGATACATCAGCAGCAGCCGTCCGCCCCCGTCCAGCTCCTCCT
GGCCCAGTCCGCCCGCGGCCTGGAACACCCGGTTGGCCAGGCGGCTGAAGGACAGCACCT
CTGCCCGCAGGGAGATCCCTGGGCCGCCAGCCTTACATAAGGCCCGCTCGGCCTCGTGGG
ACTGCTGCTCCGGGGTCATCAGCACCTGGGGCCGCTCCCCGGCGGCCTGGCACAGGCGGT
TTAGCACGGCGGTGGTCTTGCCGCTCCCGGCCCGTCCCATGAGAATGCGCAGCATCCGGC
AGCCCTGCCTCCACGTCTCCCCCGGCGTTGTTCAGCAGGATGAGCAGGCCGTGGATGTCC
GGGTCGCCCTCCACCGCCGCCAGCAAGGGCAGCACGTGCTCGTATTTGGTGCTCTTCATG
CCGGAGGGCAGCAGCTGGTGCCCCTCGATCTGCCCCACGATGGTCAGCACATAGACCCGC
CCCTCCTCCCGCTGGATGAGGGAGGTACCCAGGTCCACCAGGTTCTGCGTGTTGGCCGGG
GCCTCCGGCTGTGTGTGCTCCTGTTCCATGTCTGCCATTTCTTGCGCCTCCTAAACCGGA
ATTTCTGTTAGGGTGCGCGGAATAGGGCCGGTTTACTCACAAAAATGGGCGCGCGGCGTG
ACGAAAATTCGTCACGCCGCTGGTGTTCCGGCTAAAATCAGCCCCGCTCACGCATCTTGG
CGAAGCAGTCGCTGCAATAGACAGGACGGTCAGACTTGGGCTCGAAGGGCACCTTGGCCT
CGCCGCCGCAGGCAGCGCAGACAGCGGTGAAGTACTCACGGGGGCCACGAGCCGCGTTCT
TGCGGGCATCACGGCAGGCCTTGCAGCGCTGGGGCTCGTTCTGGAAGCCGCGCTCGGC
>MGYG000002096_86|kraken:taxid|3024
GGGGCATAGACATAGCCATACTGCGCTTTTCCGCTCAAAGAATATCCTTTGATAAAACGC
CCTTCCCAGCGGCCGTCTTTCCGTTTGCGAATATTTTCTCCTCTTCTTGGCATGCGAACT
CCTCCTTACAAGTATATGTTGTGTTTTGACCCCTGCATTTGACCACGAATAGACCCCCGC
ATGCAAACAGTCCCGTATATTGCCGATATCATCGGGAGCAGAATTGACAAATATAGCGGG
CAAGGCTTATAATTAGCATAGCTTTTGCATAGGGCAACGCATAAATAAAAGAGGTTTCAG
CAGGCAAATATCCTACCTATAGCAATATACTTTCGGAAAGAAAATTGAGCGGCTGGGTAA
TTTTGTTCTGTGGCCAAAAAAGGGAGAGGATAAAATCAGCGGTCGCTGCAAAAAGAAAAA
TGGCGCAGTGGTATTATACCATGCGGAAACAAACGGATTTTCAGGGCTCATAGTTTGCAG
AAAATGGGCAACAATATAGGAAAGGATGAACACTGTTTGAAGAGAGGAGCACGCGGATGG
GGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCA
TCGGCGGCGACGCTTCATATTCTGGGAACGGCCCCGGGCAGCAATGATCTTCAGGAAAAT
TGTGAAATTTTATTGATTTCAGATAACAAAACCGGCGCGGTTTTGACGGAAAAGAACGCG
GGCGAAAAGATGACTCTCATCGGCGGCACGGTCCGCTGGATGGTGGCTCTTACCGCCGTC
GATCACCTGGAATTATCCGAGGAGGTGACTCTGGAGGCGGCAGATCTCGAGCCCTTTCAG
GGAAACCGCAAAATCGGCCTTCGGGCGGGACAGACTCTGACGGTCAGGGACCTCCTGGCG
GCCATGATGGTGGACTGCGCCCAGGATGCCGCCGTGGCGCTGGCAAACAAGGCGGCGCAG
AAGGCGGGAGCGGATGATTTTGTCGCTCTGATGAATGAGAAGGCAGCCGAACTGGGAATG
AAGGATACAACATTTAAAAATGCCACCGGCGCGCTGGACAGCGGGCAGGTGACAACGGCG
AAATTGTTCCCTCCTATCAATGGACGGGGGAGAGCGATACCTCCCAGAAACAAGAGGCAT
AAATGCGGCAGGCTGGACCAATTGCCGGCTGCAAATGTGGCTAAAATGAAGAAAGCAGCA
ATATTTATGTGTTCTCATACAGAAAAAGAGAACCTGAATCAACGGCATTTGGAAGCAATA
TGCAATGGAGTAAAGCAGTGTGAATCTTATGTCCGCCTGTGTAGCCGCAGCCTTAAAAAG
GCCAAGTGGATTTAACCGCGGGAAGAACTGCACAAGCAGTGAGCATAGTTGGGGCCGGAT
TATGTTGGGCTATGCTGGGCGGCTAATGAACGGCGGGGGAAAATGAAAAATTTTCTATAT
AGCAGCATCTATTCGATGTAAAGAGGCGGTTTCTCAATTTCCCAATTAGGTATGTAATGA
TCTGTGGGCCAATGGCATGAAATAGTTATAAAGCGAAAATAAACTTGAGGGAAAGCAATT
CTTATGCCTCATCAGGGAAATAATTGATCTGGCGCAGTGGGGGAAAATAGCTGCAGCCAA
AGCTGGATAAATTATATAAATAAATGTACTCCGGGATGTTGGGATA
>MGYG000002720_43|kraken:taxid|3759
GTTCTCCGGCCGAGCTGCTGAAAAAGAACGGCCTGTTCCGCCATATGGCGCAGCTTCAGA
GCGAGAGCATGGAGTGGACGGCGTAGAAAGTGCTGAGATAAACAGGATTCCAAGTGCAAT
ACTGGAAAACGAAGTGTCGAAAACAAGGAGGCGTTTATGTAAATAAATGGTCTGTCAGGG
ATGTGATCACCACGGTATTGCTCTCCGCAGTACTTATCGTGATCCAGCTCGTTGTCAATA
TGGTCTGCATGGCCAACGACTTTGTCAGCATGGTGCTGTCGGTGGGGATCACCATGTTCC
TCTGCGCGCCGGTCTATATGCTCATGGTCAGCCGGATCGGCAAGCGGTTTGTGACGCTGA
TCTATATGACGCTGCTTGGCGTGATTTTCCTGCTGATGGGAAACTGGTTTTTGCTTCCCT
ACTTTATCGTTACAGGCGTCATCTGCGAGGCGATTCTCTGGAAGGAAGGCTCCTGCCAAA
AGCCGAAACGGCTGACCGCCGCCTGGACGGTGGCGAGCCTGCTGTACAACGGGGTCAATC
TGCTCCCCATCTGGTTCTTCTGGGACACCTACTACGATTTTGCACTGGCAAGCGGCATGG
AGCAGAGCTATATCGATTCCTATGTGCGTTACTACACCTCTCCCGGCTGGCTGGCCTTTA
TTCTGCTGTTTACGACGCTGATGGGCTTTTTAGGCTGCATGGTGGGCAGTCGGCTGATCC
GCAGGCATTTTAAGAAGGCCGGCGTTCTATGAGGGCGGACGCATCCTTTGCGGTGCCGGT
CAAGCTGTGGGCGCTGCTGTGCGTCTTTGCCGGAGTAACCATCGGCGGGAATGTGCTGCT
CACCTGCATCCTGACCGGCGGGGCGCTTCTTTATCTCGTCCTGCAGCGGAACTTCCGCCT
TGCCGCGTCCTATGGCTGCTTTTATCTGCTGCTGGCGCTGCTGCTGTACGGCATCCGCTT
CCACGGGCTGCACATGCCGGTCTTTTCGGAATTTTATGTTCTGATGTTCTGGAATCTGTC
ACCGATCTTTCTGGTGTCATGGGATCTGATCACCACGCCGCCGGGGATGCTCTCCGCGTT
TTTATCCCGCCTCCGTATGCCCACCCCGTTCATTCTCGGGCTTCTTGTGGTCTTTCGCTT
CTTCCCCACCATGCGGGCGGAGCTAAAGGGTGTAGGCAGGTCCATGAAAAACCGCGGGCT
GACCGCGGCGGGGCAGCTGCTCACGCATCCCGTGCAGAGCATGGAGTATGTGCTGGTTCC
GTTTCTCCTTCGGGTGTTGCAGCTTGCCGATCAGCTTTCCGTTTCCGCTGTTGCAAGAGG
CGCGGAGCGTCCGGGCGTGCGCGGCAGCTATTATGAGAAAAGGGGCGGGGCGCGGGATCA
CATCGCGGCAGCCGCATGTGCGCTGGTCACGGCTTCCTATTTGGTTTTGGAAAGGAGCAT
GGCATGATCGAGCTTACACGCGCTTCTTTTCAATACGAAAACTCAGACAGAGGCGTGCGG
GACATTTCCCTTTCGGTAAAAAGCGGCGAATGCGTCGTGCTGACGGGGCTCTCCGGCTGC
GGCAAGACCACCGTGACAAGGCTTGTCAACGGGCTTGCGCCATTCTACTATCCCGGAGTT
TTCAGCGGGAGCGTGAGGATAGACGGCAAGGACATATCAAGGCTTTCCACATGGGAGATC
GGGCGTTTGGTGGGAAGCGTTTTTCAGGACCCCAAAAGCCAGTTCTTTTCCTCCGAGCTT
CATAACACGATGCCGCCCCGATGTCAATATGCAAATTCGTCGCAAATTAGCCGCCAACGC
GGATAGTGCCCGCTCCGCTGAAATGGGCTGCGGCAGAGCTGAATGATCTGCGCTCCGTTT
TCGGCGGCAAAAAGGCAGCGGCAAATGCCAGGTGAAAAAGCGGTAGAAAAACCTGAAAAG
ATATGCTTCACACTTATCAAACCCCAAAAGCTCTTTCAGAGCAAAAGGATTTCATTACGA
ACAGCAGCTGACCCTCCGTATGTACTGAAGCATACGGAGGGTCAATCGGCTTGTGTGGCC
TGCTCAGTCGAGCAGGGGCTTTCTTTTCTGATGCGTCAACATATCCGTGCCCGGCCTATC
CCTCAAACAAATTCAGGAGGTATCCGTAGCCCTCGGCCTCCATTTTCGCGTAAGGCACGA
AGCGGAGGGCGGCGCTGTTGATGCAGTAGCGCACACCGTTGGGCGACTCCGGGTCGCCGG
TGAACACATGACCCAAATGGGAGTCACCGGCGCGGCTGCGCACCTCTGTGCGGCGCATGC
CGTGGCTCAAGTCCTCCCGTTCCACCACAGCAGGACTCTCGATGGGCTTTGTGAAGGCAG
GCCAGCCGCAGCCGCTTTCAAACTTATCCGTGGAGGAGAACAGCGGCTCGCCGGTGACGA
TGTCCACATAGATGCCCTTTTCAAACTTGTCCCAGAGCTCACCGGTAAAGGCGCGCTCCG
TGCCGCTCTCCTGCGTGATGCGGTACTGCTCCGCCGTCAGCGTGTCCCGGATGGTCTCCG
CCGCGGGCTTCCGGTAGTCGCCCGGATCGATGCGAAGCCTGGAGAACAGCTCCATCTCCG
CCCGGGGGATGTGGCAGTAGCCGTTGGGGTTCTTTTCTAGGTAGTTTTGGTGATATTCCT
CGGCGGGATAGTAATTCTTCAGCGGACCGATCTCCACAAAGAATTTCTCGCTGCGGCCCC
GCTCAATGCCTGCGATGCGCTCCACCGTTTCCCTGACCTTGTCGTTGGTGTAGTACACGC
CGGTCTGATACTGGCTGCCCCGGTCGTTGCCCTGCCGGTTCTGCACCGTGGGGTCGATCA
CATAGAAGTAAGCCAGCAGCAGCGCGTCGAGACTCACCTGCTCGGGGTCGTATTCCACCC
GCACGGTCTCCCGGAAGCCGGTGTCGCCCTTGCAGACGGTCTGGTAGGTGGCGTCTGCCT
CGCAGGTGCCGTTGGCATAGCCGCTTTGGGCGTCGATGACGCCGGGGATGGATTGCATGA
GCTGCTCCATGCCCCAGAAGCAGCCGCCTGCCAGATAGATGACATTTTCCGTGGTGTCCA
TATATATACCCTCCTCGCCGGACATATCCGAGGACTGTCCGGTTTTCTCCGCTTCCTGCG
TGTCCTTTTCCCGCGGTACACGCTGCCGTTCTGCGGCGCAGCCGCTGAGCAGCAGTACCG
AACTTCGGCAGGCTGCTGTCATCCGTCTTGCTGCCCGTGCTGCCGCTATCGGTCTTATTG
CCTGTGTTTCCGCTGTCGGTTTTATTGTCCGTACTGCCGCTGTCGGTTTTGTTGTCCTTA
TTGCCGCCGTCGGGCATCGTGTCCATATCCTTGTCCGTCGGCATATTCATATCGCTGTCC
AGAGATTTTTGCATGATCTCGGGATACTTTTCCTCCAGTGCGGCCAGCTTGTTTCCGATG
TCCACGCCGTGGATATATGCCTTTACCAGCTCCCGCGCCGCCTGCCGCATCCGGCCCGCA
GATGAACGGCAGCCAGAGCAACCACATGGCGATGGGAGGGCTTTCACGATGATCATTCAA
GCAGAACTGAAGTGTAAGCAGACCGGGTGCGAGGCAGGCCCCTGTGCCGTGGATAAGGTC
ATCGAGCTGCCAAGCCCGCGGTTTCAGCAGTTCAGCCGCGCACTGCTGGCTGATTACGAT
TTCATCGCAGAGAATAAAAACGCCATCCGGCACGATGAGGATGCCAGGCACTGTCTGCTC
ATCCTCGACGCGGATGGAACGGACGGTTTCCTTGTTGACCCACAGGGGCACAACTACGCC
CGGTACAGTGCCTTTGTCCCCAACGCCCGCAGTTTGCTGACGCCGGATATGGGGATCGAC
CGCAGCTATCTTTCGCCGGCAGAGCCTT
>MGYG000000975_36|kraken:taxid|1603
TGACGTTGTATTCTTTGTTTGTACATCTATTTAGTTTTCAAAGATCTCTTGCACTAATTG
ATTCTATCAATCTTAACGAAGAATTAATCCTTCGAACAAGTGCTTAATAATTATATCTAA
GAGATATAATCATGTCAACAAAGAAGTCTCAAATCTTGCGACAACTTCAAAATGATACAC
CTGTTCGATATTGTAAATCAAGCACAAACGAAATAGTATCTCTATCCTTATTATATATAT
AAGAGTTTTTTATTCCTGAAAAGTCCTCCAAAAAAGACAAATATTCATTTTCAATTTCTC
ATACATAATCCATGTGTCAACCGTTCCAAACATAAGTTCGCCTTTTTCAGCTCTTTCTTG
TCCATCTTCAATATTATCTAAAATAAATCTTATCTTTGAAGCACTGAAATAAGGGTTAAT
TATAAGTCCAGGCAAATTGGAGCAGATTGCCTGGATTGCCAAACAATCGCATTATAAACA
GGCATACCAGTTTCTTTTTCCCAAATAATTGTTGTTTCTCTTTGGTTTGTGATACCAACA
GAATCGATGTCATCATATGTTAAATTAGTTTTTAATAAAATATCATTGATAACATAAATT
ACACTTAAATATAAATCCAAAGCGTTTTGTTCAACCCAACCACTTCGTGGAAAAAGACAT
TCTACTTCTTGTTGAGCTTTGAAGAAAACTTTTCCTTCATGATCAAAAAGCATCGCACGA
GTTGATGTCGTGCCTTGGTCAATTGTTAAGATGTACTTAATAGTTAATTAAAACTTTAAA
TCTGAATTTAAAAACGAATAATAAGCTAAATTCCCACCTTTAGAATGTCCTAAAACTATA
CTCTACGAGGATCTTTGAGCTCTTGCTTAACATATTCATGTCCTCTTTTACATGTGTAAC
CAGTAATTTCATCATCTTTAACATGAATATGGCAGCCTTTAGGGCAAACAATACAAATAA
AATCTTTTTCCATATTAAAGTTCCTCAATTTCTAAACTTAAAGAATTCAAAGATTTATAA
TCAACATCTTTTACTGAAATCTTTTGCATTTCTGATGGAATTAAATCTACCTTCATAATT
TGCATCACTTAAATCATCAATATTATTAAGCAAAACAAAAATGTCTTCTCTATCAATAGG
AGTAATAAATTCTTTTGCTAATTGTTCTTCAACATCATCAGCTTGATGTTCAATTGCATG
AACACTATTTTTTAATTCTAATAATTTAGAATTATCAAAATTCAAAAGTCCTTCGCTTAA
TACGTCAATTTCGTTTTTAGAATATTCAATCAAATCTTTAAATGATTGATAGAAATAATT
CTTTTCTTTTTTCTTAAAAAACATATAAAACTCCTTACCTTGCTATAGCAATAAATATAA
AAAATGGCCTTGTCCTATTTCTATTCATTTTTCTACAAATTAATTGAATAGTTTTAGTAA
CTAGCCATCCTACAAAAAATCCAAATACACAACTAAGAATTAAACCATAAATTGTAAATG
CCCAGCTTGAATCTGGACCAAATTTTATAGTTGCATCTAATCCCATTGTCATTGCAGCAA
GTCCTGCTCCTGTAATTCCAGCAATTAATGCATGTGATTCACTACTTGGAATACCAAACC
ACCAAGCAAGTACTCCCCAAACGATAACACCAACCATGCCTGCACAAAGTGCACTTAAAG
CAAGATTGTAATTAGTACCAAAATCAACCATGTTAGAAATAGTTTCAGCAGTTTGCTTAG
ATATAAAAAACATAGTTAAAACACCTAAGAAATTGAATATTGCCGCCATTATAATAGCTG
TTTTTGGCTTGATACATCTTGTTGAAATGCAAGTTGCAATTGCGTTAGGTGCATCTGTCC
AACCGTTAATAAAAATTACACCAATCACTAGGATAGTTGCAACGCTAAGCATTGGATAGG
AACCAAATAACGCAAATATTTGTTCTAAGCTCATACCTAAAAATTATAAATCAGTGCTTT
TTAGAATGTAAAGTTTTTCCTTTATATTTGATAGCGTTTCATTATTGTAGGTAACAAAAC
AAAAGAGAACTCGTATAATTAAGC
>MGYG000001129_3|kraken:taxid|1782
CTCTTTAAACGATCGGCAATGGAAATGCGGCGATCAGCTGCCATGAAATAGCCGGCATGG
GTTTGTAACATTTTGGAAAGGCTTGCCGTTACAGCACCGCCGAAAATAGTAATCAATTCA
CAGCCAAGTGCAAGCCATGCCGTTTTGCTGCCTTGTTTTCCATCGGTCAGTGCCAGTACA
ACAAAATAGATTGCTGCCACTTGGAACATTTGCAGGATTGCATTTAAAAATCCAACTGCC
ACAGACTGATTGATGTTTTTTGTTTCTTTTCCGGCAAACTGCCAGATTTTACGGAATGTT
TCAATCATTTTGTGTCACCGTCCTTTGCTCCAATATGCGCCTGCCACATTTCATGGTAAA
GTGCACAATTTTTCAATAGCTCTTCATGCGTTCCCTGCGCTGCAATATGACCATCCTGTA
GGACAACAATGTTGTTGGCATCCGTGATAGTGGAAAGGCGGTGCGCAATCACAATGACGG
CGGATGTGCCGCCGTCGGAATATTCCGGCATCTGAGCAGCCTGTGCAGCACGGTCCATTT
GGTGAAAAGTTTCTTGAATCTTGAGAAAGCGTTCGATCATTTCACTGAGCACAGGCTCAA
CGATTGCGGAAAAAAGCCCTTCCTTATCTCCAAAGCGGACATAAATGGAGTTGGTACTGG
TATTGGCAGCTGTTGCAATGGTGCGCAAAGAGGCATCTACATAACCTTTTTCTAAAAATT
CCTGTTCGGCAGCTGTAAGAATGCGTTCTGCAACACCTTCAATTTGTTTTGCCATTTGTT
TCGCCTCCTTTCATTCATAACAGCGTTTCATAACACTGTTTTTATTATATACCAGCAAAT
TTAAAAAGTCAATACCTCAGGATGGATTTCGTTGGAAATCAGTTGAAAATCCGTATATGC
TATGATATAATGAGAAAAA"""

    ref_path = "test_extsim.fa"
    with open(ref_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    mock_argv = ["main.py", "-t", "dumpref", '-g', ref_path, '-k', "31",
                 '--filter-similar']
    monkeypatch.setattr(sys, 'argv', mock_argv)
    main()
    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)
    if len(outputs) != 2:
        pytest.fail(
            "Expected two separate JSON outputs, but got something else.")

    similarity = {'Similarity': {
        'ELAL0000000_131|filter_me_please:taxid|3322': {'kept': 'no',
                                                        'unique_kmers': 29,
                                                        'total_kmers': 630,
                                                        'genome_length': 660,
                                                        'similar_to': 'MGYG000004169_19|kraken:taxid|5353',
                                                        'similarity_score': 0.953968253968254},
        'MGYG000004169_19|kraken:taxid|5353': {'kept': 'yes',
                                               'unique_kmers': 89,
                                               'total_kmers': 690,
                                               'genome_length': 720,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000001129_3|kraken:taxid|1782': {'kept': 'yes',
                                              'unique_kmers': 889,
                                              'total_kmers': 889,
                                              'genome_length': 919,
                                              'similar_to': 'NA',
                                              'similarity_score': 'NA'},
        'MGYG000004832_79|kraken:taxid|6050': {'kept': 'yes',
                                               'unique_kmers': 1108,
                                               'total_kmers': 1108,
                                               'genome_length': 1138,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000002096_86|kraken:taxid|3024': {'kept': 'yes',
                                               'unique_kmers': 1636,
                                               'total_kmers': 1636,
                                               'genome_length': 1666,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000000975_36|kraken:taxid|1603': {'kept': 'yes',
                                               'unique_kmers': 2034,
                                               'total_kmers': 2034,
                                               'genome_length': 2064,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'},
        'MGYG000002720_43|kraken:taxid|3759': {'kept': 'yes',
                                               'unique_kmers': 3838,
                                               'total_kmers': 3838,
                                               'genome_length': 3868,
                                               'similar_to': 'NA',
                                               'similarity_score': 'NA'}}}

    summary_output = json.loads(outputs[0])
    similarity_output = json.loads(outputs[1])
    assert similarity_output == similarity
    os.remove(ref_path)


def test_coverage_with_overlapping_kmers(monkeypatch, capsys):
    """This test creates a FASTA AND FASTQ files and compare the resulted output
    from execute_build_ref_and_align_and_dump(). in the reads FASTQ file - some
    kmers overlap each other and supposed to be counted once"""
    k = 3
    ref_content = """>MGYG000004832_79|kraken:taxid|6050
    AAACTCAAACCC
    >MGYG000004832_78|kraken:taxid|1111
    AACCTTTTTTT
    """
    reads_content = """@MYREADS0000| read_1
    AAAAAAAAAAA
    +
    !!!!!!!!!!!
    @MYREADS0000| read_2
    CCCAAA
    +
    ??????
    @MYREADS0000| read_3
    ACCT
    +
    ?!?!"""
    reads_path = "test_coverage.fq"
    genome_path = "test_coverage.fa"
    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    execute_build_ref_and_align_and_dump(genome_path, k, reads_path, 1, 1,
                                         None, None, None, False, True, False,
                                         1, [
                                             "MGYG000004832_79|kraken:taxid|6050"])

    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)

    summary_output = json.loads(outputs[0])
    coverage_output = json.loads(outputs[1])
    expected_coverage = {
        "Coverage": {
            "MGYG000004832_79|kraken:taxid|6050": {
                "covered_bases_unique": 10,
                "covered_bases_ambiguous": 0,
                "mean_coverage_unique": 1.3,
                "mean_coverage_ambiguous": 0.0
            }
        }
    }
    expected_summary = {
        "Statistics": {
            "unique_mapped_reads": 3, "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0},
        "Summary":
            {"MGYG000004832_79|kraken:taxid|6050":
                 {"unique_reads": 2, "ambiguous_reads": 0},
             "MGYG000004832_78|kraken:taxid|1111":
                 {"unique_reads": 1, "ambiguous_reads": 0}}}
    assert coverage_output == expected_coverage
    assert expected_summary == summary_output
    os.remove(genome_path)
    os.remove(reads_path)


def test_coverage_with_overlapping_full_cov(monkeypatch, capsys):
    """This test creates a FASTA AND FASTQ files and compare the resulted output
    from execute_build_ref_and_align_and_dump(). in the reads FASTQ file - some
    kmers overlap each other and supposed to be counted once. Func run full coverage"""
    k = 3
    ref_content = """>MGYG000004832_79|kraken:taxid|6050
    AAACTCAAACCC
    >MGYG000004832_78|kraken:taxid|1111
    AACCTTTTTTT
    """
    reads_content = """@MYREADS0000| read_1
    AAAAAAAAAAA
    +
    !!!!!!!!!!!
    @MYREADS0000| read_2
    CCCAAA
    +
    ??????
    @MYREADS0000| read_3
    ACCT
    +
    ?!?!"""
    reads_path = "test_coverage.fq"
    genome_path = "test_coverage.fa"
    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    execute_build_ref_and_align_and_dump(genome_path, k, reads_path, 1, 1,
                                         None, None, None, False, True, True,
                                         1, [])

    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)

    summary_output = json.loads(outputs[0])
    coverage_output = json.loads(outputs[1])
    expected_coverage = {"Coverage": {
        "MGYG000004832_79|kraken:taxid|6050": {"covered_bases_unique": 12,
                                               "covered_bases_ambiguous": 0,
                                               "mean_coverage_unique": 1.9,
                                               "mean_coverage_ambiguous": 0.0},
        "MGYG000004832_78|kraken:taxid|1111": {"covered_bases_unique": 11,
                                               "covered_bases_ambiguous": 0,
                                               "mean_coverage_unique": 2.1,
                                               "mean_coverage_ambiguous": 0.0}},
        "Details": {"MGYG000004832_79|kraken:taxid|6050": {
            "unique_cov": [2, 3, 3, 1, 1, 1, 2, 2, 3, 2, 2, 1],
            "ambiguous_cov": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]},
            "MGYG000004832_78|kraken:taxid|1111": {
                "unique_cov": [2, 3, 3, 1, 1, 1, 2, 2, 3, 2, 2, 1],
                "ambiguous_cov": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}}}
    expected_summary = {
        "Statistics": {
            "unique_mapped_reads": 3, "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0},
        "Summary":
            {"MGYG000004832_79|kraken:taxid|6050":
                 {"unique_reads": 2, "ambiguous_reads": 0},
             "MGYG000004832_78|kraken:taxid|1111":
                 {"unique_reads": 1, "ambiguous_reads": 0}}}
    assert coverage_output == expected_coverage
    assert expected_summary == summary_output
    os.remove(genome_path)
    os.remove(reads_path)


def test_coverage_no_overlap(monkeypatch, capsys):
    """This test creates FASTA and FASTQ files and test the function
    execute_build_ref_and_align_and_dump(), the kmers in the reads FASTQ don't
    overlap."""
    k = 3
    ref_content = """>MGYG000004832_79|kraken:taxid|6050
    AAACTCAAACCC
    >MGYG000004832_78|kraken:taxid|1111
    AACCTTTTTTT
    """
    reads_content = """@MYREADS0000| read_1
    AAA
    +
    !!!
    @MYREADS0000| read_2
    CCC
    +
    ???
    @MYREADS0000| read_3
    ACCT
    +
    ?!?!"""
    reads_path = "test_coverage.fq"
    genome_path = "test_coverage.fa"
    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    mock_argv = ["main.py", "-t", "dumpalign", '-g', genome_path, '-k', "3",
                 '--reads', reads_path,
                 '--coverage', '--genomes',
                 "MGYG000004832_79|kraken:taxid|6050"]

    monkeypatch.setattr(sys, 'argv', mock_argv)
    args = readargs()
    execute_build_ref_and_align_and_dump(genome_path, k, reads_path, 1, 1,
                                         None, None, None, False, True, False,
                                         1, [
                                             "MGYG000004832_79|kraken:taxid|6050"])
    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)
    if len(outputs) != 2:
        pytest.fail(
            "Expected two separate JSON outputs, but got something else.")
    summary_output = json.loads(outputs[0])
    coverage_output = json.loads(outputs[1])
    expected_coverage = {
        "Coverage": {
            "MGYG000004832_79|kraken:taxid|6050": {
                "covered_bases_unique": 9,
                "covered_bases_ambiguous": 0,
                "mean_coverage_unique": 0.8,
                "mean_coverage_ambiguous": 0.0
            }
        }
    }
    expected_summary = {
        "Statistics": {
            "unique_mapped_reads": 3, "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0},
        "Summary":
            {"MGYG000004832_79|kraken:taxid|6050":
                 {"unique_reads": 2, "ambiguous_reads": 0},
             "MGYG000004832_78|kraken:taxid|1111":
                 {"unique_reads": 1, "ambiguous_reads": 0}}}
    assert coverage_output == expected_coverage
    assert expected_summary == summary_output
    os.remove(genome_path)
    os.remove(reads_path)


def test_coverage_extension2(monkeypatch, capsys):
    """This test creates FASTA and FASTQ files to test execute_align_and_dump(),
    it requests one genome to cover and compare expected printed result"""
    k = 3
    ref_content = """>MGYG000004832_79|kraken:taxid|6050
    AAACTCAAACCC
    >MGYG000004832_78|kraken:taxid|1111
    AACCTTTTTTT
    """
    reads_content = """@MYREADS0000| read_1
    AAA
    +
    !!!
    @MYREADS0000| read_2
    CCC
    +
    ???
    @MYREADS0000| read_3
    ACCT
    +
    ?!?!"""
    reads_path = "test_coverage.fq"
    genome_path = "test_coverage.fa"
    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")

    refs_dict = generate_ref_dict(genome_path)
    refs_kmer_db, genome_bases = build_ref_kmers_db(refs_dict, k)
    kdb_file_path = "test_coverage.kdb"
    dump_ref_kmers_db(refs_kmer_db, genome_bases, kdb_file_path)
    execute_align_and_dump(kdb_file_path, reads_path, 1, 1,
                           None, None, None, False, True, False,
                           1, [
                               "MGYG000004832_79|kraken:taxid|6050"])
    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)
    if len(outputs) != 2:
        pytest.fail(
            "Expected two separate JSON outputs, but got something else.")
    summary_output = json.loads(outputs[0])
    coverage_output = json.loads(outputs[1])
    expected_coverage = {
        "Coverage": {
            "MGYG000004832_79|kraken:taxid|6050": {
                "covered_bases_unique": 9,
                "covered_bases_ambiguous": 0,
                "mean_coverage_unique": 0.8,
                "mean_coverage_ambiguous": 0.0
            }
        }
    }
    expected_summary = {
        "Statistics": {
            "unique_mapped_reads": 3, "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0},
        "Summary":
            {"MGYG000004832_79|kraken:taxid|6050":
                 {"unique_reads": 2, "ambiguous_reads": 0},
             "MGYG000004832_78|kraken:taxid|1111":
                 {"unique_reads": 1, "ambiguous_reads": 0}}}
    assert coverage_output == expected_coverage
    assert expected_summary == summary_output
    os.remove(genome_path)
    os.remove(reads_path)
    os.remove(kdb_file_path)


def test_reverse_extension():
    """This function test the reverse extension on func map_reads. to do so func
    creates a FASTA/Q files to make an alignment object and then activate the
    map_reads() function. Func checks the mapping process resulted the expected
    mapping and that the chosen orientation for the Uniquely mapped reads are
    the reverse form"""
    ref_content = """>MGYG000004169_21|test1:taxid|51124
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
>MGYG000004169_22|test2:taxid|51315
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGA
>MGYG000004169_23|test3:taxid|51316
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCA
GCTCCGGCTGACAGCGCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGAAAAAA
ACTGGACACGGCCCA
>MGYG000004169_33|test4:taxid|51316
CCAGCGTTATGAATGGCAAATGGAACGGCAATAGTAAATCCAGCCAGTAGAACTATAGCA
TGTACAATATCAGTC
>MGYG000004169_33|test5:taxid|51316
CGGACATAAAGTTTCCAATGGTATCCGGATCAGTACTCAGTAGTTCCATCAAGACCTGAT
AACCCTCTCCGTCTC
"""
    reads_content = """@MGYG000002096_86|kraken:taxid|3024_read_0
GCCCATTGTCAAGAGCATTACGGCGGCGGCCTATCCCGGTTTTGAACTGGACACGGCCCAGCTCCGGCTGACAGC
+
DC@@B>>A>BBA??B=@CB@==?A><@??=B>><==A??>=@><;;>?>=@;?>?99>=>>:;>:8>8;87;97:
@MGYG000004832_79|kraken:taxid|6050_read_1
GCCCATTGTCATGAGCATTATGGCGGCGGCCTATCCAGGTTTTGAACTGGACACGTGCCAGCTCCGGCTGACAGT
+
B@AC??BC>>?CA@B@?A@?==BB@<@A?=@>=A?<AAA?@A;@;;><?=<=;=;9;>:<?<98;<==987887<
@MGYG000002096_86|kraken:taxid|3024_read_2
ACCCGCATCCGGAATGCGGACGATGTGGAGAAGTTACTGGATCTGCCGGTTTTGGGGCAGATTCCGCGGGAAGAT
+
B@?CCD>C>BD>ACCB?A>@B>?A?@=BABB@@A>==;=>;<@==>==><<:??;9?;?=:9=<9<:<>;<;<87
@MGYG000000975_36|kraken:taxid|1603_read_3
CAGCACTATGACCATAAACGCTTCCATAAACAATTAATGCACCTTCAACTTCAGGTGTATATGAACTCCATTTAT
+
DDE@DB@@AA>BA?=?=B=>CA?CAB=??B<B@?><;;@;@@;;??@<;><?>?=>::<>?8<=;<=;9879:8;
@MGYG000001129_3|kraken:taxid|1782_read_4
CGCTGGAAGTGTCCATGACGGTACCAAGTGGCAGACGCGGTAGCTTGGCAAGTAATTTTTCTCGAATGGTTTTCA
+
@?BBC@A@>D@@BA@@=CBA@?CB?>>?=?<BA?A=@=>==><=;=<@;>?;:>>>9>?>9>=;8<;89<9::8:
@MGYG000004787_66|kraken:taxid|6005_read_5
GGGTCGCCCTCCACCGCCGCCAGCAAGGGCAGCACGTGCTCGTATTTGGTGCTCTTCATGCCGGAGGGCAGCAGC
+
@DDEEDDD@@A>@B@>>AB?C>=@=?@<A<?<<?;A@?<A?A?=@<?=:;<;=;9>>?9<=>8>;<:8>>=9;9<
@MGYG000004787_66|kraken:taxid|6005_read_6
ACGAAAATTCGTCACGCCGCTGGTGTTCCGGCTAAAATCAGCCCCGCTCACGCATCTTGGCGAAGCAGTCGTTCA
+
BD?CCBCD>>>ADB=CC@B>@>BAB>A><>=A@?<>==A<??>?<=@;;@=>9==>><9<98<:<>;88>=9::<
@MGYG000000125_5|kraken:taxid|343_read_7
AGACTCAGATTATTTCTGTACAGTTTGGACTACATCTGCGGCATCTTTTATATTATCATAGGCAGGAGCCAACAA
+
?DABEA@BADB?A>@==>CBCBA?@<=@B=>?>?A=?<;A<=@;><<;=:<?<=?9;?>:9:9:=<8=><;9:98
@MGYG000004169_19|kraken:taxid|5353_read_8
TCAGTACCTTGGCAGAAACTTTGACGAAGGTATCCGGAACAATCATGAAGTCTGACGGAGCAGCGGAGTTAGGCA
+
BBEBDC>D?>B?@D=B@=?A@CAA@@@AA>A=><?A<><;=?@;>@;@;<:9<<==;:;<9>98;9<;;8=;<77
@MGYG000002096_86|kraken:taxid|3024_read_9
GCGGATGGGGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGG
+
?EC@@?B?CB@D@CCBBCB?=?CAB?<>B?A>A@;@>A<??=?;<@;@<@;>?9:;<:<=98=9=;8;8>;9;=;
@MGYG000004169_19|kraken:taxid|5353_read_10
ATGAATTTTTGAGAAAAGCAGGTGAACAGGTTGACTGAAAATTGTAAGAAATATGCAGTAGAGCTGCATGAACTC
+
A@?CE>A?@?B>A>BA?>=>B>?C=>B@B=>??>><?<=<;@:<;?>:=:<?=>9:9:>;<>;=9>8<:=7;=:;
@MGYG000002096_86|kraken:taxid|3024_read_11
AGAAAATGGGCAACAATATAGGAAAGGATGAACACTGTTTGAAGAGAGGAGCACGCGGATGGGGAAAAATATGGA
+
?@E@@BBB>ABCB?CB????=>>?BBA?=@A?<?<A@@@;>A??;:;@?=@?9<999<9==>>9;==>8:888<;
@MGYG000003260_60|kraken:taxid|4334_read_12
GGAAAAATATGGACGGCGCTCATTGCGGCAGCGGCGCTTTTAGGGACGACGCTGAGCGCATCGGCGGCGACGAGG
+
DBA?CBCABCA@A?@>AB==AB>?AB>@A=<<B=A?;=;;=?;=;;>>;@<>:><>;<=?><<<9;>9>;:;=77
@MGYG000002096_86|kraken:taxid|3024_read_13
GAGGGACTGGAAATTATCATTGCCCATCCCGAGCGATATCAATATGTGCCCTGGGGATCGCTCTGGACTGGCAGG
+
ABC@D@ABBCDD>B=>@C=@@>CCABB>A@@?@@=<@A=;;>@;;@=??@;?9=>??=;;:9:8:<<=;8=8<9:
@MGYG000001009_274|kraken:taxid|1639_read_14
CCTGTGGCCGGGCTTGGCGTACGAGCGCCTGGACGGTCTGCGCGAGTGGTTCTTTGGCGACTTTGAGGCCGAGCG
+
DC?DCABCDD@DA@B@B??A?C@CAAAA@AB<A@==A?;@A><@;;>?<;::;=?;=><9?=>=>9;9<9==:78
@MGYG000004832_79|kraken:taxid|6050_read_15
TGGCCCTGATACAGGCCGGTCAAAAACTCCGCCTCGGTGAGGTTGGGCACGATCACGTCGGCCACGCCGATGAGA
+
CDCC@@AD>AC>C?=>BB@?CB@=B=AA@B@BA<?<@;?<A<:==<>;:??:<9=>9=9;:>>9;98>8=<;<;9
@MGYG000000125_5|kraken:taxid|343_read_16
TCCAAAAGTGCATCCATCAACCCATAACTACCCTCCCGACAAGATGTCACAAATAGCCTTTAAATGTAGCGTTTA
+
EA?@C??DD@>CDB>@=C@>B@@A==A?A<B<@;=><=;=;A?=>>@<?@<>=<=:::;9<9>>==<>:8<8=;:
@MGYG000002096_86|kraken:taxid|3024_read_17
TGGCCGGATTTTTGGGGCTTGTCATTAGCCGCATCCTTGGCCTGGGGCTGTCCTGGGCGGTGCTTTGCATTGCCT
+
DA?C@@AC?CD?BDACCAC@B@@>B<AA<A@AB@=?A@;?>;>;??<><>;:=>??<?:999<=<8>8:9<=<;;
@MGYG000004832_79|kraken:taxid|6050_read_18
GATCTGGGCCAGACGGCTGTCGTCCCGCCCGGCGGACAGGGTGTCCACCAGCTGGGCGAAGGGCCCCTGGGGGTC
+
ECB??AD@?CDBBBAB=A?BBC??=B<AA>>@=><A<<=>A=<@==;<;:@<>?;>:==??989:=:=8:8;9=:
@MGYG000000975_36|kraken:taxid|1603_read_19
CTAATCATAACACTATCTTCGATTTCACAAATAAATCCATAATGTGAACCTAAATCAATAACTTCAACAACATTT
+
@@BCC?D>>B>@BCC?=@A>CA?>@<<=<@AB?A@A=;>;<A??:?>@;=<>>>?;;<::;9:89>8:<<;;:;8
@MGYG000002096_86|kraken:taxid|3024_read_20
GCCCATTGTTAAGAGCATTATGGCGGCGGCCTATCCCGGTTTTGTACTGGACACGGCCCATCTCCGGCTGACTGC
+
EE?@?D>@A@@A?>CBB@C?@=@C?@==<A=>A<?=@A;A??@<;@?:?=>:99?<;:=:=;=8>:;9;;;8;<:
@MGYG000004169_19|kraken:taxid|5353_read_21
CGAGAAGGAAGACGGTACAGTGGTCAACATGGGGGATCAGGTGCTGGATAATCTTCATGAGAATACCAGTCTGGG
+
AD?@D@BADC@>BC=@@BB@A?@?==A@=<>?A>A=@?<<><;=:;??>=:?<?;:=?=9<::9:<>8=978=9;
@MGYG000004832_79|kraken:taxid|6050_read_22
ATGGTCAGCACATAGACCCGCCCCTCCTCCCGCTGGATGAGGGAGGTACCCAGGTCCACCAGGTTCTGCGTGTTG
+
DAD@CB?BC>BD?DC?>A=?CA=>?B<A=<?BA<<;=>>?;>>:=@;@:>=?>:>?>=9=>:>==:=9;;<87;<
@MGYG000000125_5|kraken:taxid|343_read_23
GTTTTTTGAGAAGAGATTCTCCCATTGCCGAAGATCCCATTCCTGCATCGCAGGCAAACACTACTAGATTAATAT
+
EE?ACC?CC@DBA>>=A==?CAAA>??>B=<=AA<;=@;>=?:=:@===:?>;<?;=:9>?=>:<=9:<;<9<:9
@MGYG000004787_66|kraken:taxid|6005_read_24
CAATGCCGAGGATTTCGCCGCCGTACGCCGAAAACCCGACCTTGTCAAGCAGCGTGACGCCGTCGTCGCTGATGC
+
?@CED>BBBD@?CCBAA@ABBCB?><BA@<>?=@;<<>;<A=@;<;:@;>=9;?=>::::9=>9<;9><999879
@MGYG000000975_36|kraken:taxid|1603_read_25
CTAGGAATATTAAATAATTTTAATAATTCTTCATCATAGTCAAGTGTGTTGATATTAAAAAGCATTGTTCTAGAA
+
E?B?A>AC@C@BBCCA>=?A=>>@AB<<>@<A><>A><=?AA:;=:@@>@=<<>:9=:?;>9:><:<><=7:8:=
@MGYG000000125_5|kraken:taxid|343_read_26
GACTGATATTGTACATGCTATAGTTCTACTGGCTGGATTTACTATTGCCGTTCCATTTGCCATTCATAACGCTGG
+
@DA@@@D@>>BDDBAC=CBCB>CA>==@?BB<>A<@@<;@A><=><;<<>;9?<=?=;=<==8>98>>:99=;7;
@MGYG000004169_19|kraken:taxid|5353_read_27
GAGACGGAGAGGGTTATCAGGTCTTGATGGAACTACTGAGTACTGATCCGGATACCATTGGAAACTTTATGTCCG
+
ADACBCACDDDCAB@CA>A??=?@??@?@BA?A@A<@@;A>A=>>:;:@:=9;><=><>99;<>>==99>;<7;=
@MGYG000002720_43|kraken:taxid|3759_read_28
AGTGCGGAGATCAGCAGGGAACCCACAAGGTTCATCGCCAGCACGATGATGACCGCAACGACGACGGCGATCAGC
+
E@DBABCB>B>C?@?>?AA?AC?C=>@@=?@=<?=<<<A>@;;>;<?;>>::;9<>?9???;<9;<8:=87:8;:
@MGYG000002720_43|kraken:taxid|3759_read_29
TTTCCGCGATCCTGTGCTGTGCGCCCTCGATGGTCAGCACGCAGGGAAGTCCCAGCTCATCCACCTTCTTAGCAA
+
B?D@?CDCABD@D?@C?=CBBA>C=@>BA@B???=??A@A?><>@?:?;><?>9:>=>;?;;9;=:<8:9==<=<"""

    reads_path = "test_reverse.fq"
    genome_path = "test_reverse.fa"
    k = 31
    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    with open(reads_path, "w") as file:
        for line in reads_content.splitlines():
            file.write(line.strip() + "\n")
    refs_dict = generate_ref_dict(genome_path)
    refs_kmer_db, genome_bases = build_ref_kmers_db(refs_dict, k)
    reads_dict = generate_reads_dict(reads_path)
    aln = Alignment()
    aln.set_references(list(refs_dict.keys()))
    aln.set_reads(reads_dict)
    aln.set_ref_kmers_db(refs_kmer_db)
    aln.build_reads_kmers_db(k, None)
    aln.map_reads(1, 1, None, True)
    results = {"Unmapped": 0, "Unique": 0, "Ambiguous": 0}
    for read_header, read in aln.get_reads().items():
        status = read.get_status()
        if status == "Unique":
            assert read.get_orient() == "Reverse"
        results[status] += 1
    assert results["Unique"] == 2
    assert results["Unmapped"] == 28
    os.remove(reads_path)
    os.remove(genome_path)


def test_determine_orientation():
    """This func test the logic of __determine_orientation func() in class
    alignment. func receive few cases of different Kmers and supposed to retunn
    True/ False according to logic"""
    aln = Alignment()
    specific = ['A', 'A', 'A']
    unspecific = ['A', 'A', 'A']
    rev_specific = ['A', 'A', 'A', 'A']
    rev_unspecific = ['A']
    status_f = "Unique"
    status_r = "Unique"
    # both orientations result Unique status - choose the one with more Unique Kmers
    assert aln._Alignment__determine_orientation(status_f, status_r, specific,
                                                 unspecific,
                                                 rev_specific,
                                                 rev_unspecific) == True
    status_f = "Ambiguous"
    status_r = "Ambiguous"
    assert aln._Alignment__determine_orientation(status_f, status_r, specific,
                                                 unspecific,
                                                 rev_specific,
                                                 rev_unspecific) == False
    status_f = "Unmapped"
    status_r = "Ambiguous"
    assert aln._Alignment__determine_orientation(status_f, status_r, specific,
                                                 unspecific,
                                                 rev_specific,
                                                 rev_unspecific) == True


def test_build_kmer_db_quality():
    """This func test the build_kmers_db() function in class alignment to make
    sure it filters out low quality Kmers while providing MKQ threshold"""
    # need to make sure low quality kmers are filtered out
    read = Read("read1", "AACCAA", "IIIIKK")
    reads_dict = {"reads1": read}
    aln = Alignment()
    aln.set_reads(reads_dict)
    filtered_out = aln.build_reads_kmers_db(3, 41)
    assert filtered_out == 3


def test_similarity_extension1(capsys):
    """This func test similarity extension on execute_dumpref_and_reference()
    to make sure similar genomes are filtered out and the expected similarity
    summary is printed out"""
    genome_path = "test_sim.fa"
    ref_content = """>genome1
AAAGGGTCCCTTCTAGTCGGCA
>genome2
TTGAACCTTCGGAACCGGCAATTTGGCC
>genome3
TTGAACCTTCGGAACTGGAATTGGCCTA"""

    with open(genome_path, "w") as file:
        for line in ref_content.splitlines():
            file.write(line.strip() + "\n")
    k_size = 3
    execute_dumpref_and_reference(genome_path, k_size, True, 0.37)
    expected_dumpref = {"Kmers":
                            {"AAA": {"genome1": [0]}, "AAG": {"genome1": [1]},
                             "AGG": {"genome1": [2]}, "GGG": {"genome1": [3]},
                             "GGT": {"genome1": [4]},
                             "GTC": {"genome1": [5, 15]},
                             "TCC": {"genome1": [6]}, "CCC": {"genome1": [7]},
                             "CCT": {"genome1": [8], "genome2": [5]},
                             "CTT": {"genome1": [9], "genome2": [6]},
                             "TTC": {"genome1": [10], "genome2": [7]},
                             "TCT": {"genome1": [11]},
                             "CTA": {"genome1": [12]},
                             "TAG": {"genome1": [13]},
                             "AGT": {"genome1": [14]},
                             "TCG": {"genome1": [16], "genome2": [8]},
                             "CGG": {"genome1": [17], "genome2": [9, 15]},
                             "GGC": {"genome1": [18], "genome2": [16, 24]},
                             "GCA": {"genome1": [19], "genome2": [17]},
                             "TTG": {"genome2": [0, 22]},
                             "TGA": {"genome2": [1]},
                             "GAA": {"genome2": [2, 11]},
                             "AAC": {"genome2": [3, 12]},
                             "ACC": {"genome2": [4, 13]},
                             "GGA": {"genome2": [10]},
                             "CCG": {"genome2": [14]},
                             "CAA": {"genome2": [18]},
                             "AAT": {"genome2": [19]},
                             "ATT": {"genome2": [20]},
                             "TTT": {"genome2": [21]},
                             "TGG": {"genome2": [23]},
                             "GCC": {"genome2": [25]}},
                        "Summary": {
                            "genome1": {"total_bases": 22, "unique_kmers": 13,
                                        "multi_mapping_kmers": 7},
                            "genome2": {"total_bases": 28, "unique_kmers": 17,
                                        "multi_mapping_kmers": 9}}}
    expected_sim = {"Similarity": {
        "genome3": {"kept": "no", "unique_kmers": 2,
                    "total_kmers": 26,
                    "genome_length": 28,
                    "similar_to": "genome2",
                    "similarity_score": 0.8421052631578947},
        "genome2": {"kept": "yes", "unique_kmers": 3,
                    "total_kmers": 26,
                    "genome_length": 28, "similar_to": "NA",
                    "similarity_score": "NA"},
        "genome1": {"kept": "yes", "unique_kmers": 12,
                    "total_kmers": 20,
                    "genome_length": 22, "similar_to": "NA",
                    "similarity_score": "NA"}}}
    output: str = capsys.readouterr().out.strip()
    outputs = output.split('\n', 1)
    if len(outputs) != 2:
        pytest.fail(
            "Expected two separate JSON outputs, but got something else.")
    summary_output = json.loads(outputs[0])
    similarity_output = json.loads(outputs[1])
    assert summary_output == expected_dumpref
    assert similarity_output == expected_sim
    os.remove(genome_path)
