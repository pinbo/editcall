## modules
import std/tables
import strutils
import std/parseopt
import os
import std/re # regex, has split

proc writeHelp() =
  quit("synopsis: " & getAppFilename() & " -i=input.sam -g=reference.fa -o=output.txt")
proc writeVersion() =
  quit(getAppFilename() & " version 0.1.0")

echo paramCount(), " parameters" #, " ", paramStr(1)
if paramCount() < 1:
  writeHelp()


var 
  samFile = ""
  refFile = ""
  outfile = "out.txt"
# var outwidth = 60
# var conflag = false # print consensus line?
var p = initOptParser( commandLineParams() )
for kind, key, val in p.getopt():
  case kind
  of cmdArgument:
    samFile = key
  of cmdLongOption, cmdShortOption:
    case key
    of "help", "h": writeHelp()
    of "version", "v": writeVersion()
    of "i": samFile = val
    of "g": refFile = val
    of "o": outfile = val
    # of "w": outwidth = parseInt(val)
    # of "consensus", "c": conflag = true
  of cmdEnd: assert(false) # cannot happen
if samFile == "":
  # no samFile has been given, so we show the help
  writeHelp()


## read fasta file

# read a fasta file
proc readFasta(infile: string): OrderedTable[string, string] =
  # let contents = readFile(infile).strip() 
  # echo contents
  # let lines = contents.splitLines()
  var seqDict = initOrderedTable[string, string]()
  var seqName = ""
  for ll in lines infile:
    if contains(ll, ">"):
      seqName = ll.replace(">", "")
      seqDict[seqName] = ""
    else:
      seqDict[seqName].add(ll)
  return seqDict
## read fasta
# var testFasta = {"seq1": "MRDRTHELRQGDN", "seq2": "MKDRLEQLKAKQL", "seq3":"MRDRLPDLTACR-"}.toTable  # creates a Table
var seqDict = readFasta(refFile)

# init dict for variations
var indelDict = initTable[string, int]()
# var readDict = initTable[string, int]() # key is readID-strand

# process cigar and return
# all positions are 0-based
proc parseCigar (cigar: string, refPos: int, sameStrand: bool, readLen: int): seq[int] =
  let ss1 = cigar.split(re"\d+") #@["", "M", "D", "M", "S"]
  let ss2 = cigar.split(re"\D") #@["60", "5", "56", "26", ""]
  # all positions are 1-based
  var 
    readPos1 = 0 # left start
    readPos2 = -1 # right end
    refPos1 = refPos # left start
    refPos2 = refPos - 1 # right end
    nM = 0 # number of M, if match showed up, then no more S or H
  for i in 1 .. ss1.len-1:
    let num = ss2[i-1].parseInt
    if ss1[i] == "M" or ss1[i] == "=" or ss1[i] == "X":
      readPos2 += num
      refPos2 += num
      nM += 1
    elif ss1[i] == "S" or ss1[i] == "H":
      if nM == 0:
        readPos1 += num
        readPos2 += num
    elif ss1[i] == "I":
      readPos2 += num
    elif ss1[i] == "D" or ss1[i] == "N":
      refPos2 += num - 1
  if sameStrand:
    return @[readPos1, readPos2, refPos1, refPos2]
  else:
    return @[readLen - readPos2 - 1, readLen - readPos1 - 1, refPos1, refPos2]

## read files line by line
for line in lines samFile:
  if line[0] == '@':
    continue
  else:
    let ff = line.split("\t")
    let 
      readID = ff[0]
      flag   = ff[1].parseInt
      chrom  = ff[2]
      pos  = ff[3].parseInt - 1 # 0 based
      cigar  = ff[5] 
      readSeq= ff[9]

    let strand = if (flag and 0x10) != 0: "-" else: "+"
    echo readID
    # if readID & "-" & strand in readDict:
    #   continue
    # else:
    #   readDict[readID & "-" & strand] = 1
    var refSeq = ""
    var altSeq = ""
    var kk = "" # key of indel dict
    var readPos = 0 # relative positions for clipped reads
    var rawReadPos = 0 # positions of the raw read before clipping
    var refPos = pos # 0-based
    # check if there is insertion or deletion
    if ('I' in cigar or 'D' in cigar or 'N' in cigar) or ("SA:Z" in line and "H" notin cigar):
      let ss1 = cigar.split(re"\d+") #@["", "M", "D", "M", "S"]
      let ss2 = cigar.split(re"\D") #@["60", "5", "56", "26", ""]
      for i in 1 .. ss1.len-1:
        let num = ss2[i-1].parseInt
        if ss1[i] == "M" or ss1[i] == "=" or ss1[i] == "X":
          readPos += num
          refPos += num
          rawReadPos += num
        elif ss1[i] == "S":
          readPos += num
          rawReadPos += num
        elif ss1[i] == "H":
          rawReadPos += num
        elif ss1[i] == "I":
          refSeq = seqDict[chrom][refPos-1 .. refPos]
          altSeq = readSeq[readPos-1 .. readPos+num]
          kk = @[chrom, $refPos, refSeq, altSeq].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
          readPos += num
          rawReadPos += num
        elif ss1[i] == "D" or ss1[i] == "N":
          refSeq = seqDict[chrom][refPos-1 .. refPos+num]
          altSeq = readSeq[readPos-1 .. readPos]
          kk = @[chrom, $refPos, refSeq, altSeq].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
          refPos += num
    # check if having supplementary read
    # SA:Z:6B,185,-,54S62M26S,10,0;
    if "SA:Z" in line and "H" notin cigar: # only when the read has no hard clip, eg. read is raw
      let ff2 = line.split("SA:Z:")[1].split(",")
      let 
        saChrom  = ff2[0]
        saPos  = ff2[1].parseInt - 1 # 0-based this is close to the border on the left, may need to adjust
        saStrand = ff2[2]
        saCigar  = ff2[3]
      # var readPos2 = 0 # read position for SA
      let ss1 = cigar.split(re"\d+") #@["", "M", "D", "M", "S"]
      # let ss2 = cigar.split(re"\D") #@["60", "5", "56", "26", ""]
      let ss3 = saCigar.split(re"\d+") #@["", "M", "D", "M", "S"]
      let ss4 = saCigar.split(re"\D") #@["60", "5", "56", "26", ""]
      var rawReadPos2 = 0 # raw read position on SA
      if chrom == saChrom and strand == saStrand:# and (flag and 0x20) == 0: # big deletion, could be insertion too, but update later
        echo "potential big deletion"
        var allPos1: seq[int]
        var allPos2: seq[int]
        if (saPos > pos): # SA is on the right
          allPos2 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
          allPos1 = parseCigar(cigar, pos, true, readSeq.len)
        else:
          allPos1 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
          allPos2 = parseCigar(cigar, pos, true, readSeq.len)
          # for i in 1 .. ss3.len-1:
            # let num2 = ss4[i-1].parseInt
            # if ss3[i] == "S" or ss3[i] == "H":
            #   rawreadPos2 += num2
            # elif ss1[i] == "M" or ss1[i] == "=" or ss1[i] == "X": # first match, should be close to the deletion border
            #   let shift = if rawReadPos > rawReadPos2: rawReadPos - rawReadPos2 else: 0
            #   let delEndPos = saPos + shift
            #   echo refPos-1
            #   echo delEndPos-1
            #   refSeq = seqDict[chrom][refPos-1 .. delEndPos-1]
            #   altSeq = readSeq[rawReadPos-1 .. rawReadPos2+shift]
            #   kk = @[chrom, $refPos, refSeq, altSeq].join("\t")
            #   if kk in indelDict:
            #     indelDict[kk] += 1
            #   else:
            #     indelDict[kk] = 1
          # let allPos2 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
          # let allPos1 = parseCigar(cigar, pos, true, readSeq.len)
        echo allPos1
        echo allPos2
        let
          readPos1 = allPos1[1]
          refPos1  = allPos1[3]
          readPos2 = allPos2[0]
          refPos2  = allPos2[2]
          shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
          delEndPos = refPos2 + shift
        echo refPos1
        echo shift
        echo refPos2
        echo delEndPos
        refSeq = seqDict[chrom][refPos1 .. delEndPos-1]
        altSeq = readSeq[readPos1 .. readPos2+shift]
        kk = @[chrom, $(refPos1+1), refSeq, altSeq].join("\t")
        if kk in indelDict:
          indelDict[kk] += 1
        else:
          indelDict[kk] = 1
        # else: # SA in on the left
        #   let allPos1 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
        #   let allPos2 = parseCigar(cigar, pos, true, readSeq.len)
        #   # echo allPos1
        #   # echo allPos2
        #   let
        #     readPos1 = allPos1[1]
        #     refPos1  = allPos1[3]
        #     readPos2 = allPos2[0]
        #     refPos2  = allPos2[2]
        #     shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
        #     delEndPos = refPos2 + shift
        #   echo refPos1
        #   echo shift
        #   echo refPos2
        #   echo delEndPos
        #   refSeq = seqDict[chrom][refPos1 .. delEndPos-1]
        #   altSeq = readSeq[readPos1 .. readPos2+shift]
        #   kk = @[chrom, $(refPos1+1), refSeq, altSeq].join("\t")
        #   if kk in indelDict:
        #     indelDict[kk] += 1
        #   else:
        #     indelDict[kk] = 1
      elif chrom == saChrom and strand != saStrand: # inversions
        echo "Potential inversion"
        var allPos1: seq[int]
        var allPos2: seq[int]
        if (saPos > pos): # SA is on the right
          allPos1 = parseCigar(cigar, pos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
          allPos2 = parseCigar(saCigar, saPos, false, readSeq.len)
        else:
          allPos1 = parseCigar(saCigar, saPos, false, readSeq.len)
          allPos2 = parseCigar(cigar, pos, true, readSeq.len)
        echo allPos1
        echo allPos2
        var
          readPos1 = allPos1[1]
          refPos1  = allPos1[3]
          readPos2 = allPos2[1]
          refPos2  = allPos2[3]
          # shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
          # delEndPos = refPos2 - shift
        if (flag and 0x80) != 0: # 2nd in pair, need to count from right
          readPos1 = allPos2[0]
          refPos1  = allPos1[2]
          readPos2 = allPos1[0]
          refPos2  = allPos2[2]
        let shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
        let delEndPos = refPos2 - shift
        echo refPos1
        echo refPos2
        echo readPos1
        echo readPos2
        echo delEndPos
        refSeq = seqDict[chrom][refPos1 .. delEndPos-1]
        altSeq = "inversion"
        kk = @[chrom, $(refPos1+1), refSeq, altSeq].join("\t")
        if kk in indelDict:
          indelDict[kk] += 1
        else:
          indelDict[kk] = 1










## print at the end
for k, v in indelDict:
  echo k & "\t" & $v


# bitwise check example
# proc bitwise(a: int, b: int) =
#   echo "a and b: " , a and b
#   echo "a or b: ", a or b
#   echo "a xor b: ", a xor b
#   echo "not a: ", not a
#   echo "a << b: ", a shl b
#   echo "a >> b: ", a shr b

# bitwise(145, 0x40)

# split cigar
# import std/re # regex, has split
# let cigar = "60M5D56M26S"
# let ss1 = cigar.split(re"\d+")
# echo ss1
# let ss2 = cigar.split(re"[MIDNSHP=X]")
# echo ss2

# bitwise check example
# proc flagcheck(a: int) =
#   echo "   0x1     1  PAIRED         " , a and 0x1
#   echo "   0x2     2  PROPER_PAIR    " , a and 0x2
#   echo "   0x4     4  UNMAP          " , a and 0x4
#   echo "   0x8     8  MUNMAP         " , a and 0x8
#   echo "  0x10    16  REVERSE        " , a and 0x10
#   echo "  0x20    32  MREVERSE       " , a and 0x20
#   echo "  0x40    64  READ1          " , a and 0x40
#   echo "  0x80   128  READ2          " , a and 0x80
#   echo " 0x100   256  SECONDARY      " , a and 0x100
#   echo " 0x200   512  QCFAIL         " , a and 0x200
#   echo " 0x400  1024  DUP            " , a and 0x400
#   echo " 0x800  2048  SUPPLEMENTARY  " , a and 0x800


#    0x1     1  PAIRED         paired-end / multiple-segment sequencing technology
#    0x2     2  PROPER_PAIR    each segment properly aligned according to aligner
#    0x4     4  UNMAP          segment unmapped
#    0x8     8  MUNMAP         next segment in the template unmapped
#   0x10    16  REVERSE        SEQ is reverse complemented
#   0x20    32  MREVERSE       SEQ of next segment in template is rev.complemented
#   0x40    64  READ1          the first segment in the template
#   0x80   128  READ2          the last segment in the template
#  0x100   256  SECONDARY      secondary alignment
#  0x200   512  QCFAIL         not passing quality controls or other filters
#  0x400  1024  DUP            PCR or optical duplicate
#  0x800  2048  SUPPLEMENTARY  supplementary alignment
