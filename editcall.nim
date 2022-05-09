## version 1: inversion direction was checked by Read1 vs Read2, which is only correct for my fastqs for CRISPR NGS
## version 2: inversion direction was checked by read segment
## version 3: add depth and percent of mutation in all reads
## todo: handle big insertions and SNPs

## modules
import std/tables
import strutils
import sequtils
import os
# import std/parseopt

proc writeHelp() =
  quit("editcall v0.1.0\nUsage: ./editcall reference.fa output.txt samfile1.sam samfile2.sam [...more samfiles]")

let inputFiles = commandLineParams()
echo inputFiles
# echo paramCount(), " parameters"
if paramCount() < 3:
  writeHelp()

let
  refFile = inputFiles[0] # reference fasta file
  outfile = inputFiles[1] # output file name

# read a fasta file
proc readFasta(infile: string): Table[string, string] =
  # let contents = readFile(infile).strip() 
  # echo contents
  # let lines = contents.splitLines()
  var seqDict = initTable[string, string]()
  var seqName = ""
  for ll in lines infile:
    if contains(ll, ">"):
      seqName = ll.replace(">", "")
      seqDict[seqName] = ""
    else:
      seqDict[seqName].add(ll)
  return seqDict

## split cigar
proc mysplit (cigar: string): (seq[char], seq[int]) =
  # echo cigar # 8M1D108M
  let num = "0123456789"
  var cc: seq[char]
  for i in cigar:
    if i notin num:
      cc.add(i)
  # echo cc
  var nn: seq[int]
  var tmp = ""
  for i in cigar:
    if i in num:
      tmp = tmp & i
    else:
      nn.add(tmp.parseInt)
      tmp = ""
  # echo nn
  return (cc, nn)


# process cigar and return
# all positions are 0-based
proc parseCigar (cigar: string, refPos: int, sameStrand: bool, readLen: int): seq[int] =
  let (ss1, ss2) = cigar.mysplit #@["M", "D", "M", "S"], @["60", "5", "56", "26"]
  var 
    readPos1 = 0 # left start
    readPos2 = -1 # right end
    refPos1 = refPos # left start
    refPos2 = refPos - 1 # right end
    nM = 0 # number of M, if match showed up, then no more S or H
  for i in 0 .. ss1.len - 1:
    let num = ss2[i]
    if ss1[i] == 'M' or ss1[i] == '=' or ss1[i] == 'X':
      readPos2 += num
      refPos2 += num
      nM += 1
    elif ss1[i] == 'S' or ss1[i] == 'H':
      if nM == 0:
        readPos1 += num
        readPos2 += num
    elif ss1[i] == 'I':
      readPos2 += num
    elif ss1[i] == 'D' or ss1[i] == 'N':
      refPos2 += num - 1
  if sameStrand:
    return @[readPos1, readPos2, refPos1, refPos2]
  else:
    return @[readLen - readPos2 - 1, readLen - readPos1 - 1, refPos1, refPos2]

## read fasta
# echo "start reading fasta"
let seqDict = readFasta(refFile)
# echo "fasta file has been read!"
let f = open(outfile, fmWrite)
# echo "outfile is opened for writing!"
# Sample	Gene	POS	REF	ALT	TotalCoverage	mutCoverage	indelPercent	indelSize
f.writeLine("samFile\tChrom\trefStartPos\trefEndPos\tRef\tAlt\tindelSize\tTotalCoverage\tmutCoverage\tindelPercent")
## process all files
for j in 2 .. inputFiles.len-1:
  let samFile = inputFiles[j]
  echo samFile
  var indelDict = initTable[string, int]() # init dict for variations
  var depthDict = initTable[string, seq[int]]() # init dict for depth for all positions
  for k, v in seqDict:
    depthDict[k] = repeat(0, v.len)
  ## read files line by line
  for line in lines samFile:
    if line[0] == '@':
      continue
    else:
      let ff = line.split("\t")
      let
        flag   = ff[1].parseInt
        chrom  = ff[2]
        pos  = ff[3].parseInt - 1 # 0 based
        cigar  = ff[5] 
        readSeq= ff[9]
      let strand = if (flag and 0x10) != 0: "-" else: "+"
      let readID = ff[0]
      if cigar == "*": # no mapping
        continue
      # echo readID
      var refSeq = ""
      var altSeq = ""
      var mutSize = 0
      var kk = "" # key of indel dict
      var readPos = 0 # relative positions for clipped reads
      var rawReadPos = 0 # positions of the raw read before clipping
      var refPos = pos # 0-based
      var depStart = pos # 0-based, define the region that is mapped
      var depRange: seq[int] # all positions with coverage
      # check if there is insertion or deletion
      # if ('I' in cigar or 'D' in cigar or 'N' in cigar) or ("SA:Z" in line and 'H' notin cigar):
      let (ss1, ss2) = cigar.mysplit #@["M", "D", "M", "S"], @[60, 5, 56, 26]
      for i in 0 .. ss1.len-1:
        let num = ss2[i]
        if ss1[i] == 'M' or ss1[i] == '=' or ss1[i] == 'X':
          readPos += num
          refPos += num
          rawReadPos += num
          let tmpRange = to_seq(depStart .. depStart+num-1)
          depRange.add(tmpRange)
          depStart += num
        elif ss1[i] == 'S':
          readPos += num
          rawReadPos += num
        elif ss1[i] == 'H':
          rawReadPos += num
        elif ss1[i] == 'I':
          refSeq = seqDict[chrom][refPos-1 .. refPos]
          altSeq = readSeq[readPos-1 .. readPos+num]
          mutSize = altSeq.len - refSeq.len
          kk = @[samFile, chrom, $refPos, $(refPos+1), refSeq, altSeq, $mutSize].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
          readPos += num
          rawReadPos += num
        elif ss1[i] == 'D' or ss1[i] == 'N':
          refSeq = seqDict[chrom][refPos-1 .. refPos+num]
          altSeq = readSeq[readPos-1 .. readPos]
          mutSize = altSeq.len - refSeq.len
          kk = @[samFile, chrom, $refPos, $(refPos+num+1), refSeq, altSeq, $mutSize].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
          refPos += num
          depStart += num
      # add coverage
      for tmpPos in depRange:
        depthDict[chrom][tmpPos] += 1
      # check if having supplementary read
      # SA:Z:6B,185,-,54S62M26S,10,0;
      if "SA:Z" in line and 'H' notin cigar: # only when the read has no hard clip, eg. read is raw
        let ff2 = line.split("SA:Z:")[1].split(",")
        let 
          saChrom  = ff2[0]
          saPos  = ff2[1].parseInt - 1 # 0-based this is close to the border on the left, may need to adjust
          saStrand = ff2[2]
          saCigar  = ff2[3]
        if chrom == saChrom and strand == saStrand:# and (flag and 0x20) == 0: # big deletion, could be insertion too, but update later
          # echo "potential big deletion"
          var allPos1: seq[int]
          var allPos2: seq[int]
          if (saPos > pos): # SA is on the right
            allPos2 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
            allPos1 = parseCigar(cigar, pos, true, readSeq.len)
          else:
            allPos1 = parseCigar(saCigar, saPos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
            allPos2 = parseCigar(cigar, pos, true, readSeq.len)
          # echo "potential big deletion"
          # echo allPos1
          # echo allPos2
          let
            readPos1 = allPos1[1]
            refPos1  = allPos1[3]
            readPos2 = allPos2[0]
            refPos2  = allPos2[2]
            shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
            delEndPos = refPos2 + shift
          if refPos1 > delEndPos:
            echo "refPos1 > delEndPos-1 for ", readID
            echo "refPos1 ", refPos1
            echo "refPos2 ", refPos2
            echo "readPos1 ", readPos1
            echo "readPos2 ", readPos2
            echo "delEndPos ", delEndPos
            echo "allPos1, ", allPos1
            echo "allPos2, ", allPos2
            continue
          refSeq = seqDict[chrom][refPos1 .. delEndPos]
          altSeq = readSeq[readPos1 .. readPos2+shift]
          mutSize = altSeq.len - refSeq.len
          kk = @[samFile, chrom, $(refPos1+1), $(delEndPos+1), refSeq, altSeq, $mutSize].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
        elif chrom == saChrom and strand != saStrand:         ### inversions
          var allPos1: seq[int]
          var allPos2: seq[int]
          if (saPos > pos): # SA is on the right
            allPos1 = parseCigar(cigar, pos, true, readSeq.len) # return readPos1, readPos2, refPos1, refPos2
            allPos2 = parseCigar(saCigar, saPos, false, readSeq.len)
          else:
            allPos1 = parseCigar(saCigar, saPos, true, readSeq.len)
            allPos2 = parseCigar(cigar, pos, false, readSeq.len)
          # echo "Potential inversion"
          # echo allPos1
          # echo allPos2
          var
            readPos1 = allPos1[1]
            readPos2 = allPos2[0]
            refPos1  = allPos1[3]
            refPos2  = allPos2[3]
            shift = if readPos1 >= readPos2: readPos1 - readPos2 + 1 else: 0
            delEndPos = refPos2 - shift
          if allPos1[0] > allPos2[0]:# count from right
            readPos1 = allPos1[0]
            readPos2 = allPos2[1]
            refPos1  = allPos1[2]
            refPos2  = allPos2[2]
            shift = if readPos2 >= readPos1: readPos2 - readPos1 + 1 else: 0
            delEndPos = refPos2 + shift
          if refPos1 > delEndPos:
            echo "refPos1 > delEndPos-1 for ", readID
            echo "refPos1 ", refPos1
            echo "refPos2 ", refPos2
            echo "readPos1 ", readPos1
            echo "readPos2 ", readPos2
            echo "delEndPos ", delEndPos
            echo "allPos1, ", allPos1
            echo "allPos2, ", allPos2
            continue
          refSeq = seqDict[chrom][refPos1 .. delEndPos]
          altSeq = "inversion"
          mutSize = refSeq.len
          kk = @[samFile, chrom, $(refPos1+1), $(delEndPos+1), refSeq, altSeq, $mutSize].join("\t")
          if kk in indelDict:
            indelDict[kk] += 1
          else:
            indelDict[kk] = 1
  # write to file
  # kk = @[samFile, chrom, $refPos, $(refPos+1), refSeq, altSeq, $mutSize].join("\t")
  for k, v in indelDict:
    let ss = k.split("\t")
    let
      chrom = ss[1]
      n1 = ss[2].parseInt - 1
      n2 = ss[3].parseInt - 1
      dep1 = depthDict[chrom][n1]
      dep2 = depthDict[chrom][n2]
    let cov = if dep1 > dep2: dep1 else: dep2
    let pct = v / cov * 100
    f.writeLine(k & "\t" & $cov & "\t" & $v & "\t" & $pct)

## print at the end
# let f = open(outfile, fmWrite)
# echo "outfile is opened for writing!"
# f.writeLine("samFile\tChrom\trefStartPos\trefEndPos\tRef\tAlt\tmutSize\tnMutRead")
# for k, v in indelDict:
#   f.writeLine(k & "\t" & $v)

f.close()
echo "Everything is done!"
# use a variable
# var fileContent = "Chrom\trefStartPos\trefEndPos\tRef\tAlt\tmutSize\tnMutRead\n"
# for k, v in indelDict:
#   fileContent.add(k & "\t" & $v & "\n")
# writeFile(outfile, fileContent)
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
