# editcall
Call indels and inversions for CRISPR NGS data processed with bwa-mem

This tool can call big deletions and inversions from the split mapping (when there is `SA:Z:`tag) in the sam files from `bwa mem`.

Some notes about `bwa mem` sam file:

1. The primary alignment has soft clipping (S in the CIGAR), and the Supplementary alignment has hard clipping (H in the CIGAR).
2. Unmapped reads have CIGAR "*"

## Compile for local use
You need to install [Nim](https://nim-lang.org/) compiler first. Then compile it with command:  
`nim c -d:release nimBoxshadeMini.nim`

Then check the usage by typing:  
`./editcall`

## Compile for webassembly
I get the instruction on how to compile a Nim program to Webassembly [here](https://github.com/treeform/nim_emscripten_tutorial). You need to install [Emscripten](https://emscripten.org/docs/getting_started/downloads.html) first. Then use the command below to compile it:  
`nim c -d:emscripten -o:editcall.js editcall.nim`

Then you can use them as the way in [biowasm](https://github.com/biowasm/biowasm). You can pass more parameters to `emcc` in the last line of file `config.nims`.

For some reason, this time "std/parseopt" did not work in the webassembly, so I have to use positioned based command line arguments. I also found it was not working to run the program multiple times in the browser. So I made it to run all the sam files in the same time (give all the sam file names to the command line).

You can comment out the line `--define:release` in the file `config.nims` to debug it with `node`:

```
node editcall.js
```

### Compile to a satic library for webassembly
```sh
# 1. compile it with command:  
nim c --app:staticlib -d:emscripten editcall.nim
# 2. edit the "./tmp2/@meditcall.nim.c" and change the `main` function to some other names, for example, `main_editcall`.
# 3. then manually compile it again:  
emcc -c  -w -ferror-limit=3 -O3   -I/usr/lib/nim -I. -o tmp2/@meditcall.nim.c.o tmp2/@meditcall.nim.c

emar rcs libeditcall.a  tmp2/stdlib_assertions.nim.c.o tmp2/stdlib_dollars.nim.c.o tmp2/stdlib_formatfloat.nim.c.o tmp2/stdlib_io.nim.c.o tmp2/stdlib_system.nim.c.o tmp2/stdlib_hashes.nim.c.o tmp2/stdlib_math.nim.c.o tmp2/stdlib_tables.nim.c.o tmp2/stdlib_parseutils.nim.c.o tmp2/stdlib_strutils.nim.c.o tmp2/stdlib_sequtils.nim.c.o tmp2/stdlib_times.nim.c.o tmp2/stdlib_os.nim.c.o tmp2/@meditcall.nim.c.o
```

## Usage

```
Usage: ./editcall reference.fa output.txt samfile1.sam [...more samfiles]
```
You need to provide the reference fasta file that is used for `bwa mem`, the output file name, and all the sam files you want to check.

Right now it only support sam files. If you have bam files, please run  
`samtools view your.bam > your.sam`  
to convert the bam to sam file.