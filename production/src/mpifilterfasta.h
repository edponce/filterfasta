#ifndef FILTERFASTA_H
#define FILTERFASTA_H

/*
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE
#define _XOPEN_SOURCE 600
#define _POSIX_C_SOURCE 200112L
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>
#include "utilities.h"

// Set default options
#define OUTPUT_FILE "filter.out"  // Default output file name
#define SEQ_COUNT   LLONG_MAX  // Max number of sequences to extract
#define ANNOT_CNT   INT_MAX    // default = ALL, -# = write first # annotation fields without the sequence, 0 = NONE,  # = write first # annotation fields with sequence
#define BYTES_LIMIT LLONG_MAX  // Max number of bytes to extract
#define PIPE_MODE   0          // 0 = NONE, 1 = HMMER, 2 = MUSCLE
#define SEARCH_MODE 0          // 0 = NONE, 1 = ENABLE 
#define VERBOSE_OPT 0          // 0 = OFF, 1 = ON
#define TRACE_OPT   0          // 0 = OFF, 1 = ON

// Set internal configurations (Do not change)
#define MAXARG_CNT  5	       // Max number of (range) sequence length options
#define FILE_LEN    128        // FILENAME_MAX, Max length for filenames
#define ERROR       -1         // Error code from failed functions
#define CFGERROR    -2         // Error code from invalid configuration

#define IMAP_LIMIT  (1LL<<28)  // Memory map chunk limit for query file, 256MB
#define STRM_BUFSIZ (1LL<<22)  // Size of output stream buffer, 4MB
#define BCAST_LIMIT (1LL<<22)  // Size for broadcasting files, 4MB
#define HITS_ID_LEN 64LL       // Max length for BLAST table query and hit IDs
#define VERBOSE(ctx) if(verbose||trace) {ctx} // Verbose mode
#define TRACE(ctx)   if(trace) {ctx}   // Trace mode (debug)  
#define MIN(a,b)     ((a < b) ? a : b)
#define MAX(a,b)     ((a > b) ? a : b)
//#define STDIN       "redirStdin.txt"
//#define STDOUT      "redirStdout.txt"
//#define STDERR      "redirStderr.txt"

// Structure for command line arguments
typedef struct st_args
{
  char	         qf[FILE_LEN];           // Query file
  char	         of[FILE_LEN];           // Output file
  char	         sf[FILE_LEN];           // Search file to extract user defined sequences
  char           btable[FILE_LEN];       // BLAST table file, used to extract hit IDs
  long long int  rseqLen[MAXARG_CNT*2];  // Range sequence length to extract
  long long int  seqLen[MAXARG_CNT];     // Sequence length to search
  long long int  seqCnt;                 // Max number of sequences to extract
  long long int  bytesLimit;             // Max number of bytes to extract
  int            seqLenBuf;              // Number of sequence length options
  int            rseqLenBuf;             // Number of range sequence length options
  int            annotCnt;               // Number of annotation fields to extract
  int            pipeMode;               // Pipeline program after extracting sequences 
  int            searchMode;             // Flag for search file sequence extraction 
} args_t;

// Structure for managing I/O and memory map
typedef struct st_iomap
{
  long long int  xCnt;        // Total number of extracted sequences
  long long int  qfsz;        // Size of query file
  FILE          *qfd;         // File descriptor of query file
  FILE          *ofd;         // File descriptor of output file
  long long int *fileOffs;    // File offsets for query file memory mappings
  char          *iMap;        // Pointer to initial mapped memory
  char          *fMap;	      // Pointer to last mapped memory
} iomap_t;

// Structure for managing queries
typedef struct st_query
{
  long long int  buflen;      // Length of intermediate buffer
  char          *buf;         // Intermediate buffer to store sequences between memory map partitions
  char          *fbuf;        // Pointer to end of buffer
  char          *iaq;         // Pointer to start of current annotation
  char          *faq;         // Pointer to end of current annotation
  char          *isq;         // Pointer to start of current sequence
  char          *fsq;         // Pointer to end of current sequence
} query_t;

// Structure for BLAST table query and hit IDs
typedef struct st_hits
{
  long long int  total;        // Total number of lines in BLAST table file
  long long int  qtotal;       // Number of distinct query IDs in BLAST table file
  long long int  htotal;       // Number of distinct hit IDs in BLAST table file
  long long int *idxList;      // Array to index hit IDs corresponding to query IDs (MUSCLE pipeline)
  FILE          *tfd;          // File descriptor of BLAST table file 
  FILE          *ofd;          // File descriptor of output file for sequences not found
  int            pipeMode;     // Pipeline program after extracting sequences 
  int            searchMode;   // Flag for search file sequnce extraction
  int           *charVect;     // Characteristic vector used to determine sequences found or not 
  char          *iMap;         // Pointer to initial mapped memory
  char          *fMap;	       // Pointer to last mapped memory
  char         **queryList;    // List of queries in BLAST table file
  char         **hitList;      // List of hit IDs in BLAST table file
} hits_t;


typedef struct st_mpi
{
  int       procCnt;	      // Total number of processes in world
  int       procRank;         // Rank index of current process
  int       nameLen;          // Length of processor name
  char      procName[MPI_MAX_PROCESSOR_NAME];     // Processor name
  MPI_Comm  MPI_MY_WORLD;     // MPI communicator object
} mpi_t;


int displayHelp();
int parseCmdline(int, char **, args_t *, mpi_t *);
int getAnnot(iomap_t *, query_t *);
int getSequence(long long int *, iomap_t *, query_t *);
int initQueryMap(long long int, long long int, iomap_t *, mpi_t *);
int openQueryFile(char *, iomap_t *);
int parseAnnot(int, long long int *, query_t *);
int extractQueries(args_t *, iomap_t *, query_t *, hits_t *, mpi_t *, long long int *, int *);
int adjustMapBegin(long long int *, iomap_t *, query_t *);
int adjustMapEnd(iomap_t *, query_t *);
int combineOutputFiles(args_t *, iomap_t *, mpi_t *, long long int);
int writeHitsNotFound(char *, hits_t *, mpi_t *);
int partQueryFile(args_t *, iomap_t *, hits_t *, mpi_t *);
int parseBlastTableIDs(hits_t *, char *);
int freeHitsMemory(hits_t *);
int loadSearchIDs(char *, hits_t *);
int loadBlastTable(char *, hits_t *);
int adjustMPIProcs(mpi_t *, int);
int getInputFilesComm(mpi_t *, MPI_Comm *);
int distributeInputFiles(args_t *, mpi_t *);
int setOffs(iomap_t *, mpi_t *);

#endif
