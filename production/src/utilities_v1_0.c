#include "utilities.h"

// Compute wall time of code region
double get_wtime()
{
  struct timeval tm;
  double wtm;

  if( gettimeofday(&tm, NULL) != 0 )
  {
    fprintf(stdout, "Warning: wall time function failed\n");
    return 0;
  }

  wtm = (double)tm.tv_sec + (double)tm.tv_usec * 0.000001;
  
  return wtm;
}


// Function that takes a file descriptor, its data size, and partition count to compute offsets for balanced partitions.
// Partitions are guaranteed to be independent based on a separator symbol supplied by the user.
// Function returns the partition count considered and an array where every 3 contigous elements correspond to the offset values for a particular partition considered as,
// [i]=file_offset from beginning of file, [i+1]=map offset from file offset, [i+2]=total map size
// Element 1 = page-based offset from initial file position
// Element 2 = offset from element 1 where the actual independent data begins
// Element 3 = size of independent data for this partition
// For programs that mmap before, use "element 1 + element 2" as your initial address for processing. For programs that will mmap partitions after, need to respect system's page size, mmap from "element 1" and use "element 1 + element 2" as your initial address for processing.
int computePartitionOffsets(long long int **offs, int *parts, int fd, long int sz, char sym)
{
  char *buffer;
  int err;
  int i;
  int foundSymbol;
  int lparts;
  long long int readOffs;
  long long int partSz; 
  long long int c;
  long long int j;
  long long int offset;
  long long int chunks;
  long int multiplier;
  long long int chunkOffs;
  long long int prevOffs;
  long long int bytesRead;

  // Invalid inputs
  if( *parts < 1 || sz < 1 )
  {
    *offs = NULL;
    fprintf(stdout, "Invalid values for partition offsets, check partition count or data size\n");
    return ERROR;
  }
 
  // Loop until valid offsets are computed
  lparts = *parts;
  while( 1 )
  {
    err = 0;
  
    // Allocate array offsets 
    *offs = realloc(*offs, (lparts * 3) * sizeof offs);

    // Single partition
    if( lparts == 1 )
    {
      (*offs)[0] = (*offs)[1] = 0;
      (*offs)[2] = (long long int)sz;
      break;
    }
  
    // Compute page size multiplier for partitions
    // Adjustment occurs if partitions are smaller than page size 
    chunks = (long long int)sysconf(_SC_PAGESIZE);
    multiplier = (long long int)round((double)sz / (lparts * chunks));
    if( multiplier < 1 )
    {
      lparts--;
      fprintf(stdout, "Warning: adjusted partitions to (%d) based on size\n", lparts);
      continue;
    }
    else
    {
      // Allocate buffer for adjusting partitions
      buffer = malloc(chunks * sizeof buffer);
     
      // Compute partition size aligned to page size 
      partSz = chunks * multiplier;
      
      // Find offsets for each partition
      for(i = 0; i < lparts && err == 0; i++)
      {
        // First and middle partitions
        if( i < (lparts - 1) ) 
        { 
          // First partition
          if( i == 0 )
          {
            (*offs)[i*3] = 0;
            (*offs)[i*3+1] = 0;
          }
          // Middle partitions
          else
          {
            // Compute a page size offset that is not after the end of the previous partition
            prevOffs = (*offs)[(i-1)*3] + (*offs)[(i-1)*3+1] + (*offs)[(i-1)*3+2];
            chunkOffs = (long long int)floor((double)prevOffs/chunks);
            (*offs)[i*3] = chunks * chunkOffs;
            (*offs)[i*3+1] = prevOffs - (*offs)[i*3];
          }

          // Loop backwards current partition in chunks of page size to find symbol
          j = 0;
          offset = 0;
          foundSymbol = 0;
          while( foundSymbol == 0 && err == 0 )
          {
            j++;

            // Compute offset from beginning of data to read at beginning of chunk
            readOffs = (*offs)[i*3] + partSz - (chunks * j);

            // Set file position at beginning of current chunk
            if( lseek(fd, readOffs, SEEK_SET) == -1 )
              fprintf(stdout, "Warning: failed lseek() in partition offsets\n");
      
            // Read a chunk of data, store in buffer
            // If data is mmap, we should read data directly
            bytesRead = read(fd, buffer, chunks);
            if( bytesRead != chunks )
              fprintf(stdout, "Warning: data read does not matchn in partition offsets\n");
   
            // Iterate through chunk to find symbol and its offset, this marks the end of the partition
            for(c = chunks - 1; c >= 0; c--)
            {
              offset++;
              (*offs)[i*3+2] = partSz - offset - (*offs)[i*3+1];
              if( (*offs)[i*3+2] < 1 )
              {
                err = 2;
                lparts--;
                fprintf(stdout, "Warning: adjusted partitions to (%d) based on data\n", lparts);
                break;
              }
              else if( buffer[c] == sym )
              {
                foundSymbol = 1;
                break;
              }
            }
          } 
        }
        // Last partition
        else 
        {
          // Compute a page size offset that is not after the end of the previous partition
          prevOffs = (*offs)[(i-1)*3] + (*offs)[(i-1)*3+1] + (*offs)[(i-1)*3+2];
          chunkOffs = (long long int)floor((double)prevOffs/chunks);
     
          (*offs)[i*3] = chunks * chunkOffs;
          (*offs)[i*3+1] = prevOffs - (*offs)[i*3];
          (*offs)[i*3+2] = sz - ((*offs)[i*3] + (*offs)[i*3+1]);
        }
      }
    }

    free(buffer);
  
    // Adjust occurred
    if( err == 2 )
      continue;
    // Error occurred
    else if( err != 0 )
    {
      lparts = ERROR;
      free(*offs);
      *offs = NULL;
      return ERROR;
    }
   
    // We are done 
    break;
  }

  *parts = lparts;  

  return 0;
}

