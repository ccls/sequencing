/*
PRICE was written by Graham Ruby.
This file is part of PRICE.

PRICE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PRICE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PRICE.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef READFILE_CPP
#define READFILE_CPP

#include "ReadFile.h"
#include "AssemblyException.h"
#include "ReadFileFastaSingle.h"
#include "ReadFileFastqSingle.h"
#include "ReadFilePriceqSingle.h"
#include "ReadFileDual.h"
#include "ReadFilePaired.h"

ReadFile* ReadFile::makeBasicReadFile(string filename, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rf;
  if (fileType == FASTQFILE){ rf = new ReadFileFastqSingle(filename,ReadFileFastqSingle::fastqSystem,countFactor); }
  else if (fileType == PRICEQFILE){ rf = new ReadFilePriceqSingle(filename); }
  else if (fileType == FASTAFILE){ rf = new ReadFileFastaSingle(filename,countFactor); }
  else if (fileType == ILLUMINAFILE){ rf = new ReadFileFastqSingle(filename,ReadFileFastqSingle::illuminaSystem,countFactor); }
  else { throw AssemblyException::ArgError("read file is of an inappropriate type."); }
  return rf;
}
ReadFile* ReadFile::makeInvertedReadFile(string filename, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rf;
  if (fileType == FASTQFILE){ rf = new ReadFileFastqSingle(true,filename,ReadFileFastqSingle::fastqSystem,countFactor); }
  else if (fileType == PRICEQFILE){ rf = new ReadFilePriceqSingle(true,filename); }
  else if (fileType == FASTAFILE){ rf = new ReadFileFastaSingle(filename,countFactor,true); }
  else if (fileType == ILLUMINAFILE){ rf = new ReadFileFastqSingle(true,filename,ReadFileFastqSingle::illuminaSystem,countFactor); }
  else { throw AssemblyException::ArgError("read file is of an inappropriate type."); }
  return rf;
}


ReadFile* ReadFile::makePairedReadFile(string filename, long ampliconSize, float countFactor){
  ReadFile* rfh = makeBasicReadFile(filename, countFactor);
  ReadFile* rf = new ReadFilePaired(rfh, ampliconSize);
  delete rfh;
  return rf;
}
ReadFile* ReadFile::makePairedReadFile(string filenameA, string filenameB, long ampliconSize, float countFactor){
  ReadFile* rfhA = makeBasicReadFile(filenameA, countFactor);
  ReadFile* rfhB = makeBasicReadFile(filenameB, countFactor);
  ReadFile* rf = new ReadFileDual(rfhA, rfhB, ampliconSize);
  delete rfhA;
  delete rfhB;
  return rf;
}


ReadFile* ReadFile::makeMatePairFile(string filename, long ampliconSize, float countFactor){
  ReadFile* rfh = makeInvertedReadFile(filename, countFactor);
  ReadFile* rf = new ReadFilePaired(rfh, ampliconSize);
  delete rfh;
  return rf;
}
ReadFile* ReadFile::makeMatePairFile(string filenameA, string filenameB, long ampliconSize, float countFactor){
  ReadFile* rfhA = makeInvertedReadFile(filenameA, countFactor);
  ReadFile* rfhB = makeInvertedReadFile(filenameB, countFactor);
  ReadFile* rf = new ReadFileDual(rfhA, rfhB, ampliconSize);
  delete rfhA;
  delete rfhB;
  return rf;
}


ReadFile* ReadFile::makeFalsePairFile(string filename, long readSize, long ampliconSize, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rfA;
  ReadFile* rfB;
  if (fileType == FASTQFILE){
    rfA = new ReadFileFastqSingle(false, readSize, filename,ReadFileFastqSingle::fastqSystem, countFactor);
    rfB = new ReadFileFastqSingle(true, readSize, filename,ReadFileFastqSingle::fastqSystem, countFactor);
  } else if (fileType == FASTAFILE){
    rfA = new ReadFileFastaSingle(filename, countFactor, readSize, false);
    rfB = new ReadFileFastaSingle(filename, countFactor, readSize, true);
  } else if (fileType == ILLUMINAFILE){
    rfA = new ReadFileFastqSingle(false, readSize, filename,ReadFileFastqSingle::illuminaSystem, countFactor);
    rfB = new ReadFileFastqSingle(true, readSize, filename,ReadFileFastqSingle::illuminaSystem, countFactor);
  } else { throw AssemblyException::ArgError("read file is of an inappropriate type for a false-pair file."); }
  ReadFile* rf = new ReadFileDual(rfA, rfB, ampliconSize);
  delete rfA;
  delete rfB;
  return rf;
}


ReadFile::FileType ReadFile::getFileType(string filename){
  size_t lastDot = filename.rfind('.');
  if (lastDot == string::npos){ return NONEFILE; }
  else {
    string append = filename.substr(lastDot);
    if (append==".fa" or append==".fna" or append==".ffn" or append==".frn" or append==".fasta"){ return FASTAFILE; }
    else if (append == ".fq" or append == ".fastq"){ return FASTQFILE; }
    else if (append == ".pq" or append == ".priceq"){ return PRICEQFILE; }
    else if (append == ".txt"){
      size_t underscore = filename.rfind('_');
      if (underscore == string::npos){ return NONEFILE; }
      else {
	string uAppend = filename.substr(underscore);
	if (uAppend == "_sequence.txt"){ return ILLUMINAFILE; }
	else { return NONEFILE; }
      }
    }
    else { return NONEFILE; }
  }
}


ReadFile::~ReadFile(){}


ReadFile::StringInterpreter::StringInterpreter(char* rawSeq, long rawLen, long maxLen, set<char>* okToSkip, bool revComp){
  char* newSeqExtra = new char[maxLen+1];
  long n2 = 0;
  if (revComp){
    // remember: this IS the RC case
    long n = rawLen;
    while (n > 0 and n2 < maxLen){
      --n;
      switch(rawSeq[n]) {
      // allowed nucleotide characters
      case 'A': newSeqExtra[n2] = 'T'; break;
      case 'C': newSeqExtra[n2] = 'G'; break;
      case 'G': newSeqExtra[n2] = 'C'; break;
      case 'T': newSeqExtra[n2] = 'A'; break;
      case 'U': newSeqExtra[n2] = 'A'; break;
      case 'R': newSeqExtra[n2] = 'N'; break;
      case 'Y': newSeqExtra[n2] = 'N'; break;
      case 'S': newSeqExtra[n2] = 'N'; break;
      case 'W': newSeqExtra[n2] = 'N'; break;
      case 'K': newSeqExtra[n2] = 'N'; break;
      case 'M': newSeqExtra[n2] = 'N'; break;
      case 'B': newSeqExtra[n2] = 'N'; break;
      case 'D': newSeqExtra[n2] = 'N'; break;
      case 'H': newSeqExtra[n2] = 'N'; break;
      case 'V': newSeqExtra[n2] = 'N'; break;
      case 'N': newSeqExtra[n2] = 'N'; break;
      case 'a': newSeqExtra[n2] = 'T'; break;
      case 'c': newSeqExtra[n2] = 'G'; break;
      case 'g': newSeqExtra[n2] = 'C'; break;
      case 't': newSeqExtra[n2] = 'A'; break;
      case 'u': newSeqExtra[n2] = 'A'; break;
      case 'r': newSeqExtra[n2] = 'N'; break;
      case 'y': newSeqExtra[n2] = 'N'; break;
      case 's': newSeqExtra[n2] = 'N'; break;
      case 'w': newSeqExtra[n2] = 'N'; break;
      case 'k': newSeqExtra[n2] = 'N'; break;
      case 'm': newSeqExtra[n2] = 'N'; break;
      case 'b': newSeqExtra[n2] = 'N'; break;
      case 'd': newSeqExtra[n2] = 'N'; break;
      case 'h': newSeqExtra[n2] = 'N'; break;
      case 'v': newSeqExtra[n2] = 'N'; break;
      case 'n': newSeqExtra[n2] = 'N'; break;
      case '.': newSeqExtra[n2] = 'N'; break;
      default:
	if (okToSkip->count(rawSeq[n]) == 0){
	  throw AssemblyException::ArgError("nucleotide seq contains invalid char");
	} else { n2--; }
      // no default; any other letter is skipped
      }
      n2++;
    }
  } else {
    // remember: this is NOT the RC case
    long n = 0;
    while (n < rawLen and n2 < maxLen){
      switch(rawSeq[n]) {
      // allowed nucleotide characters
      case 'A': newSeqExtra[n2] = 'A'; break;
      case 'C': newSeqExtra[n2] = 'C'; break;
      case 'G': newSeqExtra[n2] = 'G'; break;
      case 'T': newSeqExtra[n2] = 'T'; break;
      case 'U': newSeqExtra[n2] = 'T'; break;
      case 'R': newSeqExtra[n2] = 'N'; break;
      case 'Y': newSeqExtra[n2] = 'N'; break;
      case 'S': newSeqExtra[n2] = 'N'; break;
      case 'W': newSeqExtra[n2] = 'N'; break;
      case 'K': newSeqExtra[n2] = 'N'; break;
      case 'M': newSeqExtra[n2] = 'N'; break;
      case 'B': newSeqExtra[n2] = 'N'; break;
      case 'D': newSeqExtra[n2] = 'N'; break;
      case 'H': newSeqExtra[n2] = 'N'; break;
      case 'V': newSeqExtra[n2] = 'N'; break;
      case 'N': newSeqExtra[n2] = 'N'; break;
      case 'a': newSeqExtra[n2] = 'A'; break;
      case 'c': newSeqExtra[n2] = 'C'; break;
      case 'g': newSeqExtra[n2] = 'G'; break;
      case 't': newSeqExtra[n2] = 'T'; break;
      case 'u': newSeqExtra[n2] = 'T'; break;
      case 'r': newSeqExtra[n2] = 'N'; break;
      case 'y': newSeqExtra[n2] = 'N'; break;
      case 's': newSeqExtra[n2] = 'N'; break;
      case 'w': newSeqExtra[n2] = 'N'; break;
      case 'k': newSeqExtra[n2] = 'N'; break;
      case 'm': newSeqExtra[n2] = 'N'; break;
      case 'b': newSeqExtra[n2] = 'N'; break;
      case 'd': newSeqExtra[n2] = 'N'; break;
      case 'h': newSeqExtra[n2] = 'N'; break;
      case 'v': newSeqExtra[n2] = 'N'; break;
      case 'n': newSeqExtra[n2] = 'N'; break;
      case '.': newSeqExtra[n2] = 'N'; break;
      default:
	if (okToSkip->count(rawSeq[n]) == 0){
	  throw AssemblyException::ArgError("nucleotide seq contains invalid char");
	} else { n2--; }
      // no default; any other letter is skipped
      }
      n++;
      n2++;
    }
  }

  // I am hard-coding this decision about how much memory saving is worth the
  // time saved by not copying
  if (float(n2) / float(maxLen) > 0.9){ _seqString = newSeqExtra; }
  else {
    _seqString = new char[n2+1];
    for (long n = 0; n < n2; ++n){ _seqString[n] = newSeqExtra[n]; }
    delete [] newSeqExtra;
  }
  _seqString[n2] = '\0';
  _seqLen = n2;
}
ReadFile::StringInterpreter::~StringInterpreter(){}


#endif
