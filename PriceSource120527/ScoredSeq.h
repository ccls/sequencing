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

/*
This is the representation for sequences (reads and/or contigs) in PRICE.  Abstractly, a ScoredSeq
has three properties: 1) a series of nucleotide identities, 2) a series of confidence scores for 
the indicated nucleotide identities, one per nucleotide, and 3) a series of confidence scores for th
adjacency of two nucleotides (essentially a score for the existance of the indicated 3p-5p phosphodiester
bond).  Both of these scores are decimal numbers that are supposed to reflect the amount of data
in support of the nucleotide identity/phosphodiester bond, i.e. a nucleotide at a position covered by
5 sequences, with all of the sequences in agreement for that nucleotide's identity, would have a score
of 5.  The decimal nature of the score allows quality scores (such as those provided by .fastq files)
to be taken into account.  Scores represent the possibility of a nucleotide identity being mis-called.
So in the example above, if another read was added to support the nucleotide identity, but that read had
a 1% scored possibility of being incorrect, the resulting contig's score at that position would become
5.99.  Ambiguous nucleotides (N's) have a score of zero by definition.  Intermediate ambiguities are not
allowed (i.e. Y for pyrimidines).  All nucleotides are represented as upper-case DNA letters (A,T,C,G - U
is not recognized by a ScoredSeq).</li>

SPEC FIELDS:
seq: the sequence in nucleic acids (A,T,C,G, or N)
scoreList: the scores of the nucleotide identities
linkerList: the scores of the phosphodiester bond
 */

#ifndef SCOREDSEQ_H
#define SCOREDSEQ_H
# include <vector>
# include <string>
using std::string;
using std::vector;


class ScoredSeq {

 public:
  enum Nuc{ BASE=0, A=1, C=2, G=3, T=4, N=5 };

  static ScoredSeq* getScoredSeq(char* seq, float score, long size);
  static ScoredSeq* getScoredSeq(char* seq, float* scores, long size);
  static ScoredSeq* getScoredSeq(char* seq, float* scores, float* links, long size);
  static ScoredSeq* copyShallowSeq(ScoredSeq* seq, char sense);
  // uses the input values directly as the stored values
  static ScoredSeq* repExposedSeq(char* seq, float score, long size);
  static ScoredSeq* repExposedSeq(char* seq, float* scores, long size);
  static ScoredSeq* repExposedSeq(char* seq, float* scores, float* links, long size);
  virtual ~ScoredSeq();

  /* returns the nucleotide identity at the indicated position
     on the indicated strand.  The position is zero-indexed from the
     5p end of the sequence; i.e. for the nucleotide at position n
     on the '+' strand, its reverse complement would be obtained by
     getting the nucleotide at position (lengthOfSeq - n - 1) from the
     '-' strand.  Remember: numbers always ascend in the 5p->3p direction.
     RETURNS: nucleotide identity
     REQUIRES: sense is '+' or '-'
     THROWS: ??? if position >= length of this
     PARAMS: position: zero-indexed, nums ascend 5p->3p
             sense: '+' or '-'
  */
  virtual char nucAtPosition(long position, char sense) = 0;
  // optimized to not have to pass the sense char and not have to ask which sense it is
  virtual char nucAtPosPlus(long position) = 0;
  virtual char nucAtPosMinus(long position) = 0;

  /* returns the score for the nucleotide called at the
     indicated position on the indicated strand, numbered 5p->3p
     on that strand
     RETURNS: nucleotide confidence score
     REQUIRES: sense is '+' or '-'
     THROWS: ??? if position >= length of this
     PARAMS: position: zero-indexed, nums ascend 5p->3p
             sense: '+' or '-'
  */
  virtual float scoreAtPosition(long position, char sense) = 0;
  virtual float scoreAtPosPlus(long position) = 0;
  virtual float scoreAtPosMinus(long position) = 0;

  virtual float linkAfterPosition(long position, char sense) = 0; // must be <= nucleotide counts
  virtual float linkAfterPosPlus(long position) = 0; // must be <= nucleotide counts
  virtual float linkAfterPosMinus(long position) = 0; // must be <= nucleotide counts

  /* gets the full nucleotide sequence as a string, presented 5p->3p
     RETURNS: the nucelotide DNA sequence (composed of A,T,C,G,N)
     REQUIRES: sense is '+' or '-'
     PARAMS: sense: '+' or '-'
  */
  virtual char* getSeq(char sense) = 0;
  virtual float* getScores(char sense) = 0;
  virtual float* getLinks(char sense) = 0;

  virtual char* getSubseq(long position, long length, char sense) = 0;
  virtual float* getSubScores(long position, long length, char sense) = 0;
  virtual float* getSubLinks(long position, long length, char sense) = 0;

  /* has no effects on state; purely an optimization tool that can make
     retrieval of the sequence or scores using the above methods more efficient
     later in some implementations.
     EFFECTS: none
  */
  virtual void buffer() = 0;
  /* for a multi-nested seq, allows just the bottom level
   * to be buffered */
  virtual void bottomBuffer() = 0;

  /* has no effects on state; purely an optimization tool that can save
     memory in some implementations if called after sequence/score info
     is accessed.
     EFFECTS: none
  */
  virtual void unbuffer() = 0;

  /* RETURNS: length of the sequence in number of nucleotides */
  virtual long size() = 0;

  /* Identifies whether or not the sequence has a paired end.  Paired-ends
     are defined as sequences that derive from the same DNA template as one 
     another and derive from opposing strands.  they may overlap, but the 
     non-sequencing-primer-derived sequences will always be mutually downstream
     of one another (at worst, perfectly overlapping). */
  virtual bool hasPairedEnd() = 0;

  /* if there is a paired end, it can be provided as a ScoredSeq by this method.
     See hasPairedEnd() comments for a description of paired ends.
     THROWS: ??? if this does not have a paired-end.
  */
  virtual ScoredSeq * getPairedEnd() = 0;
  virtual ScoredSeq * getTempPairedEnd() = 0;

  /* returns a copy of the ScoredSeq with shallow scope and so can be destroyed without adverse effects */
  virtual ScoredSeq * shallowCopy() = 0;

  /* returns a copy of the ScoredSeq with shallow scope and also with opposite sense as the original */
  virtual ScoredSeq * flipCopy() = 0;

  /* for wrapper classes */
  virtual bool isNested() = 0;

  /* these methods are supported in all cases even if they are only meaningfully different from their
   * shallow counterparts (destructor and unbuffer) for nested classes */
  virtual void deepDelete() = 0;
  virtual void deepUnbuffer() = 0;

  static char* reverseComplement(char* seq, long length);
 private:
  static Nuc _revComp[6];

};

#endif
