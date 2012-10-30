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


/* this is for when the assembly job is for a series of linked contigs,
 * i.e. a scaffold or super-contig.  it is sort of a localized meta-assembly.
 * in the base case where one contig has not been linked to anything by 
 * overwhelming evidence, it combines the two edge assemblies for that contig.
 */

#ifndef ASSEMBLYJOBGRAPH_H
#define ASSEMBLYJOBGRAPH_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "ScoredSeqNormalized.h"
#include "Alignment.h"
#include "AssemblyJob.h"
#include "AssemblyJobFactory.h"
#include "AssemblyGraphNode.h"
#include "ScoredSeqCollectionBwt.h"
using namespace::std;

class AssemblyJobGraph: public AssemblyJob {

 public:

  static AssemblyJob* AssemblyJobGraphFactory(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode,
					      AssemblyJobFactory* jobFactory);

  // REQUIRES: all of the nodes must be concordant in terms of their senses (they will not be flipped over)
  // REQUIRES: by following plus linkages, the beginNode must eventually point to (or be) the endNode
  AssemblyJobGraph(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode,
		   AssemblyJobFactory* jobFactory);
  ~AssemblyJobGraph();

  /* see specs from parent class */
  void makeShallow(bool shallowDelete);
  void runJob(AssemblyJob::AssemblyThreadedness threadedness);
  bool wasRun();

  void allInputSeqs( set<ScoredSeq*>* );
  void remainingInputSeqs( set<ScoredSeq*>* );
  void discardedInputSeqs( set<ScoredSeq*>* );
  void newlyAssembledSeqs( set<ScoredSeq*>* );
  void allAssembledSeqs( set<ScoredSeq*>* );

  // class-specific (and really just here for ExtendJobCreator testing purposes)
  void OK();

 private:
  // set-up functions
  void constructorHelper(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode, AssemblyJobFactory* ajf);
  void createNormSets(AssemblyJob::AssemblyThreadedness threadedness);
  void createNormSetsHelper(long seqN, char sense, ScoredSeq** seqArray, long* numInstances, bool* inArray, long* indexArray);
  // delegation of running assembly jobs
  void collapseOneNode(long tipIndex, long edgeIndex, set<ScoredSeq*>* localProducts, AssemblyJob::AssemblyThreadedness threadedness);
  void runCenterAssemblyJob(long tipIndex, long edgeIndex, AssemblyJob::AssemblyThreadedness threadedness);
  void assembleOneEdge(long edgeIndex, AssemblyJob::AssemblyThreadedness threadedness);
  void insertFullDiscard(ScoredSeq* seq, bool delSeq);

  bool _wasRun;


  // for navigation through the arrays below
  long _numNodes;

  // for more permanent storage/return
  set<ScoredSeq*> _rawInputTotalPlus;
  set<ScoredSeq*> _rawInputTotalMinus;
  bool _isShallow;

  // temporary storage
  set<ScoredSeq*>* _rawInputSetsPlus;
  set<ScoredSeq*>* _rawInputSetsMinus;

  // these contain only the normalized sequences
  set<ScoredSeq*> _justNormPlus;
  set<ScoredSeq*> _justNormMinus;
  map<ScoredSeq*,long> _normToLivingCount;
  // these contain all of the sequences, but only the ones
  // that appear in multiple bins in the bottom set arrays are normalized
  set<ScoredSeq*>* _normInputSetsPlus;
  set<ScoredSeq*>* _normInputSetsMinus;


  // arrays of sets of scoredseqs: len = num nodes + 1
  set<ScoredSeq*>* _edgeAssembledPlus;
  set<ScoredSeq*>* _edgeAssembledMinus;
  // as product seqs are removed 
  set<ScoredSeq*> _fullNovelSet;
  set<ScoredSeq*> _fullRemainingSet;
  set<ScoredSeq*> _fullDiscardSet;

  // applies to the meta-assembly of two adjacent edges from above; len = num nodes
  // index n applies to the combination of _plus/minusSets[n] with _plus/minusSets[n+1]
  long* _minOverlaps;

  // applies to the edge sets; true if the contents need to be individually assembled before they
  // are included in a 2-edge meta-assembly
  bool* _needsToBeAssembled;

  AssemblyJobFactory* _jobFactory;

  bool _inputRemoved;

};

#endif
