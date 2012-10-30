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

#ifndef EXTENDJOBCREATOR_H
#define EXTENDJOBCREATOR_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "ReadFile.h"
#include "DynamicProgrammingAligner.h"
#include "ScoredSeqCollectionBwt.h"
#include "AssemblyJobFactory.h"
#include "AssemblyJob.h"
#include "ScoredSeqNormalized.h"
#include "ScoredSeqWithCarriers.h"
#include "ScoredSeqFlip.h"
#include "ParamsMapping.h"
#include "AssemblerListener.h"
#include "ExtendCycle.h"
#include "AssemblyGraphNode.h"
#include "AssemblyJobGraph.h"
using namespace::std;

class ExtendJobCreator {

 public:
  // no constructors, destructors or fields - this is a class of static methods

  // the linked contig (should be the SSNormalized, I will maybe require that casting later)
  // the sense is the sense index for the collection in which the center contig (input for
  // the helper method) is supposed to show up in contig's two edge jobs.
  struct Linkage {
    ScoredSeq* contig;
    int senseIndex;
    float score;
  };

  // the method for public use (all are public for testing purposes)
  static vector<AssemblyJob*>* createJobs(set<ScoredSeqWithCarriers*>* mapSets, set<ScoredSeqNormalized*>* doNotExtendContigs,
					  AssemblyJobFactory* jobFactory, long maxContigLinkages,
					  set<ScoredSeqWithCarriers*>* garbageDump);

  // the meat of the process - guarantees that all nodes will be included in the 
  static vector<AssemblyJob*>* makeJobsFromGraph(AssemblyGraphNode** sortedNodeArray, long numNodes,
						 set<ScoredSeqNormalized*>* doNotExtendContigs,
						 AssemblyJobFactory* jobFactory);

  // returned array is the same size as the specified contigArray size, with contigs in the same order.
  static AssemblyGraphNode** makeGraphStructure(ScoredSeqWithCarriers** contigArray, long numContigs,
						long maxContigLinkages, set<ScoredSeqWithCarriers*>* allReadsWithCarriers);

  // a helper function; sorts the nodeArray in place
  static void sortGraphNodes(AssemblyGraphNode** nodeArray, long numNodes);

  // returns the linked contig within a flip, the orientation of which is the orientation that the
  // contig would be in to be linked in the indicated direction.  if there is no linkage, NULL is returned.
  // the nested seq is the SSNormalized contig.
  static Linkage createExtendJobHelper(set<ScoredSeqFlip*>* subJob, ScoredSeqWithCarriers* contigWc, int senseIndex,
				    long maxContigLinkages,
				    set<ScoredSeqWithCarriers*>* allReadsWithCarriers);

  // REQUIRES: pe carries SSWC objects
  // REQUIRES: numReadTargetMatches > 0
  static float measureLinkageWeight(ScoredSeqWithCarriers* pe, ScoredSeq* contig,
				    int senseIndex, long numReadTargetMatches,
				    map<ScoredSeq*,float>* contigsToFlip,
				    float totalLinkCount);

  // helper; makes sure that the same contig is not added to the same job
  // twice in two orientations.  modifies the components of the contig array
  // and their component reads.  deletes reads that are no longer accessible
  // from the contig array.
  //static void removeDsLinks(ScoredSeqWithCarriers** contigArray, long numContigs);
};

#endif
