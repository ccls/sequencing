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


/* this is a node on an assembly job graph.  it is centered around a
 * contig that is to be extended on either end (the two collections of
 * sequences are the sequences that are to be used for each of those
 * two jobs), but it can be linked to other such nodes to form a graph
 * that represents a series of linked, related jobs (and, in a way, also
 * represents a super-contig).
 */

#ifndef ASSEMBLYGRAPHNODE_H
#define ASSEMBLYGRAPHNODE_H

#include <set>
#include <map>
#include "ScoredSeqFlip.h"
#include "Alignment.h"
#include "AssemblyJob.h"
#include "AssemblyJobFactory.h"
using namespace::std;

class AssemblyGraphNode {

 public:
  //AssemblyGraphNode();
  // REQUIRES: SSFlip senses must be worked out so that the two seq sets are concordant in terms of their sense.
  // (would jive if combined into a single-stranded assembly job)
  AssemblyGraphNode(ScoredSeq* centerContig, set<ScoredSeqFlip*>* plusEdgeSeqs, set<ScoredSeqFlip*>* minusEdgeSeqs);
  ~AssemblyGraphNode();

  void visit();
  bool wasVisited();
  void unvisit();

  void setPlusLink(AssemblyGraphNode* plusLink, float linkScore, bool sameDir);
  void setMinusLink(AssemblyGraphNode* minusLink, float linkScore, bool sameDir);
  AssemblyGraphNode* getPlusLink();
  AssemblyGraphNode* getMinusLink();
  float getPlusLinkScore();
  float getMinusLinkScore();
  bool plusLinkAgrees();
  bool minusLinkAgrees();
  void setPlusAgreement(bool sameDir);
  void setMinusAgreement(bool sameDir);

  ScoredSeq* getCenterContig();
  bool centerContigIsPlus();
  void flipOrientation();

  // does not return the original SSFlips, but the ones it does return are equivalent to the input
  // also, the flip carriers can (and should) be deleted in the scope where they were obtained
  void getPlusEdgeSeqs(set<ScoredSeqFlip*>* edgeSeqs);
  void getMinusEdgeSeqs(set<ScoredSeqFlip*>* edgeSeqs);

 private:
  void flipOrientationHelper(set<ScoredSeqFlip*>* setToFlip, set<ScoredSeqFlip*>* garbageBin);

  bool _wasVisited;
  ScoredSeq* _centerContig;
  bool _centerIsPlus;
  set<ScoredSeqFlip*> _plusEdgeSeqs;
  set<ScoredSeqFlip*> _minusEdgeSeqs;

  AssemblyGraphNode* _plusLink;
  AssemblyGraphNode* _minusLink;
  float _plusLinkScore;
  float _minusLinkScore;
  bool _plusLinkAgrees;
  bool _minusLinkAgrees;

};

#endif
