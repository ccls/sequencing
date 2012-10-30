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


/* this is for when there are few/no input seqs to be assembled, so
 * the input should just be returned.
 */


#ifndef ASSEMBLYJOBSSDEBRUIJN_H
#define ASSEMBLYJOBSSDEBRUIJN_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "AssemblyJob.h"
#include "ParamsDeBruijn.h"
using namespace::std;

class AssemblyJobSsDeBruijn: public AssemblyJob {

 public:
  AssemblyJobSsDeBruijn();
  AssemblyJobSsDeBruijn(set<ScoredSeq*>* seqs, ParamsDeBruijn* paramsDeBruijn);
  ~AssemblyJobSsDeBruijn();

  /* see specs from parent class */
  void makeShallow(bool shallowDelete);
  void runJob(AssemblyJob::AssemblyThreadedness threadedness);
  bool wasRun();

  void allInputSeqs( set<ScoredSeq*>* );
  void remainingInputSeqs( set<ScoredSeq*>* );
  void discardedInputSeqs( set<ScoredSeq*>* );
  void newlyAssembledSeqs( set<ScoredSeq*>* );
  void allAssembledSeqs( set<ScoredSeq*>* );

  void OK();
  void testTreeConstruction(set<ScoredSeq*>* input, int kmerSize, set<ScoredSeq*>* foundKmers);

 private:
  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _novelSeqs;
  set<ScoredSeq*> _remainingSeqs;
  set<ScoredSeq*> _usedSeqs;
  bool _wasRun;
  int _kmerSize;
  ParamsDeBruijn* _paramsDeBruijn;

  // the internal representation of the de Bruijn graph is these three classes
  class SeqNode{
  public:
    virtual ~SeqNode();
    virtual void deepDelete() = 0;
    // ---for finding the sequence---
    // required true in order to add/get a branch
    virtual bool isBranchable() = 0;
    // returns NULL pointer (0) if there is no node
    virtual SeqNode* getBranch(char nuc) = 0;
    virtual void addBranch(SeqNode* sn) = 0;
    // ---for reconstructing the sequence---
    // returns NULL pointer (0) if there is no parent
    virtual SeqNode* getParent() = 0;
    // returns the null char '\0' if the top level was reached
    virtual char get3pNuc() = 0;
    // DOES NOT REFERENCE DATA FIELDS; JUST USES PUBLIC METHODS
    string getKmer();
    void getKmer(char* seqToFillIn, long startCoord);
    char get5pNuc();
  };

  // this is just for construction of the graph
  class NucNode : public SeqNode {
  public:
    NucNode();
    NucNode(SeqNode* parent, char nuc);
    ~NucNode();
    void deepDelete();
    // inherited
    // this class is branchable
    bool isBranchable();
    SeqNode* getBranch(char nuc);
    void addBranch(SeqNode* sn);
    SeqNode* getParent();
    char get3pNuc();
  private:
    int getNucIndex(char nuc);
    SeqNode* _branches[4];
    SeqNode* _parent;
    char _nuc;
  };

  // these are nodes on the actual De Bruijn graph
  class DeBruijnNode : public SeqNode {
  public:
    DeBruijnNode();
    DeBruijnNode(SeqNode* parent, char nuc);
    ~DeBruijnNode();
    void deepDelete();
    // inherited
    SeqNode* getBranch(char nuc);
    // this class is NOT branchable
    bool isBranchable();
    void addBranch(SeqNode* sn);
    SeqNode* getParent();
    char get3pNuc();
    // class-specific
    void add5pLink(DeBruijnNode* node, float linkScore);
    void add3pLink(DeBruijnNode* node, float linkScore);
    DeBruijnNode* get5pLink(char nuc);
    DeBruijnNode* get3pLink(char nuc);
    float get5pLinkScore(char nuc);
    float get3pLinkScore(char nuc);
    // properties of the kmer overall
    float getKmerScore();
    void addKmerScore(float score);
    // for keeping track during assembly to avoid circular contigs
    void visit();
    bool wasVisited();
    void setBest5pLink(DeBruijnNode* link);
    void setBest3pLink(DeBruijnNode* link);
    DeBruijnNode* best5pLink();
    DeBruijnNode* best3pLink();
  private:
    int getNucIndex(char nuc);
    SeqNode* _parent;
    char _nuc;
    DeBruijnNode* _5pLinks[4];
    DeBruijnNode* _3pLinks[4];
    float _5pLinkScores[4];
    float _3pLinkScores[4];
    float _kmerScore;
    bool _wasVisited;
    DeBruijnNode* _best5pLink;
    DeBruijnNode* _best3pLink;
  };

  // this is for use during the actual assembly step
  class DeBruijnGraph{
  public:
    DeBruijnGraph();
    DeBruijnGraph(SeqNode* topNode, long nodeCount, int kmerSize);
    ~DeBruijnGraph();
    void sortNodes();

    class Iterator{
    public:
      Iterator();
      Iterator(DeBruijnGraph* source, long nodeNum);
      Iterator(const Iterator& it);
      Iterator operator++();
      bool operator==(const Iterator& it);
      bool operator!=(const Iterator& it);
      DeBruijnNode* operator*();
    private:
      long _nodeNum;
      DeBruijnGraph* _source;
    };
    Iterator begin();
    Iterator end();

  protected:
    friend class Iterator;
    SeqNode* _topNode;
    DeBruijnNode** _nodeArray;
    long _nodeCount;
    int _kmerSize;
    long nodeArrayHelper(SeqNode* node, int stepsToBottom, long currentIndex);
  };


  // returns the top node of the tree and the number of bottom nodes
  pair<NucNode*,long> makeGraphTree(set<ScoredSeq*>* seqs, int kmerSize);
  void findBestLinks(DeBruijnGraph* graph);
  void makeSeqsFromGraph(DeBruijnGraph* graph, set<ScoredSeq*>* productSeqs);

  // requires that the start node will eventually meet the endNode when moving in the 5p->3p direction
  ScoredSeq* makeSeqFromPath(DeBruijnNode* node5p, DeBruijnNode* node3p, long numNodes);
  void testTreeConstructionHelper(set<ScoredSeq*>* foundKmers, SeqNode* node, int stepsToBottom);

  // these methods assist with threading
  long combineTwoTrees(SeqNode* topNodeA, SeqNode* topNodeB, int kmerSize);
  long combineTreesMidHelp(SeqNode* topNodeA, int kmerSize,
			   SeqNode* branchA, SeqNode* branchB, int remainingKmer);
  long combineTreesBottomHelp(SeqNode* topNodeA, int kmerSize,
			      DeBruijnNode* tipA, DeBruijnNode* tipB);
  pair<DeBruijnNode*,long> combineTreesAddNodeHelp(SeqNode* topNodeA, DeBruijnNode* bottomNodeB, int kmerSize);

};

#endif
