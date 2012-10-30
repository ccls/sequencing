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

#ifndef EXTENDJOBCREATOR_CPP
#define EXTENDJOBCREATOR_CPP

#include "ExtendJobCreator.h"
#include "ExtendJobMapper.h"
#include "Read.h"
#include "ScoredSeqSubseq.h"
#include "ScoredSeqNested.h"
#include "ScoredSeqMonoScore.h"
#include "AssemblyJobNullCopy.h"
#include <typeinfo>
#include <omp.h>
using namespace::std;



vector<AssemblyJob*>* ExtendJobCreator::createJobs(set<ScoredSeqWithCarriers*>* mapSets, set<ScoredSeqNormalized*>* doNotExtendContigs,
						   AssemblyJobFactory* jobFactory, long maxContigLinkages,
						   set<ScoredSeqWithCarriers*>* garbageDump){

  ScoredSeqWithCarriers** contigArray = new ScoredSeqWithCarriers*[ mapSets->size() ];
  long numContigsTotal = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = mapSets->begin(); it != mapSets->end(); ++it){
    contigArray[numContigsTotal] = *it;
    numContigsTotal++;
  }
  // has been threaded
  // THIS WAS MADE OBSOLETE BY MOVING TO THE READ MAPPER
  //removeDsLinks(contigArray, numContigsTotal);
  // could be threaded???
  AssemblyGraphNode** agnArray = ExtendJobCreator::makeGraphStructure(contigArray, numContigsTotal, maxContigLinkages, garbageDump);
  // could be threaded???
  ExtendJobCreator::sortGraphNodes(agnArray, numContigsTotal);
  // cannot be usefully threaded
  vector<AssemblyJob*>* ajSet = ExtendJobCreator::makeJobsFromGraph(agnArray, numContigsTotal, doNotExtendContigs, jobFactory);
  delete [] contigArray;
  for (long n = 0; n < numContigsTotal; ++n){ delete agnArray[n]; }
  delete [] agnArray;

  return ajSet;
}


vector<AssemblyJob*>* ExtendJobCreator::makeJobsFromGraph(AssemblyGraphNode** sortedNodeArray, long numNodes,
							  set<ScoredSeqNormalized*>* doNotExtendContigs,
							  AssemblyJobFactory* jobFactory){

  // i need to re-cast, otherwise this is just too annoying
  set<ScoredSeq*> dneContigs;
  set<ScoredSeq*>::iterator dcIt = dneContigs.begin();
  for (set<ScoredSeqNormalized*>::iterator it = doNotExtendContigs->begin(); it != doNotExtendContigs->end(); ++it){
    dcIt = dneContigs.insert(dcIt, *it);
  }

  vector<AssemblyJob*>* jobsToBeRun = new vector<AssemblyJob*>;

  // this loop cannot be threaded because the "visited" status is updated in the loop,
  // and ignoring that would allow for contigs to be inappropriately/redundantly added
  // to multiple jobs
  for (long nodeN = 0; nodeN < numNodes; ++nodeN){
    AssemblyGraphNode* beginNode = sortedNodeArray[nodeN];

    // do not seed a path with this contig if it was visited or if it is in the do-not-extend bin
    if ( beginNode->wasVisited() ){
      // nothing happens
    } else if ( dneContigs.count( beginNode->getCenterContig() ) == 1){
      beginNode->visit();
      // create a null job
      set<ScoredSeq*> justTheContig;
      justTheContig.insert( beginNode->getCenterContig() );
      jobsToBeRun->push_back( new AssemblyJobNullCopy( &justTheContig ) );
    } else {
      beginNode->visit();
      AssemblyGraphNode* endNode = beginNode;
      bool pushBegin = true;
      bool pushEnd = true;
      while(pushBegin or pushEnd){
	if (pushBegin){
	  AssemblyGraphNode* nextNode = beginNode->getMinusLink();
	  if (nextNode == NULL or nextNode->wasVisited() or dneContigs.count( nextNode->getCenterContig() ) == 1){
	    pushBegin = false;
	  } else {
	    nextNode->visit();
	    if (nextNode->getPlusLink() != beginNode){ nextNode->flipOrientation(); }
	    if (nextNode->getPlusLink() != beginNode){
	      throw AssemblyException::LogicError("EJC:mJFG begin links should reciprocate at this point.");
	    }
	    beginNode = nextNode;
	  }
	}
	if (pushEnd){
	  AssemblyGraphNode* nextNode = endNode->getPlusLink();
	  if (nextNode == NULL or nextNode->wasVisited() or dneContigs.count( nextNode->getCenterContig() ) == 1){
	    pushEnd = false;
	  } else {
	    nextNode->visit();
	    if (nextNode->getMinusLink() != endNode){ nextNode->flipOrientation(); }
	    if (nextNode->getMinusLink() != endNode){
	      throw AssemblyException::LogicError("EJC:mJFG begin links should reciprocate at this point.");
	    }
	    endNode = nextNode;
	  }
	}
      }
      // this is just a temp fill-in since I don't have the graph job implemented yet
      // but I want to test to make sure that the correct number of jobs is created
      jobsToBeRun->push_back( AssemblyJobGraph::AssemblyJobGraphFactory(beginNode, endNode, jobFactory) );
    }
  }

  // delete the unvisited nodes... or delete all of the nodes???
  /*
  for (long nodeN = 0; nodeN < numNodes; ++nodeN){
    if (! sortedNodeArray[nodeN]->wasVisited()){ delete sortedNodeArray[nodeN]; }
  }
  */
  return jobsToBeRun;
}



AssemblyGraphNode** ExtendJobCreator::makeGraphStructure(ScoredSeqWithCarriers** contigArray, long numContigs,
							 long maxContigLinkages, set<ScoredSeqWithCarriers*>* allReadsWithCarriers){

  map<ScoredSeq*,AssemblyGraphNode*> contigToNode;
  AssemblyGraphNode** nodeArray = new AssemblyGraphNode*[ numContigs+1 ];
  Linkage* plusLinkArray = new Linkage[ numContigs+1 ];
  Linkage* minusLinkArray = new Linkage[ numContigs+1 ];

  // make nodes for ALL contigs, even doNotExtend contigs
  int plusIndex = ExtendCycle::plusIndex();
  int minusIndex = ExtendCycle::minusIndex();

  // these data structures avoid thread conflicts, but must be resolved after the threaded loop
  set<ScoredSeqWithCarriers*>* allReadsWcForThreads = new set<ScoredSeqWithCarriers*>[numContigs+1];
  ScoredSeq** contigNormArray = new ScoredSeq*[numContigs+1];

  // this loop calls methods that sort through all of the paired-end reads that
  // link contigs and organize that data; it is threaded by contig; it can be threaded
  // because it only collects information and creates new objects, it doesn't modify any
  // existing objects (except the thread-safe sets from the array created above)

  #pragma omp parallel for schedule(dynamic)
  for (long cN = 0; cN < numContigs; ++cN){
    ScoredSeqWithCarriers* contigWc = contigArray[cN];
    ScoredSeqNormalized* contigNorm = dynamic_cast<ScoredSeqNormalized*>( contigWc->getNested() );
    set<ScoredSeqFlip*> plusEdgeSeqs;
    set<ScoredSeqFlip*> minusEdgeSeqs;

    // make the node
    Linkage plusLink = createExtendJobHelper(&plusEdgeSeqs, contigWc, plusIndex, maxContigLinkages, &allReadsWcForThreads[cN]);
    Linkage minusLink = createExtendJobHelper(&minusEdgeSeqs, contigWc, minusIndex, maxContigLinkages, &allReadsWcForThreads[cN]);
    AssemblyGraphNode* contigNode = new AssemblyGraphNode(contigNorm, &plusEdgeSeqs, &minusEdgeSeqs);
    // now that they have been deposited, i can delete the flips
    for (set<ScoredSeqFlip*>::iterator flipIt = plusEdgeSeqs.begin(); flipIt != plusEdgeSeqs.end(); ++flipIt){ delete *flipIt; }
    for (set<ScoredSeqFlip*>::iterator flipIt = minusEdgeSeqs.begin(); flipIt != minusEdgeSeqs.end(); ++flipIt){ delete *flipIt; }

    // prepare for linkage later
    contigNormArray[cN] = contigNorm;
    nodeArray[cN] = contigNode;
    plusLinkArray[cN] = plusLink;
    minusLinkArray[cN] = minusLink;
  }

  // not-threaded recovery from the steps above
  for (long cN = 0; cN < numContigs; ++cN){
    contigToNode.insert( pair<ScoredSeq*,AssemblyGraphNode*>(contigNormArray[cN], nodeArray[cN]) );
    allReadsWithCarriers->insert(allReadsWcForThreads[cN].begin(), allReadsWcForThreads[cN].end());
  }
  delete [] allReadsWcForThreads;
  delete [] contigNormArray;

  // now create the links
  for (long contigIndex = 0; contigIndex < numContigs; ++contigIndex){
    AssemblyGraphNode* contigNode = nodeArray[contigIndex];

    // the plus link
    {
      Linkage link = plusLinkArray[contigIndex];
      // do not create the linkage if this node already has a superior linkage
      if (link.contig != NULL and (contigNode->getPlusLink() == NULL or contigNode->getPlusLinkScore() < link.score) ){
	AssemblyGraphNode* linkNode = contigToNode.find( link.contig )->second;
	if (link.senseIndex == plusIndex){
	  // do not create the linkage if the link already has a superior linkage
	  if (linkNode->getPlusLink() == NULL or linkNode->getPlusLinkScore() < link.score){
	    contigNode->setPlusLink(linkNode, link.score, true);
	    linkNode->setPlusLink(contigNode, link.score, true);
	  }
	} else {
  	  // do not create the linkage if the link already has a superior linkage
	  if (linkNode->getMinusLink() == NULL or linkNode->getMinusLinkScore() < link.score){
	    contigNode->setPlusLink(linkNode, link.score, false);
	    linkNode->setMinusLink(contigNode, link.score, false);
	  }
	}
      }
    }

    // the minus link
    {
      Linkage link = minusLinkArray[contigIndex];
      // do not create the linkage if this node already has a superior linkage
      if (link.contig != NULL and (contigNode->getMinusLink() == NULL or contigNode->getMinusLinkScore() < link.score) ){
	AssemblyGraphNode* linkNode = contigToNode.find( link.contig )->second;
	if (link.senseIndex == plusIndex){
	  // do not create the linkage if the link already has a superior linkage
	  if (linkNode->getPlusLink() == NULL or linkNode->getPlusLinkScore() < link.score){
	    contigNode->setMinusLink(linkNode, link.score, false);
	    linkNode->setPlusLink(contigNode, link.score, false);
	  }
	} else {
  	  // do not create the linkage if the link already has a superior linkage
	  if (linkNode->getMinusLink() == NULL or linkNode->getMinusLinkScore() < link.score){
	    contigNode->setMinusLink(linkNode, link.score, true);
	    linkNode->setMinusLink(contigNode, link.score, true);
	  }
	}
      }
    }
  }

  // resolve best-edge conflicts
  for (long contigIndex = 0; contigIndex < numContigs; ++contigIndex){
    AssemblyGraphNode* contigNode = nodeArray[contigIndex];
    AssemblyGraphNode* bestPlusLink = contigNode->getPlusLink();

    if (bestPlusLink != NULL){
      if (contigNode->plusLinkAgrees()){
      //Linkage link = plusLinkArray[contigIndex];
      //if (link.senseIndex == plusIndex){
	if (bestPlusLink->getPlusLink() != contigNode){
	  AssemblyGraphNode* bestRecip = bestPlusLink->getPlusLink();
	  if (bestRecip == NULL){ bestPlusLink->setPlusLink(contigNode, contigNode->getPlusLinkScore(), true); }
	  else { contigNode->setPlusLink( NULL, 0, true); }
	}
      } else {
	if (bestPlusLink->getMinusLink() != contigNode){
	  AssemblyGraphNode* bestRecip = bestPlusLink->getMinusLink();
	  if (bestRecip == NULL){ bestPlusLink->setMinusLink(contigNode, contigNode->getPlusLinkScore(), false); }
	  else { contigNode->setPlusLink( NULL, 0, true); }
	}
      }
    }

    AssemblyGraphNode* bestMinusLink = contigNode->getMinusLink();
    if (bestMinusLink != NULL){
      if (contigNode->minusLinkAgrees()){
	//Linkage link = minusLinkArray[contigIndex];
	//if (link.senseIndex == plusIndex){
	  if (bestMinusLink->getMinusLink() != contigNode){
	    AssemblyGraphNode* bestRecip = bestMinusLink->getMinusLink();
	    if (bestRecip == NULL){ bestMinusLink->setMinusLink(contigNode, contigNode->getMinusLinkScore(), true); }
	    else { contigNode->setMinusLink( NULL, 0, true ); }
	  }
      } else {
	  if (bestMinusLink->getPlusLink() != contigNode){
	    AssemblyGraphNode* bestRecip = bestMinusLink->getPlusLink();
	    if (bestRecip == NULL){ bestMinusLink->setPlusLink(contigNode, contigNode->getMinusLinkScore(), false); }
	    else { contigNode->setMinusLink( NULL, 0, true ); }
	}
      }
    }
  }

  delete [] plusLinkArray;
  delete [] minusLinkArray;
  return nodeArray;
}






void ExtendJobCreator::sortGraphNodes(AssemblyGraphNode** nodeArray, long numNodes){
  if (numNodes > 1){
    // this is going to make sortNodes O( n * log(m) ) where n in the num of nodes but m is the num
    // of unique scores, which is smaller.
    map<long,long> nodeLenToCount;
    // figure out how many instances of each score there are (many scores will be 1 or 0)
    for (long n = 0; n < numNodes; ++n){
      long nodeLen = nodeArray[n]->getCenterContig()->size();
      map<long,long>::iterator lengthBin = nodeLenToCount.find( nodeLen );
      if (lengthBin == nodeLenToCount.end()){ nodeLenToCount.insert( pair<long,long>(nodeLen,1) ); }
      else { lengthBin->second++; }
    }

    // figure out what index each score will start at and set up an array of counts for the score bin
    long numLengths = nodeLenToCount.size();
    long rankToStart[numLengths];
    long rankToCount[numLengths];
    // same keys as nodeToLenCount (i.e. sorted the same)
    map<long,long> nodeLenToRank;
    map<long,long>::iterator nltrIt = nodeLenToRank.begin();
    long rank = numLengths;
    long start = numNodes;
    for (map<long,long>::iterator it = nodeLenToCount.begin(); it != nodeLenToCount.end(); ++it){
      rank--;
      start -= it->second;
      rankToStart[rank] = start;
      rankToCount[rank] = 0;
      nltrIt = nodeLenToRank.insert( nltrIt, pair<long,long>( it->first, rank ) );
    }

    // now place the nodes appropriately in a new array
    AssemblyGraphNode** sortedArray = new AssemblyGraphNode*[numNodes];
    for (long n = 0; n < numNodes; ++n){
      long rank = nodeLenToRank.find( nodeArray[n]->getCenterContig()->size() )->second;
      sortedArray[rankToStart[rank] + rankToCount[rank]] = nodeArray[n];
      rankToCount[rank]++;
    }

    // now replace the contents of the old array with the sorted contents
    for (long n = 0; n < numNodes; ++n){ nodeArray[n] = sortedArray[n]; }
    delete [] sortedArray;
  }
}



// this function makes sure that the linked contigs are not the contig itself (prevents
// the contig from being introduced in the RC orientation or multiple times).
float ExtendJobCreator::measureLinkageWeight(ScoredSeqWithCarriers* pe, ScoredSeq* contig,
					     int senseIndex, long numReadTargetMatches,
					     map<ScoredSeq*,float>* contigsToFlip,
					     float totalLinkCount){

  if (numReadTargetMatches < 1){ throw AssemblyException::ArgError("EJC::mLW nRTM cannot be < 1"); }
  set<ScoredSeq*> replacementContigs;
  pe->getCarriedSeqs(&replacementContigs, senseIndex);
  if (replacementContigs.size() > 0){

    // normalize the linkage count to the product of the number of places that the read and its paired-end map
    float countsPerLink = float(1.0) / (float(replacementContigs.size()) * float(numReadTargetMatches));

    for (set<ScoredSeq*>::iterator linkIt = replacementContigs.begin(); linkIt != replacementContigs.end(); ++linkIt){
      ScoredSeq* localNested = dynamic_cast<ScoredSeqNested*>(*linkIt)->getNested();
      if ( localNested != contig ){
	map<ScoredSeq*,float>::iterator seqToCountIt = contigsToFlip->find( localNested );
	if (seqToCountIt == contigsToFlip->end()){
	  contigsToFlip->insert( pair<ScoredSeq*,float>(localNested,countsPerLink) );
	} else {
	  seqToCountIt->second += countsPerLink;
	}
	totalLinkCount += countsPerLink;
      }
    }
  }
  return totalLinkCount;
}

// MODIFIES allReadsWithCarriers
// MODIFIES subJob
// MODIFIES NOTHING ELSE!!!!!
ExtendJobCreator::Linkage ExtendJobCreator::createExtendJobHelper(set<ScoredSeqFlip*>* subJob, ScoredSeqWithCarriers* contigWc,
								  int senseIndex, long maxContigLinkages,
								  set<ScoredSeqWithCarriers*>* allReadsWithCarriers){

  ScoredSeqNormalized* contig = dynamic_cast<ScoredSeqNormalized*>( contigWc->getNested() );
  set<ScoredSeq*> matches;
  contigWc->getCarriedSeqs( &matches, senseIndex );

  // the index of the sense that 'plus' input seqs should be given in their flip objects
  int plusIndex = ExtendCycle::plusIndex();
  int minusIndex = ExtendCycle::minusIndex();
  char peSense; char peAntis;
  int pePlusIndex; int peMinusIndex;
  if (senseIndex == plusIndex){
    peSense = '-'; peAntis = '+';
    pePlusIndex = minusIndex; peMinusIndex = plusIndex;
  } else {
    peSense = '+'; peAntis = '-';
    pePlusIndex = plusIndex; peMinusIndex = minusIndex;
  }

  // calculate the max num linkages that a read may have in order for its paired end to be included;
  // it is the inverse of the average of the inverses of the match count numbers.
  long maxReadMatchNum = 0;
  if (matches.size() > 0){
    float sumOfCounts = 0;
    for (set<ScoredSeq*>::iterator mrIt = matches.begin(); mrIt != matches.end(); ++mrIt){
      ScoredSeqWithCarriers* read = dynamic_cast<ScoredSeqWithCarriers*>(*mrIt);
      sumOfCounts += float(1.0) / float(read->numCarriedSeqs(plusIndex) + read->numCarriedSeqs(minusIndex));
    }
    maxReadMatchNum = long(float(matches.size()) * (float(maxContigLinkages) / float(sumOfCounts)));
  }


  // contigs must be added non-redundantly, and the evidence supporting their linkage
  // must be kept track of so that a graph can be built for the best-supported linkages
  map<ScoredSeq*,float> contigsToFlipPlus;
  map<ScoredSeq*,float> contigsToFlipMinus;
  // this will include linked AND unlinked PE reads, but none that map to self
  float totalLinkCount = 0.0;

  // the reads are sorted for efficient set addition, but not their paired ends
  set<ScoredSeqWithCarriers*> readsToAdd;
  set<ScoredSeqWithCarriers*>::iterator rtaIt = readsToAdd.begin();
  set<ScoredSeqWithCarriers*> pairsToAdd;

  for (set<ScoredSeq*>::iterator mrIt = matches.begin(); mrIt != matches.end(); ++mrIt){
    ScoredSeqWithCarriers* read = dynamic_cast<ScoredSeqWithCarriers*>(*mrIt);
    ScoredSeqWithCarriers* pe = dynamic_cast<ScoredSeqWithCarriers*>( read->getPairedEnd() );
    rtaIt = readsToAdd.insert(rtaIt, read);
    pairsToAdd.insert( pe );

    //if (! read->hasCarriedSeq(contigWc,senseIndex) ){ throw AssemblyException::LogicError("in EJC, link isn't reciprocal"); }

    // i need to decide if the read will even be used based on the number of different contigs
    // to which it maps
    long numReadTargetMatches = read->numCarriedSeqs(plusIndex) + read->numCarriedSeqs(minusIndex);
    if (numReadTargetMatches <= maxReadMatchNum){

      // add pair-mapping contigs if there are not too many
      // check the sum of the mapped contigs is within the number of allowed contig replacements
      bool peReplaced = false;
      long numPeTargetMatches = pe->numCarriedSeqs(plusIndex) + pe->numCarriedSeqs(minusIndex);
      if (numPeTargetMatches <= maxContigLinkages){
	// plus strand matches
	if ( pe->numCarriedSeqs(plusIndex) > 0 ){
	  peReplaced = true;
	  totalLinkCount = measureLinkageWeight(pe, contig, plusIndex, numReadTargetMatches,
						&contigsToFlipPlus, totalLinkCount);
	}
	// minus strand matches
	if ( pe->numCarriedSeqs(minusIndex) > 0 ){
	  peReplaced = true;
	  totalLinkCount = measureLinkageWeight(pe, contig, minusIndex, numReadTargetMatches,
						&contigsToFlipMinus, totalLinkCount);
	}
      }
      // PE was retained
      if (! peReplaced){
	ScoredSeqNormalized* innerPe = dynamic_cast<ScoredSeqNormalized*>( pe->getNested() );
	subJob->insert( ScoredSeqFlip::getFlip( innerPe, peSense ) );
	totalLinkCount += float(1.0) / float(numReadTargetMatches);
      }
    }
  }
  // now add the reads and pairs to the deletion scheduler
  allReadsWithCarriers->insert( readsToAdd.begin(), readsToAdd.end() ); // schedules deletion of WithCarriers nest
  allReadsWithCarriers->insert( pairsToAdd.begin(), pairsToAdd.end() );

  // insert the contig itself
  subJob->insert( ScoredSeqFlip::getFlip(contig,'+') );


  // these will be updated as the linked contigs compete
  Linkage link = { NULL, 0, 0 };
  ScoredSeq* maxLinkSeq = 0;
  //map<ScoredSeq*,float>* maxLinkMap = NULL;
  float maxLinkCount = totalLinkCount / float(2);

  // make sure that the same contig is not being included in both senses,
  // EVEN IF it is because of two independent reads
  map<ScoredSeq*,float>::iterator plusIt = contigsToFlipPlus.begin();
  map<ScoredSeq*,float>::iterator minusIt = contigsToFlipMinus.begin();
  map<ScoredSeq*,float>::iterator plusEnd = contigsToFlipPlus.end();
  map<ScoredSeq*,float>::iterator minusEnd = contigsToFlipMinus.end();
  // these seqs should not be added to the final contig set; the iterator is for fast addition
  set<ScoredSeq*> removeThese;
  set<ScoredSeq*>::iterator remIt = removeThese.begin();

  // iterate through the two maps, adding seqs to be removed to another set
  while (minusIt != minusEnd){
    while ( plusIt != plusEnd and plusIt->first < minusIt->first ){ ++plusIt; }
    if (plusIt == plusEnd){ minusIt = minusEnd; }
    else if (plusIt->first == minusIt->first){
      remIt = removeThese.insert(remIt, plusIt->first);
      ++plusIt;
      ++minusIt;
    } else { ++minusIt; }
  }


  // insert the now non-redundant contigs that don't show up in both orientations
  set<ScoredSeq*>::iterator removeIt = removeThese.begin();
  set<ScoredSeq*>::iterator removeEnd = removeThese.end();
  for (map<ScoredSeq*,float>::iterator cIt = contigsToFlipPlus.begin(); cIt != contigsToFlipPlus.end(); ++cIt){
    // makes sure that contigs aren't added in both senses
    if (removeIt==removeEnd or *removeIt != cIt->first){
      subJob->insert( ScoredSeqFlip::getFlip( cIt->first, peSense ) );
      if (cIt->second >= maxLinkCount){
	link.contig = cIt->first; link.score = cIt->second;
	if (pePlusIndex == senseIndex){ link.senseIndex = minusIndex; }
	else { link.senseIndex = plusIndex; }
	maxLinkSeq = cIt->first;
	maxLinkCount = cIt->second;
      }
    } else { ++removeIt; }
  }

  // reset for minus contig addition
  removeIt = removeThese.begin();
  removeEnd = removeThese.end();
  for (map<ScoredSeq*,float>::iterator cIt = contigsToFlipMinus.begin(); cIt != contigsToFlipMinus.end(); ++cIt){
    // makes sure that contigs aren't added in both senses
    if (removeIt==removeEnd or *removeIt != cIt->first){
      subJob->insert( ScoredSeqFlip::getFlip( cIt->first, peAntis ) );
      if (cIt->second >= maxLinkCount){
	link.contig = cIt->first; link.score = cIt->second;
	if (peMinusIndex == senseIndex){ link.senseIndex = minusIndex; }
	else { link.senseIndex = plusIndex; }
	maxLinkSeq = cIt->first;
	maxLinkCount = cIt->second;
      }
    } else { ++removeIt; }
  }

  return link;
}



#endif
